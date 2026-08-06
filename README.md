# Germline Variant Filtering Pipeline

A trio/duo germline variant-filtering pipeline for **clinical candidate triage** — it
annotates a VCF, then keeps rare, gene-panel variants that show damaging evidence
(missense pathogenicity, high CADD, splicing impact, or ClinVar) and labels them with
inheritance and recessive context for **downstream manual curation**.

> ⚠️ **Patient data never lives in this repository.** The pipeline runs on patient
> germline VCFs (PHI), but the `.gitignore` is an allow-list that tracks *only* code and
> non-patient reference config. Do not commit `*.vcf.gz`, `*.candidatos`, or any
> run artifact. See [Data privacy](#data-privacy).

---

## Contents

| File | Purpose |
|------|---------|
| `vep_annotate.sh` | Annotate a germline VCF with Ensembl VEP (plugins + custom gnomAD & ClinVar), splitting multiallelics first. Produces `*.germline.vep.vcf.gz`. |
| `filtering_r.pl` | The filtering algorithm. Reads the annotated VCF, applies gates, writes `<proband>.<panel>.candidatos`. Also the **single-variant consult** entry point (`-v`, or `--lookup` for a pre-annotated VCF): annotate one or a few variants (coords or HGVS) from scratch and report everything, gates bypassed. |
| `parse_pangolin.pl` | Convert Pangolin output into a per-variant splice-score map (`<proband>.<panel>.pangolin.tsv`). |
| `run_filtering.sh` | End-to-end driver: emit candidates → score with Pangolin → final filtering. |
| `site.sh` | One place for every external path (VEP, plugin data, Pangolin, ClinVar AA tables). Override in an untracked `site.env` — see [Setup](#setup). |
| `test/test_filtering.sh` | Regression test on synthetic data: self-tests, reference-file integrity, and end-to-end gating. No VEP, no data, no GPU, ~5s. |
| `LICENSE` | MIT, with a note that this is a triage tool requiring professional review. |
| `g4e-2026.txt` | Gene panel: `gene⇥Association⇥MOI⇥GDV`. Restricts output to panel genes; supplies MOI. Source: **Genes4Epilepsy v2026-03** (bahlolab/Genes4Epilepsy), 1078 genes; provenance header in the file. |
| `typevar.txt` | Consequence whitelist (atomic terms; matched per `&`-separated sub-term). |
| `mane-plus-clinical-names.txt` | MANE Select + MANE Plus Clinical transcript IDs; only these transcripts are considered. |
| `acmg_sf_v3.2.txt` | ACMG SF v3.2 secondary-findings genes (81): `gene⇥condition⇥MOI⇥report_category`. Always scanned. |
| `gnomad-mis-constraint.txt` | gnomAD v4.1.1 missense constraint per MANE gene: `gene⇥mis.oe⇥mis.z⇥flags`. Drives the ACMG **PP2** criterion. |

---

## Pipeline overview

```
            patient germline VCF (per individual: proband, mother, father)
                                   │
                  vep_annotate.sh  │   bcftools norm -m-any  +  VEP
                                   ▼
                *.germline.vep.vcf.gz   (CSQ: VEP + gnomAD + ClinVar + LOFTEE …)
                                   │
        run_filtering.sh ─────────┼─────────────────────────────────────────────
                                   │
   Pass 1  filtering_r.pl  ──►  <proband>.<panel>.pangolin_input.csv   (structural-pass variants)
                                   │
   Pangolin (GPU, de novo)  ──►  <proband>.<panel>.pangolin.csv  ──►  <proband>.<panel>.pangolin.tsv   (scratch, auto-removed)
                                   │
   Pass 2  filtering_r.pl  ──►  <proband>.<panel>.candidatos   (final, curation-ready)
```

`filtering_r.pl` auto-discovers families by **filename**, using an explicit **role-suffix
convention**: `<FAMILY>-P` = proband, `<FAMILY>-M` = mother, `<FAMILY>-F` = father (e.g.
`EPID107-P`, `EPID107-M`, `EPID107-F` form one trio; `EPIC280-P` + `EPIC280-M` a duo). It
globs `*.germline.vep.vcf.gz`, groups by the shared `<FAMILY>` prefix, analyzes each `-P`
sample as a proband, and pairs it with its `-M`/`-F` parents. A name not ending in
`-P`/`-M`/`-F` is ignored by auto-discovery (still usable via `--proband`). The discovery
logic has a built-in self-test: `perl filtering_r.pl --selftest`.

> **Singleton cohorts (no family structure):** if the inputs carry **no** `-P`/`-M`/`-F`
> suffixes at all (e.g. `EPIGEN01…20.germline.vep.vcf.gz`), auto-discovery would find zero
> probands — so a **singleton fallback** kicks in and analyzes **every** sample as a standalone
> proband (`inheritance = NA`), printing a `NOTE:` to that effect. This means a plainly-named
> singleton cohort "just works" with no `--proband` needed. (To force trio/duo analysis, name
> files with the `-P`/`-M`/`-F` convention or pass `--proband`.)

It is a two-pass design — if the Pangolin score map is missing it emits the candidate list
and stops; once scores exist it produces the final table.

You can **override** which sample is the proband (see [Forcing a proband](#forcing-a-proband)).

---

## Annotation (`vep_annotate.sh`)

```
bash vep_annotate.sh <input.vcf[.gz]> <output.germline.vep.vcf.gz>
```

- Splits multiallelic sites **and left-aligns** (`bcftools norm -m-any -f $REF_FASTA`) so downstream
  `chr-pos-ref-alt` keys are unambiguous **and minimal**. Left-alignment is not cosmetic: both custom
  sources below join with `type=exact`, which matches on position *and* allele, and both store indels
  minimally. Splitting without `-f` leaves the split alleles in the parent record's padded
  representation (`chr1 100 AT ATT,A` → `AT>ATT`, where gnomAD holds the same allele at 101 as
  `T>TT`), so the join silently misses — the variant then reports `AC=0`, which this pipeline reads as
  *absent from gnomAD*, and collects a spurious **PM2**. Indels were affected ~1.6× as often as SNVs
  before this was fixed. The script warns loudly and falls back to split-only if `REF_FASTA` is unset.
- Ensembl VEP (offline cache, GRCh38) with plugins **LOFTEE, REVEL, AlphaMissense, EVE,
  CADD**, `--mane_select`, HGVS, etc.
- `--custom` **gnomAD v4.1 joint** → `gnomADmin_AC_joint / AN_joint / AF_joint /
  nhomalt_joint / FILTER`.
- `--custom` **ClinVar** (NCBI weekly VCF, chr-prefixed) → `ClinVar_CLNSIG / CLNREVSTAT /
  CLNDN`, kept current independently of the VEP cache's older bundled `CLIN_SIG`.

All annotations land inside the `CSQ` INFO field; `filtering_r.pl` resolves them **by
name** from the CSQ header (no hard-coded column indices).

---

## The filtering algorithm (`filtering_r.pl`)

Logic runs **per transcript annotation** of each variant, then collapses to **one MANE
row per variant *per gene*** (a variant hitting >1 MANE transcript of the **same** gene — MANE Select
+ MANE Plus Clinical — keeps a single row: panel-primary first, then MANE Select, then the most
evidence arms). The collapse is keyed on **(variant, gene)**, not on coordinates alone: overlapping
MANE gene models are common (MYH11+NDE1, HPDL+MUTYH, COL4A1+COL4A2, SETD1A+STX1B), and keying on
position let a panel candidate silently delete a reportable **ACMG-SF incidental** at the same
position — and corrupted the other gene's comp-het tally by reassigning the row's gene.
`--lookup` consults still report every annotation.

### Stage 1 — Structural gates (ALL required, AND)

| Gate | Source | Rule |
|------|--------|------|
| MANE transcript | `mane-plus-clinical-names.txt` | CSQ `Feature` ∈ MANE set |
| Consequence | `typevar.txt` | consequence split on `&`; kept if **any** sub-term is whitelisted |
| Gene panel | `g4e-2026.txt` (default) or a custom genes-of-interest file | CSQ `SYMBOL` ∈ panel |
| Rarity (MOI-aware) | gnomAD joint AC/AN | AF = AC/AN×100 ≤ threshold: **dominant `$FREQ_AD`=0.01%**, **recessive `$FREQ_AR`=1.0%** (MOI contains AR/XLR) |

> **ClinVar exemption.** A panel-gene variant classified **Pathogenic/Likely-pathogenic in ClinVar
> with ≥1 review star** satisfies Stage 1 regardless of the **rarity** and **consequence** gates
> (MANE and panel membership still apply). An established classification is evidence in its own
> right, and both gates would otherwise discard it before any evidence is read: a **founder allele**
> in a bottlenecked population can sit above the frequency ceiling, and a **pathogenic non-coding or
> synonymous** variant carries a consequence term that is not whitelisted. The variant is reported
> with whatever frequency-based benign criteria it earns (`BS1`/`BA1`), so the curator sees the
> tension rather than losing the variant. This mirrors the ACMG-SF path, which exempts the same
> tier for the same reason.

### Stage 2 — Inclusion / rescue gate (at least ONE, OR)

A surviving variant must trip **one or more** of these. Each is independent; a `kept_by`
column records which fired.

The **`kept_by` token** column is the literal string written to the output — grep on that, not on
the arm's prose name.

| Arm | `kept_by` token | Threshold |
|-----|-----------------|-----------|
| CADD | `CADD` | `$CADD_MIN` = 25.3 |
| AlphaMissense | `AM` | `am_score` ≥ `$AM_MIN` = 0.792 (ClinGen PP3) |
| EVE pathogenic | `EVE` | `eve_class` is Pathogenic |
| REVEL | `REVEL` | `$REVEL_MIN` = 0.644 (ClinGen PP3) |
| Pangolin (splice) | `Pangolin` | `$SPLICE_MIN` = 0.5 (max \|Δscore\|) |
| ClinVar P/LP | `ClinVar` | `ClinVar_CLNSIG` Pathogenic/Likely_pathogenic (excludes Conflicting & Benign) |
| PS1 / PM5 | `PS1` / `PM5` | ClinVar amino-acid match (≥1★): **PS1** = a *different* variant giving the same AA change is P/LP, **PM5** = a different change at the same residue is P/LP. A **single-codon in-frame deletion** of the residue also triggers PM5 (a different protein change at the same P/LP residue; tagged `(in-frame del)`). Rescues the variant even when CADD/AM/REVEL miss it; the `clinvar_aa` column carries the detail (and any `(conflicting)` flag). |
| LoF | `LoF` | LOFTEE `LoF=HC`, or a high-impact truncating consequence (frameshift / stop_gained / splice_donor / splice_acceptor / start_lost) unless LOFTEE downgraded it to `LC`. Covers truncating indels that CADD (SNV-only) and the missense predictors miss. |
| AR_hom / AR_hem | `AR_hom` / `AR_hem` | **Homozygous or hemizygous, protein-altering** (`missense` / `inframe_*` / `stop_lost` / `start_lost` / `protein_altering_variant`) variant in a gene whose panel **MOI contains AR or XLR** — this includes dual `AD, AR` genes, unlike the *pure*-recessive carrier logic below — with clean allele balance (**AB > 0.75**). A **hemizygous** male call on non-PAR chrX qualifies on the same rationale and is tagged `AR_hem`. Rescued even without in-silico/ClinVar support. *Rationale:* under recessive inheritance a **biallelic (homozygous) genotype in a disease gene is itself pathogenicity evidence**, independent of missense predictors — which are calibrated largely on dominant/heterozygous effects and can miss true recessive alleles. **Restricted to coding changes** (so it doesn't flood on benign homozygous intronic/polypyrimidine variants — those use the Pangolin arm; truncating LoF uses the LoF arm) **and to clean homozygous/hemizygous calls** (AB > 0.75 guards against false-hom artifacts). Already rare (AR freq gate), MANE, in-panel by this point; `BS1`/`BS2`/`BA1` still flag benign-leaning ones. |

All thresholds are single constants at the top of `filtering_r.pl`.

### Genotype-aware annotation (not gates)

- **Zygosity / GT / DP / GQ / AB** are read from the proband `FORMAT`/sample column.
  **Haploid genotypes** (`GT=1`, emitted by DRAGEN for non-PAR chrX/chrY in males and for chrM)
  report `zygosity=hem` and are treated as **biallelic-equivalent**: a hemizygous call is a complete
  genotype, not half a het. It earns the `HEM` flag and the `AR_hem` rescue arm, gets the `AB_hom` QC
  check, and is never dropped as a carrier. A **hemizygous parent counts as a carrier**, so an
  inherited X-linked variant is not misreported as de novo.
- **Inheritance** uses *parental genotype* (carrier = non-ref GT, not mere site presence):
  `IB / IM / IF / DN` in a full trio; duo-ambiguous `DN/IF` (mother-only) or `DN/IM`
  (father-only); `NA` for a singleton.
- **`flags`** — one `;`-separated column carrying the recessive verdict **first**, then the QC /
  artifact warnings (they used to be two columns, `recessive_flag` and `qc_flag`, four columns apart).
- **The recessive verdict** — the recessive picture is worked out **per gene** (for **recessive-capable
  genes only** — AR/XLR/dual; a purely dominant gene never gets one, so two independent hets are not
  mislabeled comp-het), but the label written to each row **describes that row**: `HOM` on a
  homozygous row, `HEM` on a hemizygous row, `CompHet(trans)` (≥2 het variants phaseable to opposite parents — trio only) or
  `CompHet?` (≥2 het, unphaseable — e.g. duo/singleton) on the heterozygous rows that constitute it,
  or `carrier-only` (see below). A heterozygous variant sitting in a gene that is homozygous for some
  *other* variant is therefore left blank rather than labelled `HOM` — the flag never contradicts the
  row's own `zygosity`. The gene-level verdict still governs which rows survive the carrier drop.
- **Dual-inheritance genes** (panel MOI lists **both** AD and AR, e.g. `AD, AR`, **and plain `XL`**)
  are treated as **dominant** for the carrier logic: a **solitary het passes through as a normal
  candidate** (no recessive flag), while a genuine `HOM`/`HEM`/comp-het still gets the recessive flag.
  Plain `XL` is the Genes4Epilepsy vocabulary for an X-linked gene with no XLD/XLR split — **72 of the
  1078** g4e-2026 genes, including CDKL5, MECP2, ARX, IQSEC2, PCDH19, DDX3X, ATRX, SLC6A8 and FLNA.
  Before 2026-08 it matched *neither* MOI predicate, so those genes silently got no HOM/CompHet flag,
  no `AR_hom` rescue, and the strict dominant AF ceiling. `XLR` / `XLD` keep their specific meanings. This
  prevents dropping a dominant-acting variant (e.g. a LoF) just because the gene *also* has a
  recessive mechanism. Only **pure** AR/XLR genes use the carrier path below.
- **Recessive carrier drop (default):** in a **pure** recessive (AR/XLR) panel gene, a **solitary het**
  that is not biallelic (neither `HOM` nor comp-het) is **dropped** — a single het cannot explain a
  recessive disease, and carrier states are noise for clinical interpretation. **This does not affect
  true compound hets:** a gene with ≥2 gate-passing hets is biallelic (`CompHet` flag) and every such
  row is kept regardless. The same drop covers recessive ACMG SF genes
  (see [Secondary findings](#secondary-findings-acmg-sf-v32)).
  - **`--keep-ar-carriers`** (or env `KEEP_AR_CARRIERS=1`) is the **opt-in** for a targeted second-hit
    hunt: it **surfaces the strong** solitary AR/XLR carriers — the **carrier-only tier**, kept only if
    they clear a strong-evidence bar (ClinVar P/LP ≥1★, LOFTEE-HC, or ≥2 strong predictors: AM ≥ 0.906,
    CADD ≥ 28.1, EVE pathogenic, REVEL ≥ 0.773) **and** are not Benign/Likely-benign — flagged
    `flags=carrier-only`. Use it when you suspect the exome missed a second allele (deep-intronic,
    CNV, regulatory) in a specific gene/case.
    ⚠️ A **common** SNP (high gnomAD AF) is not a valid second hit — an apparent "comp-het" of one rare
    plus one common variant is really a **solitary carrier** of the rare allele, not a biallelic genotype.

### Stage 2b — Splice discovery probes

`typevar.txt` contains no bare `intron_variant` and no bare `synonymous_variant`. Those variants were
therefore dropped at Stage 1, never reached Pangolin, and the splice rescue arm could only ever
**upgrade** a variant that was already whitelisted as splice-region — never **discover** one. The
classic pathogenic deep-intronic alleles (CFTR `c.3718-2477C>T`, USH2A `c.7595-2144A>G`) were
structurally unreachable, and BP7's entire target population never reached the classifier.

A **probe** is an intronic or synonymous variant in a panel or ACMG-SF gene that is sent to Pangolin
purely so the model can look at it. A probe is **not a candidate**: it earns a row only if the splice
arm fires (`pangolin_score` ≥ `$SPLICE_MIN`), and otherwise disappears. The probe set widens what can
be *found* without widening the table.

Three bounds keep it affordable and honest:

| Bound | Value | Why |
|---|---|---|
| Intron distance | `$INTRON_MAX_DIST` = 300 bp | Read from the HGVSc offset (`c.1234+56` → 56; a range keeps its closest endpoint). Pangolin scores ~8 variants/s on one GPU and a WGS proband carries ~155k rare intronic variants in panel genes, so distance is the primary volume bound |
| Rarity | `$PROBE_FREQ_MAX` = 0.01% | Applied to **every** gene. A pathogenic splice variant is rare regardless of mode of inheritance, so the permissive `$FREQ_AR` carrier ceiling (1%, meant for recessive coding candidates) is deliberately not used |
| gnomAD coverage | `AN > 0` (`$PROBE_REQUIRE_GNOMAD`, default on) | See below — this one is load-bearing. `--probe-uncovered` turns it off |

> ⚠️ **The custom gnomAD resource is MANE-restricted.** `gnomAD.joint.v4.1.mane.all.vcf.gz` covers
> MANE transcripts and their flanks, *not* deep intronic sequence. An uncovered position yields
> `AC=""`/`AN=""`, which the code coerces to 0 — so `$freq` computes as 0, passes **every** rarity
> ceiling, and **PM2 fires on what is really an annotation gap**. Probing uncovered territory would
> rescue variants with no working frequency filter *and* a manufactured pathogenic criterion.
> Requiring `AN > 0` bounds the probe set to where the resource can actually answer the question.
>
> **Measured consequence.** On a real WGS proband the 300 bp window holds ~22,700 rare intronic
> variants, of which only **~8** are gnomAD-covered. So with the default the probe set costs almost
> nothing (8 probes/proband) — and discovery stays confined to roughly the MANE footprint, meaning a
> classic deep-intronic allele such as CFTR `c.3718-2477C>T` is **still out of reach**.
> `--probe-uncovered` (env `PROBE_UNCOVERED=1`) probes all ~22,700: ~48 min of GPU time per proband,
> in territory where the rarity gate cannot work. PM2 is withheld on those rows so the run does not
> invent a pathogenic criterion, but the frequency filter is genuinely blind there.
>
> **The real unlock is rebuilding the custom gnomAD VCF with genome-wide coverage.** That is a
> pending resource task, not a code change; once done, this flag stops mattering.

ACMG-SF genes are also sent to Pangolin now. Previously only panel candidates were scored, so every
`GDV=Incidental` row carried a blank `pangolin_score` and could never earn BP7 or the splice rescue —
a secondary finding was structurally denied evidence a primary candidate got.

Disable the whole probe set with `--no-splice-discovery` (or `NO_SPLICE_DISCOVERY=1`); the EMIT
summary reports how many probes were added and the final summary how many Pangolin rescued.

### Stage 3 — Cohort recurrent-artifact filter (internal panel-of-normals)

When a **real cohort** is auto-analyzed together (≥ `$COHORT_MIN` = **5 samples**, ≥2 probands), the pipeline
builds an internal *panel of normals*: a one-pass genotype tally (no CSQ, so it is cheap) over all
input samples counting, per `chr-pos-ref-alt`, how many carry the ALT and their zygosity breakdown.
A candidate is then **dropped** when it is **both**:

1. **carried by ≥ `$COHORT_MIN_CARRIERS` (3) samples** — an absolute floor, and
2. **cohort-recurrent** — carried by ≥ `$COHORT_MAX_FRAC` (25%) of the cohort, and
3. **gnomAD-absent** — gnomAD joint **AC = 0**, i.e. the allele is wholly unobserved.

This targets systematic technical artifacts — reference/mapping errors in paralog-rich or
low-complexity genes (e.g. the recurrent `SYNE1` / `KMT2C` sites that appear, gnomAD-absent, in a
large fraction of a cohort). **Both conditions are required, and that conjunction is the
discriminator:** recurrence alone would remove true **founder / bottleneck** alleles, but a real
founder allele frequent enough to reach a quarter of the cohort leaves a footprint in gnomAD's large
Admixed-American sample — so requiring gnomAD-absence keeps the filter founder-safe (important for
under-represented ancestries). Absence is tested as **AC = 0** rather than against a frequency
ceiling, because a ceiling anywhere near `$FREQ_AD` is already implied by the Stage 1 rarity gate for
dominant genes: recurrence would then be the only condition still doing work, and the founder
protection would be nominal. Requiring a literal zero keeps the second condition meaningful under
every mode of inheritance — **any** gnomAD observation, even a single allele, spares the variant.
The **absolute carrier floor** is what makes this safe at real batch sizes. At N=8 the 25% fraction
alone means *two* carriers, and two unrelated probands sharing a gnomAD-absent allele is an ordinary
founder/relatedness event, not evidence of a technical artifact. The floor binds on small batches;
the fraction binds on large cohorts (at N=56, 25% = 14 carriers, well above the floor).

Every drop is logged (`cohort_artifact drop:` lines with carrier count + hom/het breakdown, and a
per-proband total).

> **Thresholds revised 2026-08.** `$COHORT_MIN` was **10 probands** — a bar no internal batch reaches,
> since they run 6–9 samples — so the filter had **never fired on a real clinical run**, and the
> recurrent KMT2C / SYNE1 mismapping artifacts reached every delivered table. Measured on batch4
> (8 singleton probands) the carrier distribution is cleanly bimodal: the four artifacts sit at 5, 6, 7
> and 8 of 8 carriers (63–100%), and the next most recurrent candidate is 2 of 8 (25%). Both new
> thresholds fall inside that gap. Eligibility and the denominator are now both counted in **samples**;
> previously eligibility counted probands while the denominator counted all VCFs, which diluted the
> acting fraction ~3× on trio cohorts (an artifact in 4 of 10 probands scored 4/30 = 13%).
>
> ⚠️ **The filter only sees the samples in ONE run.** A batch split across per-sample directories has
> a cohort of 1 and the filter stays off. Put the whole batch's VCFs in one `$WORKDIR`.

- **OFF** for single-variant (`-v` / `--lookup`), forced/single proband (`--proband`), and any run of
  fewer than `$COHORT_MIN` samples — a per-variant or per-proband consult has no cohort to compare
  against.
- **`--keep-cohort-artifacts`** (or env `KEEP_COHORT_ARTIFACTS=1`) keeps them, tagged
  `flags=cohort_artifact`, instead of dropping — for a founder-enriched cohort, review the drop log
  (a genuinely private founder allele would surface there and can be kept this way).
- Self-test (synthetic cohort, no reference files needed): `perl filtering_r.pl --selftest-cohort`.
- Thresholds are the `$COHORT_MIN` / `$COHORT_MAX_FRAC` / `$COHORT_MIN_CARRIERS` constants at the
  top of `filtering_r.pl`.

A per-proband **run summary** prints counts (read / multiallelic-skipped / structural-pass
/ candidates / cohort-artifacts dropped) and breakdowns by `kept_by` and inheritance.

### Output columns (`<proband>.<panel>.candidatos`, TSV)

`chr, start, end, ref, alt, gene, strand, consequence, hgvs,
revel, eve_class, eve_score, cadd, am_class, am_score, pangolin_score,
clinvar_sig, clinvar_stars, clinvar_disease, clinvar_aa, loftee,
gnomAD_ac, gnomAD_an, gnomAD_af, gnomAD_nhomalt, gnomAD_filter,
zygosity, GT, DP, GQ, AB, GT_SOURCE, NCALLERS, CONF,
inheritance, kept_by,
acmg_class, acmg_criteria, flags, Association, MOI, GDV`

A run also writes **`batch.<panel>.candidatos`** — every proband's rows in one table, prefixed with a
`sample` column; columns 2..N are byte-identical to the per-proband header, so anything that reads one
reads the other. Written once, after all probands are processed.

`start`/`end` are 1-based and span the **REF allele** (`end` = `start` + len(`ref`) − 1), so an indel
reports its full footprint and the row round-trips into the 5-field `chr-start-end-ref-alt` form the
consult mode accepts. For an SNV the two are equal.

`hgvs` combines HGVSc and HGVSp as `TRANSCRIPT:c.… (p.…)` (the `ENSP…:` protein-accession
prefix is stripped; non-coding/synonymous variants show only the `c.` part).

- **`clinvar_stars`** — review-status stars (0–4) of the exact variant's ClinVar classification.
- **`clinvar_aa`** — the PS1/PM5 amino-acid evidence string (see below), e.g.
  `PS1:BRCA1 p.R1699W (3*)` or `PM5:… [P/LP at residue: …] |conflicting`; empty if none.
- **`GT_SOURCE` / `NCALLERS` / `CONF`** — consensus provenance, populated only when the input
  came from `consensus.sh` (the Sarek union consensus); empty for single-source (e.g. DRAGEN) VCFs.
  `GT_SOURCE` names the caller the genotype was taken from (`deepvariant`, `strelka`,
  `haplotypecaller`); a non-DeepVariant value also means **no `VAF`** (allele balance is computed
  from `AD`) and raises a `GT_rescued` QC flag.

### Automated ACMG/AMP classification & QC flags

- **`acmg_class` / `acmg_criteria`** — a **triage** classification per variant
  (Pathogenic / Likely_pathogenic / VUS / Likely_benign / Benign / Conflicting), combined per the
  **categorical ACMG 2015 rules** from the criteria the pipeline evaluates automatically:

  **Pathogenic**

  | Criterion | What triggers it | Source |
  |---|---|---|
  | **PVS1** | LoF: LOFTEE = HC, or a truncating consequence with LOFTEE ≠ LC | VEP / LOFTEE |
  | **PS1** | A **different** variant producing the same amino-acid change is ClinVar P/LP (≥1★). The variant's own ClinVar record is excluded, so a variant that is itself P/LP does not earn PS1 from its own submission | ClinVar MANE-missense |
  | **PS2** | De novo in a **full trio** (`inheritance=DN`, clean proband genotype); relatedness is assumed confirmed. Structurally unreachable without both parents — `inheritance` is only ever `DN` when both are present, a singleton gets `NA` and a duo gets `DN/IM`–`DN/IF` | parental GT |
  | **PM2** | Absent or singleton in gnomAD (AC ≤ 1), counted at **Supporting** (`PM2_Supporting`) per [ClinGen SVI 2020](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/), not ACMG 2015's Moderate. See the caveat below | gnomAD v4.1 |
  | **PM4** | Protein length change (in-frame indel / `stop_lost`). **Not counted when PVS1 fired** — VEP compound terms (`start_lost&inframe_deletion`, `frameshift_variant&stop_lost`) otherwise yielded two ACMG lines for one protein-terminus effect, pushing an LP call to Pathogenic | consequence |
  | **PM5** | Different change — **or a single-codon in-frame deletion** — at a residue carrying a P/LP missense (≥1★) | ClinVar MANE-missense |
  | **PM6** | **Assumed** de novo: a trio `DN` whose genotype isn't clean, or a duo-ambiguous `DN/IF`–`DN/IM` **in a gene whose panel MOI contains AD or XL**. A duo-ambiguous call in a pure-AR gene — or under any panel with `MOI = NA`, e.g. a plain-symbol custom list — never earns PM6 | parental GT + panel MOI |
  | **PP2** | **Missense** in a gene with low benign-missense variation — gnomAD v4.1.1 missense constraint `mis.oe < 0.6` (MANE; constraint outliers excluded). Counts **independently of PP3** (both are legitimate separate ACMG lines — gene-level intolerance vs variant-level prediction), but **suppressed when BP4 fires** (a benign-predicted variant gets no gene-level pathogenic support). | gnomAD v4.1.1 constraint |
  | **PP3** | Computational damaging, graded Supporting/Moderate/Strong (see below) | AlphaMissense / REVEL |
  | **PP5** | This variant is reported pathogenic in ClinVar **with ≥1 review star**. The star gate matches every other ClinVar consumer in the pipeline; without it a single 0-star "no assertion criteria provided" submission (~16% of the P/LP corpus) supplied the criterion that lifts LP to Pathogenic | ClinVar |

  **Benign**

  | Criterion | What triggers it | Source |
  |---|---|---|
  | **BA1** | gnomAD AF ≥ 5% | gnomAD v4.1 |
  | **BS1** | gnomAD AF ≥ 1% (and < 5%) | gnomAD v4.1 |
  | **BS2** | ≥ 10 homozygotes in gnomAD | gnomAD v4.1 |
  | **BP4** | Computational benign, graded (see below) | AlphaMissense / REVEL |
  | **BP6** | This variant is reported benign in ClinVar (≥1★) | ClinVar |
  | **BP7** | Synonymous **and scored** by Pangolin at < 0.2. The score must exist: an unscored variant is unknown, not benign, so BP7 is withheld rather than assumed | Pangolin |

  **Not evaluated (manual curation only):** PS3/BS3 (functional), PS4 (case-control), PM1
  (hotspot/domain), PM3 (in trans), PP1/BS4 (segregation), PP4 (phenotype specificity),
  BP1/BP2/BP3/BP5.

  **PS1/PM5** use the ClinVar MANE-missense resource (`clinvar.MANE_missense.{PLP,BLB}.tsv`), matched on
  gene + protein residue + amino-acid change, requiring **≥1 review star**. The resource is indexed by
  the **source variant** (`chr-pos-ref-alt`) behind each amino-acid change, which is what lets PS1
  exclude the record belonging to the variant being classified — PS1 rests on a *previously
  established* variant, so a self-match would count one submission twice, as PS1 and again as PP5.
  PM5 is unaffected: a different amino-acid change is necessarily a different variant. A match that is
  also reported B/LB is tagged **`(conflicting)`** (still counted — flag for manual review) and
  detailed in the `clinvar_aa` column.
  **PP3/BP4 come from a single calibrated predictor**, graded **Supporting/Moderate/Strong**:
  **AlphaMissense** primary ([Bergquist 2025](https://doi.org/10.1016/j.gim.2025.101402): PP3
  supp ≥0.792 / mod ≥0.906 / strong ≥0.990; BP4 supp ≤0.169 / mod ≤0.099), **REVEL** fallback
  ([Pejaver 2022](https://doi.org/10.1016/j.ajhg.2022.10.013): PP3 supp ≥0.644 / mod ≥0.773 /
  strong ≥0.932; BP4 supp ≤0.290 / mod ≤0.183 / strong ≤0.016) — with a **REVEL direction-conflict
  veto**, mapped to the 2015 tiers (BP4_Moderate → supporting-benign, since 2015 has no benign-Moderate).
  **Not a final clinical call**: PM1 not assessed; PP2 is gene-level constraint only (no domain/hotspot
  reasoning); PVS1 doesn't verify gene mechanism/NMD; PS1/PM5 rely on ClinVar AA matching (no independent
  re-curation; PS1 and PP5 can still co-occur when a *different* variant supplies the amino-acid match).
- **QC / artifact components of `flags`:** `lowDP` (<`$QC_MIN_DP`), `lowGQ` (<`$QC_MIN_GQ`),
  `AB_het`/`AB_hom` (skewed allele balance; `AB_hom` also covers hemizygous calls),
  `homopolymer` (indel in **or adjacent to** a ≥5 bp
  homopolymer — the reference is scanned ±12 bp around the position, so a nearby run also flags —
  error-prone),
  `GT_rescued` (genotype borrowed from a non-DeepVariant caller via `consensus.sh`; no VAF),
  `inh_lowqual` (carrying-parent genotype is weak), `DN_unconfirmed`, `cohort_artifact` (recurrent
  gnomAD-absent cohort artifact, present only under `--keep-cohort-artifacts` — otherwise dropped),
  `clinvar_conflict` (see below).
- **`clinvar_conflict`** — the auto-class reached Pathogenic/Likely pathogenic **while a hard benign
  line fired** (`BP6`/`BS1`/`BS2`/`BA1`). The categorical `Conflicting` verdict requires *both* sides to
  reach a 2-tier threshold independently, so a single benign criterion never blocks a pathogenic call:
  `PVS1,PM2,BP6` reads as Likely_pathogenic on a variant ClinVar calls **Benign with review stars**,
  and a curator sorting by `acmg_class` sees a clean LP. This is a triage tool, so the class is left
  alone and the contradiction is made visible instead. `BP4` is deliberately excluded — a computational
  prediction disagreeing with PVS1/PM2 is routine, not a contradiction.
- **De-novo confidence [#6]:** parent VCFs here are *variant-only* (no reference depth at non-variant
  sites), so de-novo cannot be confirmed from parental coverage — `DN` rows are flagged
  `DN_unconfirmed`. Inherited rows instead get `inh_lowqual` when the parental call is low quality.
  True de-novo confirmation needs parental gVCFs/BAMs.

---

## Secondary findings (ACMG SF v3.2)

The **81 ACMG SF v3.2 genes** (`acmg_sf_v3.2.txt`) are **always** scanned, independent of the
candidate `-l`/`--list` panel, with a **stricter** gate than candidates. Findings are written into
the **same** `.candidatos` output, flagged **`GDV = Incidental`** (with `Association`/`MOI` from the
ACMG table and `kept_by` = the evidence tier). Curators split primary vs secondary on the GDV column.

Inclusion (any one):
- **`ClinVar_P/LP`** — ClinVar Pathogenic/Likely-pathogenic with **≥1 review star** (frequency-agnostic,
  so known founder alleles are not lost). *Known / directly reportable.*
- **`LoF`** — novel LOFTEE-HC. *Expected pathogenic (review-queue; verify gene mechanism).*
- **`Computational`** — **≥2 of** AM ≥ 0.906, CADD ≥ 28.1, EVE Pathogenic, REVEL ≥ 0.773
  (rarity-capped). *Candidate SF requiring expert classification — not auto-reportable.*

Gene-specific rules from the ACMG table are honored: `TTN` truncating-only, `HFE` C282Y-homozygotes-only,
and recessive (AR) genes report biallelic (hom or comp-het) findings; a solitary het in a recessive SF
gene is dropped by default (surfaced via `--keep-ar-carriers`, the **carrier-only** tier, like primary
genes). Thresholds are `$SF_*` constants in `filtering_r.pl`.

> ⚠️ Secondary findings carry distinct **consent / reporting** obligations — handle per your lab policy.

---

## Splice scoring (Pangolin)

There is no official Pangolin VEP plugin, so splice scores are computed **de novo** by
standalone [Pangolin](https://github.com/tkzeng/Pangolin) on GPU and merged by position.
Only the proband's **structural-pass** variants are scored (a few hundred), not the whole
VCF. `parse_pangolin.pl` reduces each variant to `max(|increase|, |decrease|)`.

> All Pangolin scratch (`*.pangolin_input.csv`, `*.pangolin.csv`, `*.pangolin.tsv`) is
> **deleted after Pass 2 writes the final table** — the pipeline keeps only
> `<proband>.<panel>.candidatos` and the annotated VCFs (`*.germline.vep.vcf.gz` + `.tbi` +
> `_summary.html`). Pangolin is therefore recomputed on every run; it is cheap because only
> the few hundred structural-pass candidates are scored. (Cleanup runs only on success, so a
> failed run leaves intermediates in place for debugging.)

---

## Setup

> **Setting this up on a new machine? Follow 0.1 → 0.8 below in order.** It is a long
> one-time job: the annotation resources total roughly **200 GB**. Nothing here is tied to one
> machine — every path is configurable (step 0.2), so the repo can live anywhere and the data
> wherever you have room.

### 0.1 — Clone and check the parts that need nothing

```bash
git clone https://github.com/edoper/candidate-filtering.git
cd candidate-filtering
./test/test_filtering.sh          # ~5s: needs only perl + bgzip, no VEP, no data
```

That exercises the filtering logic on synthetic data. It passing means the algorithm and the
tracked gene panels are intact — you can then add the annotation resources below.

### 0.2 — Tell the repo where your data lives

All paths come from `site.sh`. Override any of them by exporting, or by creating an untracked
**`site.env`** beside it — that file is where your own layout belongs, and it is never committed:

```bash
# site.env
VEP=/opt/ensembl-vep/vep                 # the vep executable
VEP_DATA=/data/vep_cache                 # --dir_cache (the offline GRCh38 cache)
VEP_REFS=/data/vep_refs                  # plugin data + custom VCFs (steps 0.4-0.5)
VEP_PLUGINS=$HOME/.vep/Plugins           # --dir_plugins (must contain loftee/)
CLINVAR_AA_DIR=/data/clinvar-aa          # step 0.6 (optional; enables PS1/PM5)
CONDA_BASE=$HOME/miniconda3              # step 0.7
```

Defaults (`$HOME/ensembl-vep`, `$HOME/vep_data`, `$HOME/vep_refs`, `$HOME/.vep/Plugins`) are the
layout this pipeline was developed against.

### 0.3 — Ensembl VEP + the GRCh38 cache

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep && perl INSTALL.pl        # choose the homo_sapiens_vep GRCh38 cache (~25 GB)
```

`INSTALL.pl` also installs plugins; select **LoF (LOFTEE), REVEL, AlphaMissense, EVE, CADD**. The
filter needs `bcftools` too (`conda install -c bioconda bcftools htslib`).

**LOFTEE also needs `Bio::DB::HTS` and the htslib shared library it links against.** `vep_annotate.sh`
puts `$PERL5LIB_EXTRA` on `PERL5LIB` and `LD_PRELOAD`s `$HTSLIB_SO`; both default to the developed-against
layout (`$VEP_PLUGINS/loftee:$HOME/perl5/lib/perl5` and `$HOME/htslib/libhts.so`). If your htslib lives
elsewhere, set `HTSLIB_SO` in `site.env` — the preload is skipped silently when the file is missing, and
LOFTEE then fails at run time rather than at startup.

### 0.4 — Plugin data files

Each plugin needs its own dataset, all GRCh38, under `$VEP_REFS`. They are large and each has its
own licence — check the terms for your use (AlphaMissense and REVEL are free for academic use;
CADD requires a licence for commercial use).

| Under `$VEP_REFS/` | What | Where from |
|---|---|---|
| `CADD/whole_genome_SNVs.tsv.gz` (+ `.tbi`) | CADD v1.7 SNVs (~80 GB) | [cadd.gs.washington.edu/download](https://cadd.gs.washington.edu/download) |
| `CADD/gnomad.genomes.r4.0.indel.tsv.gz` (+ `.tbi`) | CADD v1.7 indels | same |
| `REVEL/new_tabbed_revel_grch38.tsv.gz` (+ `.tbi`) | REVEL | [sites.google.com/site/revelgenomics](https://sites.google.com/site/revelgenomics/downloads) |
| `AlphaMissense/AlphaMissense_hg38.tsv.gz` (+ `.tbi`) | AlphaMissense | [Zenodo 10813168](https://zenodo.org/records/10813168) |
| `EVE/eve_merged.vcf.gz` (+ `.tbi`) | EVE | [evemodel.org](https://evemodel.org/) |
| `loftee/GRCh38/` | `human_ancestor.fa.gz`, `loftee.sql`, `gerp_conservation_scores…bw` | [LOFTEE GRCh38 branch](https://github.com/konradjk/loftee/tree/grch38) |

### 0.5 — The two custom VCFs (gnomAD + ClinVar)

These are `--custom` annotations, not plugins, and both must be **`chr`-prefixed, bgzipped and
tabixed** — the annotation step asserts they exist and fails early if not.

```bash
mkdir -p $VEP_REFS/gnomAD_min $VEP_REFS/clinvar

# ClinVar (small, ~100 MB) — add the 'chr' prefix the pipeline expects
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
bcftools annotate --rename-chrs <(for i in $(seq 1 22) X Y MT; do echo -e "$i\tchr$i"; done) \
  clinvar.vcf.gz -Oz -o $VEP_REFS/clinvar/clinvar.chr.vcf.gz
tabix -p vcf $VEP_REFS/clinvar/clinvar.chr.vcf.gz

# gnomAD v4.1 joint frequencies, reduced to the fields the filter reads
#   (AC_joint, AN_joint, AF_joint, nhomalt_joint, FILTER) — the full release is ~2 TB, so
#   subset to MANE regions and strip everything else before saving.
#   -> $VEP_REFS/gnomAD_min/gnomAD.joint.v4.1.mane.all.vcf.gz  (+ .tbi)
```

> **`AC=0 / AN=0` in the output means "absent from gnomAD"** (the custom VCF is sites-only), **not**
> an uncallable region. Never read `AN=0` as evidence of a technical artifact.

### 0.6 — ClinVar amino-acid tables (optional — enables ACMG PS1/PM5)

PS1/PM5 need per-residue P/LP and B/LB missense evidence in two TSVs:

```
$CLINVAR_AA_DIR/clinvar.MANE_missense.PLP.tsv
$CLINVAR_AA_DIR/clinvar.MANE_missense.BLB.tsv
```

Build them by splitting the ClinVar VCF's MANE missense records by clinical significance, keyed by
`gene:protein_position`. **If `CLINVAR_AA_DIR` is unset or the files are missing, filtering still
runs normally and PS1/PM5 are simply skipped** with a warning — so you can defer this.

### 0.7 — Pangolin splice scoring (optional but recommended)

Needs a GPU-capable PyTorch environment:

```bash
conda create -y -n pangolin -c conda-forge python=3.10 pip
conda activate pangolin
pip install torch --index-url https://download.pytorch.org/whl/cu121
pip install pyvcf3 gffutils biopython pandas pyfastx "setuptools<81"
pip install git+https://github.com/tkzeng/Pangolin.git
```

Plus a chr-named GRCh38 primary-assembly FASTA (samtools-indexed) and the GENCODE annotation DB, by
default at `$VEP_REFS/pangolin/GRCh38.primary_assembly.genome.fa` and `…/gencode.v38.annotation.db`
(override with `PANGOLIN_FASTA` / `PANGOLIN_DB`). Without Pangolin the splice rescue arm and the
`pangolin_score` column are unavailable; everything else works — run `filtering_r.pl` directly
instead of `run_filtering.sh`.

The same FASTA doubles as `REF_FASTA` for the homopolymer QC flag; if absent, that flag is skipped.

### 0.8 — Where the input comes from

This repo starts from an **annotated** VCF. To produce the calls in the first place, see the
companion repo **[sarek-clinical](https://github.com/edoper/sarek-clinical)** (four-caller
consensus germline calling on Google Cloud). It hands `<sample>.consensus.vcf.gz` straight to
`vep_annotate.sh` here. Single-source VCFs (e.g. DRAGEN) work equally well — the filter picks up
the consensus `GT_SOURCE`/`NCALLERS`/`CONF` tags when present and ignores them otherwise.

Name inputs `<FAMILY>-P/-M/-F` (proband/mother/father) — the filename drives trio/duo
auto-discovery. Plainly-named singletons also work (each is analysed as its own proband).

---

## Usage

```bash
# 1) Annotate each family member, naming outputs with the role suffix
#    (-P proband, -M mother, -F father)
bash vep_annotate.sh EPIC280.raw.vcf.gz    EPIC280-P.germline.vep.vcf.gz
bash vep_annotate.sh EPIC280M.raw.vcf.gz   EPIC280-M.germline.vep.vcf.gz

# 2) Run the full filtering pipeline (emit → Pangolin → final)
bash run_filtering.sh
#    → EPIC280-P.g4e.candidatos
```

### `filtering_r.pl` command-line flags

The complete accepted set — **anything else is a hard error** (there is no positional argument).
Every value-taking flag also accepts the `--flag=value` form.

| Flag | Value | What it does |
|------|-------|--------------|
| `-l`, `--list` | genes file | Candidate-gene panel, replacing `g4e-2026.txt`. **This is the only way to set the panel**, and it also sets the output `<panel>` tag. |
| `-v`, `--variant` | variant | Single-variant consult; **repeatable** for several variants. Coords or `ENST…` HGVS. |
| `--lookup` | annotated VCF | Consult a pre-annotated `*.germline.vep.vcf.gz` directly. Mutually exclusive with `-v`. |
| `-p`, `--proband` | sample base-name | Force a sample as proband, overriding filename auto-discovery. **Repeatable.** |
| `--all-transcripts` | — | Consult mode: report every transcript, not just MANE. |
| `--keep-vcf` | — | Consult mode: keep the annotated VCF instead of deleting it. |
| `--no-splice` | — | Consult mode: skip the inline Pangolin run. |
| `--keep-ar-carriers` | — | Surface strong solitary AR/XLR carriers (carrier-only tier) instead of dropping them. |
| `--keep-cohort-artifacts` | — | Tag cohort recurrent artifacts `flags=cohort_artifact` instead of dropping them. |
| `--no-splice-discovery` | `NO_SPLICE_DISCOVERY=1` | Skip the intronic/synonymous Pangolin probe set (Stage 2b). Faster; loses splice discovery. |
| `--probe-uncovered` | `PROBE_UNCOVERED=1` | Probe intronic variants outside the gnomAD footprint. **~48 min GPU per proband**; the rarity gate cannot work there, and PM2 is withheld on those rows. |
| `--selftest` | — | Family-discovery self-test; exits. Needs no data. |
| `--selftest-cohort` | — | Cohort recurrent-artifact self-test; exits. Needs no data. |

### Custom gene list (genes of interest)

By default the panel is `g4e-2026.txt`. To restrict to a different gene set, pass a
genes-of-interest file (one gene symbol per line; `#` comments and blanks ignored) with
`-l`/`--list` — it is forwarded to both passes:

```bash
bash run_filtering.sh my_genes.txt        # full pipeline with the custom list
perl filtering_r.pl -l my_genes.txt       # filtering only
```

- The file may be **plain symbols** (Association/MOI/GDV columns are filled with `NA`) or the
  full 4-column g4e format (`gene⇥Association⇥MOI⇥GDV`, in which case those values are used).
- With `MOI = NA`, genes are treated as **dominant** for the rarity gate (`$FREQ_AD` = 0.01%).
  If a custom gene has recessive forms, supply its MOI (column 3 = `AR`) to get the
  recessive threshold, or relax `$FREQ_AD`.
- Outputs are **namespaced by panel** (`<proband>.<panel>.candidatos`, where `<panel>` is the
  panel-file basename with any trailing year dropped, e.g. `EPIC280-P.g4e.candidatos` vs `EPIC280-P.Hyperparathyroidism.candidatos`),
  so different gene lists produce **side-by-side** results instead of overwriting. Pangolin
  scratch is namespaced the same way but deleted after each run (see [Splice scoring](#splice-scoring-pangolin)).

### Single-variant lookup (`filtering_r.pl -v` / `--lookup`)

To **consult one (or a few) variants** and see *everything the pipeline can say about each* —
every predictor, ClinVar, gnomAD, PS1/PM5, the triage ACMG class, QC flags — in the **same
`.candidatos` format**, without any panel / rarity / consequence / evidence gating:

```bash
# GRCh38 genomic coordinates (100% offline) — chr-pos-ref-alt, or :/space separated
perl filtering_r.pl -v 'chr17-7675088-C-T'
perl filtering_r.pl -v '2:166073617:T:G'

# HGVS on an Ensembl transcript (resolved to coordinates via the Ensembl REST API)
perl filtering_r.pl -v 'ENST00000269305.9:c.524G>A'

# several at once: repeat -v
perl filtering_r.pl -v 'chr17-7675088-C-T' -v 'ENST00000269305.9:c.524G>A'

# options
perl filtering_r.pl -v 'chr17-7675088-C-T' -l my_genes.txt     # custom candidate-gene panel
perl filtering_r.pl -v 'chr17-7675088-C-T' --all-transcripts   # report every transcript, not just MANE
perl filtering_r.pl -v 'chr17-7675088-C-T' --keep-vcf          # keep the annotated VCF
```

Output: the **transposed, human-readable view only** — one `field <TAB> value` line per column
(the transposed view, not the TSV cohort table) — written to `Lookup.<tag>.<panel>.candidatos` and echoed to stdout.
For `-v`, `<tag>` is the variant id (`chr-pos-ref-alt`) for a single variant; for several, it is
`<first-id>_<N>` where **`N` is the number of *additional* variants** (`-v` count − 1). The tag is
sanitized to `[A-Za-z0-9._-]`, so two variants starting at `chr9-6644629-T-C` give
`Lookup.chr9-6644629-T-C_1.g4e.candidatos`.

A pre-annotated VCF can also be analyzed directly with `--lookup <file.germline.vep.vcf.gz>`, in
which case `<tag>` is the VCF's basename: `--lookup EPIC280-P.germline.vep.vcf.gz` writes
`Lookup.EPIC280-P.g4e.candidatos`. **Consult output always carries the `Lookup.` prefix**, so it
occupies a separate namespace from the `<proband>.<panel>.candidatos` tables a cohort run produces
and can never overwrite one.

- The variant(s) are built **sites-only** (no sample), so genotype columns (zygosity/GT/DP/GQ/AB)
  are blank and `inheritance = NA`; every annotation-derived field is still computed.
- **MANE-only** by default (use `--all-transcripts` to see all transcripts). A variant with no
  MANE annotation yields an empty file; re-run with `--all-transcripts`.
- `kept_by` lists whichever evidence arms fire (or `none`); off-panel genes get
  `Association/MOI/GDV = NA`, ACMG-SF genes get their condition + `GDV = Incidental`.
- **Coordinates resolve 100% offline.** **HGVS** requires transcript→genomic mapping, which VEP
  cannot do offline, so it is resolved via the **Ensembl REST API** (GRCh38) — only the variant
  notation is sent (a public variant string, **never patient data**); override the endpoint with
  `$ENSEMBL_REST`. The local cache is Ensembl (not RefSeq), so use `ENST…` HGVS, not `NM_…`.
  The HGVS path needs `curl` + `jq`; the coordinate path needs neither.
- **Splicing is scored too.** Because a single-variant consult should report *everything*,
  **`-v`** runs **Pangolin** on the variant inline (from the normalized annotated VCF) and
  fills `pangolin_score` + the splice rescue arm — no separate two-pass step needed. **`--lookup` on a
  pre-annotated VCF does not** — it reads an existing
  `<base>.<panel>.pangolin.tsv` from the **current working directory** if one is there, and otherwise
  leaves `pangolin_score` blank. It **degrades
  gracefully**: if the `pangolin` conda env or references are missing, or Pangolin fails, it warns
  and leaves `pangolin_score` blank rather than aborting. Pass **`--no-splice`** to skip it (faster,
  and avoids the conda/GPU dependency on air-gapped hosts).
- Mechanically: build a sites-only VCF → `vep_annotate.sh` (full annotation) → Pangolin (unless
  `--no-splice`) → report-everything readable output. The splice scratch and the temp dir are
  **always** removed afterward; **`--keep-vcf` keeps only the annotated VCF**
  (`Lookup.<tag>.germline.vep.vcf.gz` + index + `_summary.html`). The splice score itself is already
  in the `.candidatos` output, so nothing is lost.

### Forcing a proband

By default the proband is auto-detected from filenames (only `-P` samples are analyzed; `-M`/`-F`
are locked in as parents). To analyze a specific sample — e.g. the mother — override it by its
full base-name:

```bash
PROBAND="EPIC280-M" bash run_filtering.sh             # analyze the mother
PROBAND="EPIC280-P EPIC280-M" bash run_filtering.sh   # analyze both
perl filtering_r.pl --proband EPIC280-M              # filtering only
```

The forced sample must have a `<name>.germline.vep.vcf.gz`. Its parents are still derived from
the family prefix (`<family>-M` / `<family>-F`, stripping a trailing `-P`); if they are absent (as for a mother whose own parents aren't in
the dataset) the sample is analyzed as a **singleton** — `inheritance = NA`, no compound-het
*trans* phasing (HOM and `CompHet?` flags still apply from the sample's own genotypes). Each
proband writes its own `<name>.<panel>.candidatos`, so forcing one does not overwrite another.

### Configuration

Every external path lives in **`site.sh`**. Override any of them by exporting the variable, or by
creating an untracked **`site.env`** beside it (see [Setup 0.2](#02--tell-the-repo-where-your-data-lives)). Precedence: explicit export > `site.env` > the defaults below.

| Variable | Default | Used by |
|----------|---------|---------|
| `VEP` | `$HOME/ensembl-vep/vep` | annotation |
| `VEP_DATA` | `$HOME/vep_data` | annotation (`--dir_cache`) |
| `VEP_REFS` | `$HOME/vep_refs` | annotation (plugin data + custom VCFs) |
| `VEP_PLUGINS` | `$HOME/.vep/Plugins` | annotation (`--dir_plugins`, incl. `loftee/`) |
| `VEP_FORKS` | `4` | annotation |
| `CLINVAR_AA_DIR` | *(empty — PS1/PM5 skipped)* | filtering |
| `REF_FASTA` | `$PANGOLIN_FASTA` | filtering (homopolymer QC flag) |
| `CONDA_BASE` | `$HOME/miniconda3` | Pangolin |
| `PANGOLIN_ENV` | `pangolin` | Pangolin |
| `PANGOLIN_FASTA` | `$VEP_REFS/pangolin/GRCh38.primary_assembly.genome.fa` | Pangolin |
| `PANGOLIN_DB` | `$VEP_REFS/pangolin/gencode.v38.annotation.db` | Pangolin |
| `PERL5LIB_EXTRA` | `$VEP_PLUGINS/loftee:$HOME/perl5/lib/perl5` | annotation — extra Perl libs VEP/LOFTEE need |
| `HTSLIB_SO` | `$HOME/htslib/libhts.so` | annotation — `LD_PRELOAD`ed so LOFTEE's `Bio::DB::HTS` links (skipped silently if the file is absent) |

Those are the variables `site.sh` defines and to which the precedence rule above applies. The
remaining ones are read directly by the individual scripts. The four resource paths (`GNOMAD_VCF`,
`CLINVAR_VCF`, `CADD_SNV`, `CADD_INDEL`) can also be set in `site.env`, since it is sourced into the
same shell; `WORKDIR`, `PROBAND`, `ENSEMBL_REST` and the `KEEP_*` toggles belong on the command line:

| Variable | Default | Used by |
|----------|---------|---------|
| `WORKDIR` | the repo directory | `run_filtering.sh` — directory holding the `*.germline.vep.vcf.gz` inputs |
| `PROBAND` | *(empty)* | `run_filtering.sh` — space-separated sample(s) forwarded as `--proband` |
| `GNOMAD_VCF` | `$VEP_REFS/gnomAD_min/gnomAD.joint.v4.1.mane.all.vcf.gz` | `vep_annotate.sh` |
| `CLINVAR_VCF` | `$VEP_REFS/clinvar/clinvar.chr.vcf.gz` | `vep_annotate.sh` |
| `CADD_SNV` | `$VEP_REFS/CADD/whole_genome_SNVs.tsv.gz` | `vep_annotate.sh` |
| `CADD_INDEL` | `$VEP_REFS/CADD/gnomad.genomes.r4.0.indel.tsv.gz` | `vep_annotate.sh` |
| `ENSEMBL_REST` | `https://rest.ensembl.org` | `filtering_r.pl` — HGVS→coordinate recoding; point at a private mirror on an air-gapped host |
| `KEEP_AR_CARRIERS` | *(unset)* | `filtering_r.pl` — same as `--keep-ar-carriers` |
| `KEEP_COHORT_ARTIFACTS` | *(unset)* | `filtering_r.pl` — same as `--keep-cohort-artifacts` |

Filtering thresholds (`$FREQ_AD`, `$FREQ_AR`, `$CADD_MIN`, `$REVEL_MIN`, `$AM_MIN`,
`$SPLICE_MIN`, and the cohort-artifact `$COHORT_MIN` / `$COHORT_MAX_FRAC`) are
edited directly in `filtering_r.pl`. `--keep-ar-carriers` / `KEEP_AR_CARRIERS` and
`--keep-cohort-artifacts` / `KEEP_COHORT_ARTIFACTS` toggle the two drop rules.

---

## Notes & limitations

- The rarity gate binds only variants with **no established ClinVar classification**; a P/LP call
  with ≥1 review star is exempt (Stage 1), so a founder allele above the ceiling is still reported.
- `$FREQ_AR` = 1% is deliberately permissive (sensitivity); solitary het carriers in pure
  recessive genes are dropped by default (surfaced only via `--keep-ar-carriers`), so the permissive threshold
  matters only for variants that pair up. `gnomAD_nhomalt` surfaces high-homozygote variants
  for quick triage. Tighten if noisy.
- `REVEL ≥ 0.644` matches the ClinGen PP3 calibration. The AlphaMissense rescue uses a
  *score* threshold (`am_score ≥ 0.792`), not the categorical `am_class`.
- Compound-het *trans* confirmation needs a full trio; duos report `CompHet?`.
- De-novo calls rely on parent VCF genotypes; a parental no-call (uncovered site) can
  masquerade as de novo — verify against parental depth before reporting.
- The **cohort recurrent-artifact filter** drops variants that are simultaneously cohort-recurrent and
  gnomAD-absent (Stage 3). It is deliberately conservative (both conditions, high 25% threshold) so
  founder / bottleneck alleles — which carry a gnomAD footprint — are preserved. For a strongly
  founder-enriched cohort where a *private* founder allele could plausibly reach 25% while being
  gnomAD-absent, review the per-run drop log or run with `--keep-cohort-artifacts`.
- **Without Pangolin, splice evidence is simply absent, never assumed.** Running `filtering_r.pl`
  directly (see [Setup 0.7](#07--pangolin-splice-scoring-optional-but-recommended)) leaves
  `pangolin_score` blank; the splice rescue arm and BP7 both stay silent rather than defaulting either
  way. Synonymous variants therefore remain unclassified on splicing instead of being labelled benign
  on no evidence — use `run_filtering.sh` when that distinction matters.
- **Known-open triage limitations** (deliberate, not defects — they change *class*, not *coverage*):
  - **PM2 is now Supporting** (`$PM2_STRENGTH`), per ClinGen SVI 2020. ⚠️ **Read this before
    interpreting classes:** the downgrade interacts badly with *categorical* combining. ACMG 2015
    has no "PVS1 + 1 supporting" pathway, so a gnomAD-absent nonsense or frameshift variant in a
    disease gene (`PVS1,PM2_Supporting`) now lands on **VUS**, not Likely pathogenic. Under the
    ClinGen/Tavtigian **Bayesian points** framework — the framework PM2_Supporting was calibrated
    in — the same variant scores 8+1 = 9 points and stays Likely pathogenic. The two are not
    interchangeable, and this pipeline combines categorically. Set `$PM2_STRENGTH = 'moderate'` to
    restore the 2015 reading. Moving to a points-based combiner is the real fix and is not done.
    Separately, a **missing** gnomAD annotation is still coerced to `AC=0`, so PM2 cannot distinguish
    "gnomAD never saw this allele" from "the position is outside the MANE-restricted resource".
  - **PS2 vs PM6.** A trio `DN` earns **PS2** (Strong) when the *proband's* genotype is clean, but the
    discriminator carries no information about parentage or parental coverage — which is what actually
    separates PS2 from PM6 — and the same rows are stamped `DN_unconfirmed`. Trios are assumed
    confirmed; treat PS2 rows as PM6 unless relatedness and parental coverage were verified.
  - **BA1/BS1 use global joint AF**, not popmax or filtering AF, and have no ClinGen exception list, so
    a known pathogenic founder allele that survives the Stage-1 ClinVar exemption can auto-classify as
    Benign. `BA1` also fires at exactly 5% (ACMG specifies *>* 5%).
- This is a **triage tool to feed manual curation**, not an automated classifier.

## Changelog — 2026-08

Correctness fixes from a full audit. Everything here changes **which variants reach the curator** or
**what class they carry**, so tables produced before this point are not comparable.

### Second pass (same day)

| Change | Effect |
|---|---|
| **PM2 → Supporting** (`$PM2_STRENGTH`) | ClinGen SVI 2020. ⚠️ Under categorical combining this drops gnomAD-absent LoF variants from Likely pathogenic to VUS — see the limitations section |
| **Splice discovery probes** (Stage 2b) | Deep-intronic (≤300 bp) and synonymous variants are scored by Pangolin and rescued if they disrupt splicing. Previously unreachable — the splice arm could only upgrade, never discover |
| **ACMG-SF genes reach Pangolin** | Incidental rows can finally carry a real `pangolin_score` and earn BP7 / the splice rescue |
| Probes require gnomAD coverage (`AN > 0`) | The custom gnomAD VCF is MANE-restricted; probing uncovered introns would rescue variants with no frequency filter and a manufactured PM2. Measured: only ~8 of ~22,700 in-window intronic variants per proband are covered — so discovery is real but resource-limited. `--probe-uncovered` widens it |
| PM2 withheld on gnomAD-uncovered probe rows | `AN=0` outside the resource means "not looked at", not "unobserved" |
| `--no-splice-discovery` | Opt out of the probe set |

### First pass

| Change | Effect |
|---|---|
| Plain `XL` MOI now matches both MOI predicates | 72 panel genes (CDKL5, MECP2, ARX, IQSEC2, PCDH19, DDX3X, ATRX…) had matched *neither*, so they got no HOM/CompHet flag, no `AR_hom` rescue, and the strict dominant AF ceiling |
| Haploid `GT=1` parsed as `hem` | Male non-PAR chrX/chrY and chrM calls had no zygosity; a **hemizygous parent was invisible**, so inherited X-linked variants were reported de novo (and fed PS2/PM6) |
| One-MANE-row collapse keyed on (variant, **gene**) | A panel candidate could silently delete a reportable ACMG-SF incidental at the same position |
| **PP5** requires ≥1 review star | A single 0-star ClinVar submission (~16% of the P/LP corpus) was lifting LP calls to Pathogenic |
| **PM4** suppressed when PVS1 fired | VEP compound consequence terms yielded two ACMG lines for one protein-terminus effect |
| `bcftools norm -f` left-alignment | Non-minimal indels missed the gnomAD/ClinVar `type=exact` join → spurious `AC=0` → spurious PM2 |
| Cohort artifact thresholds (`$COHORT_MIN` 10→5, new `$COHORT_MIN_CARRIERS`=3, samples as the unit) | The filter had **never fired on a real clinical run**; KMT2C/SYNE1 artifacts reached every delivered table |
| New `clinvar_conflict` flag | Surfaces P/LP auto-classes that carry a hard benign line (e.g. `PVS1,PM2,BP6` on a ClinVar-Benign variant) |
| `flags` replaces `recessive_flag` + `qc_flag` | One column instead of two, four apart |
| `batch.<panel>.candidatos` | Batch-level table, `sample` column first |
| Five HGNC aliases for the constraint file + startup coverage line | PP2 was silently dead for renamed genes (GBA1, BMAL1, AFG2A, AFG2B, BLTP1) |
| Reference files resolve from the repo dir | `$WORKDIR` now works without symlinking the repo into every run directory |
| Empty Pangolin score map aborts the run | A silent Pangolin failure produced a complete-looking table with the splice arm dead |

## Data privacy

Patient VCFs and all run outputs (`*.candidatos`, `*.pangolin*`, `*_summary.html`, logs) are PHI.

`.gitignore` is an **allow-list**: it ignores `*` and then un-ignores only pipeline code and
non-patient reference config, so patient data cannot be committed by accident. **Never remove the
leading `*` rule**; to track a new file add an explicit `!<file>` exception.

That design is what makes it safe for **this repository to be public** — only code and public
reference data are tracked. Your own paths belong in `site.env` (untracked), never in a tracked
file. If you fork this for a deployment where run outputs might land inside the working tree, verify
`git status` shows nothing patient-derived before your first push.
