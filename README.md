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
| `filtering_r.pl` | The filtering algorithm. Reads the annotated VCF, applies gates, writes `<proband>.<panel>.candidatos`. Also the **single-variant consult** entry point (`-v`/`-l`): annotate one or a few variants (coords or HGVS) from scratch and report everything, gates bypassed. |
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

- Splits multiallelic sites (`bcftools norm -m-any`) so downstream `chr-pos-ref-alt`
  keys are unambiguous.
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
row per variant** (a variant hitting >1 MANE transcript — MANE Select + MANE Plus Clinical,
or overlapping gene models — keeps a single row: panel-primary first, then MANE Select, then
the most evidence arms). `--lookup` consults still report every annotation.

### Stage 1 — Structural gates (ALL required, AND)

| Gate | Source | Rule |
|------|--------|------|
| MANE transcript | `mane-plus-clinical-names.txt` | CSQ `Feature` ∈ MANE set |
| Consequence | `typevar.txt` | consequence split on `&`; kept if **any** sub-term is whitelisted |
| Gene panel | `g4e-2026.txt` (default) or a custom genes-of-interest file | CSQ `SYMBOL` ∈ panel |
| Rarity (MOI-aware) | gnomAD joint AC/AN | AF = AC/AN×100 ≤ threshold: **dominant `$FREQ_AD`=0.01%**, **recessive `$FREQ_AR`=1.0%** (MOI contains AR/XLR) |

### Stage 2 — Inclusion / rescue gate (at least ONE, OR)

A surviving variant must trip **one or more** of these. Each is independent; a `kept_by`
column records which fired.

| Arm | Threshold |
|-----|-----------|
| CADD | `$CADD_MIN` = 25.3 |
| AlphaMissense | `am_score` ≥ `$AM_MIN` = 0.792 (ClinGen PP3) |
| EVE pathogenic | `eve_class` is Pathogenic |
| REVEL | `$REVEL_MIN` = 0.644 (ClinGen PP3) |
| Pangolin (splice) | `$SPLICE_MIN` = 0.5 (max \|Δscore\|) |
| ClinVar P/LP | `ClinVar_CLNSIG` Pathogenic/Likely_pathogenic (excludes Conflicting & Benign) |
| PS1 / PM5 | ClinVar amino-acid match (≥1★): **PS1** = same AA change is P/LP, **PM5** = a different change at the same residue is P/LP. A **single-codon in-frame deletion** of the residue also triggers PM5 (a different protein change at the same P/LP residue; tagged `(in-frame del)`). Rescues the variant even when CADD/AM/REVEL miss it; the `clinvar_aa` column carries the detail (and any `(conflicting)` flag). |
| LoF | LOFTEE `LoF=HC`, or a high-impact truncating consequence (frameshift / stop_gained / splice_donor / splice_acceptor / start_lost) unless LOFTEE downgraded it to `LC`. Covers truncating indels that CADD (SNV-only) and the missense predictors miss. |
| AR_hom | **Homozygous, protein-altering** (missense / inframe / stop_lost / start_lost) variant in a **recessive (AR/XLR) panel gene**, with clean allele balance (**AB > 0.75**). Rescued even without in-silico/ClinVar support. *Rationale:* under recessive inheritance a **biallelic (homozygous) genotype in a disease gene is itself pathogenicity evidence**, independent of missense predictors — which are calibrated largely on dominant/heterozygous effects and can miss true recessive alleles. **Restricted to coding changes** (so it doesn't flood on benign homozygous intronic/polypyrimidine variants — those use the Pangolin arm; truncating LoF uses the LoF arm) **and to clean homozygous calls** (AB > 0.75 guards against false-hom artifacts). Already rare (AR freq gate), MANE, in-panel by this point; `BS1`/`BS2`/`BA1` still flag benign-leaning ones. |

All thresholds are single constants at the top of `filtering_r.pl`.

### Genotype-aware annotation (not gates)

- **Zygosity / GT / DP / GQ / AB** are read from the proband `FORMAT`/sample column.
- **Inheritance** uses *parental genotype* (carrier = non-ref GT, not mere site presence):
  `IB / IM / IF / DN` in a full trio; duo-ambiguous `DN/IF` (mother-only) or `DN/IM`
  (father-only); `NA` for a singleton.
- **`recessive_flag`** per gene (computed for **recessive-capable genes only** — AR/XLR/dual;
  a purely dominant gene never gets one, so two independent hets are not mislabeled comp-het):
  `HOM` (homozygous), `CompHet(trans)` (≥2 het variants phaseable to opposite parents — trio
  only), `CompHet?` (≥2 het, unphaseable — e.g. duo/singleton), or `carrier-only` (see below).
- **Dual-inheritance genes** (panel MOI lists **both** AD and AR, e.g. `AD, AR`) are treated as
  **dominant** for the carrier logic: a **solitary het passes through as a normal candidate**
  (`recessive_flag` empty), while a genuine `HOM`/comp-het still gets the recessive flag. This
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
    `recessive_flag=carrier-only`. Use it when you suspect the exome missed a second allele (deep-intronic,
    CNV, regulatory) in a specific gene/case.
    ⚠️ A **common** SNP (high gnomAD AF) is not a valid second hit — an apparent "comp-het" of one rare
    plus one common variant is really a **solitary carrier** of the rare allele, not a biallelic genotype.

### Stage 3 — Cohort recurrent-artifact filter (internal panel-of-normals)

When a **real cohort** is auto-analyzed together (≥ `$COHORT_MIN` = 10 probands), the pipeline
builds an internal *panel of normals*: a one-pass genotype tally (no CSQ, so it is cheap) over all
input samples counting, per `chr-pos-ref-alt`, how many carry the ALT and their zygosity breakdown.
A candidate is then **dropped** when it is **both**:

1. **cohort-recurrent** — carried by ≥ `$COHORT_MAX_FRAC` (25%) of the cohort, and
2. **gnomAD-absent** — gnomAD AF < `$COHORT_ART_FREQ` (0.01%, i.e. ~1 in 10,000).

This targets systematic technical artifacts — reference/mapping errors in paralog-rich or
low-complexity genes (e.g. the recurrent `SYNE1` / `KMT2C` sites that appear, gnomAD-absent, in a
large fraction of a cohort). **Both conditions are required, and that conjunction is the
discriminator:** recurrence alone would remove true **founder / bottleneck** alleles, but a real
founder allele frequent enough to reach a quarter of the cohort leaves a footprint in gnomAD's large
Admixed-American sample — so requiring gnomAD-absence keeps the filter founder-safe (important for
under-represented ancestries). Every drop is logged (`cohort_artifact drop:` lines with carrier
count + hom/het breakdown, and a per-proband total).

- **OFF** for single-variant (`-v` / `--lookup`), forced/single proband (`--proband`), and any run of
  fewer than `$COHORT_MIN` probands — a per-variant or per-proband consult has no cohort to compare
  against.
- **`--keep-cohort-artifacts`** (or env `KEEP_COHORT_ARTIFACTS=1`) keeps them, tagged
  `qc_flag=cohort_artifact`, instead of dropping — for a founder-enriched cohort, review the drop log
  (a genuinely private founder allele would surface there and can be kept this way).
- Self-test (synthetic cohort, no reference files needed): `perl filtering_r.pl --selftest-cohort`.
- Thresholds are the `$COHORT_MIN` / `$COHORT_MAX_FRAC` / `$COHORT_ART_FREQ` constants at the top of
  `filtering_r.pl`.

A per-proband **run summary** prints counts (read / multiallelic-skipped / structural-pass
/ candidates / cohort-artifacts dropped) and breakdowns by `kept_by` and inheritance.

### Output columns (`<proband>.<panel>.candidatos`, TSV)

`chr, start, end, ref, alt, gene, strand, consequence, hgvs,
revel, eve_class, eve_score, cadd, am_class, am_score, pangolin_score,
clinvar_sig, clinvar_stars, clinvar_disease, clinvar_aa, loftee,
gnomAD_ac, gnomAD_an, gnomAD_af, gnomAD_nhomalt, gnomAD_filter,
zygosity, GT, DP, GQ, AB, GT_SOURCE, NCALLERS, CONF,
inheritance, recessive_flag, kept_by,
acmg_class, acmg_criteria, qc_flag, Association, MOI, GDV`

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
  | **PS1** | Same amino-acid change is ClinVar P/LP (≥1★) | ClinVar MANE-missense |
  | **PS2** | De novo **confirmed** in a trio (`inheritance=DN`, clean genotype) | parental GT |
  | **PM2** | Absent or singleton in gnomAD (AC ≤ 1) | gnomAD v4.1 |
  | **PM4** | Protein length change (in-frame indel / `stop_lost`) | consequence |
  | **PM5** | Different change — **or a single-codon in-frame deletion** — at a residue carrying a P/LP missense (≥1★) | ClinVar MANE-missense |
  | **PM6** | **Assumed** de novo (DN unconfirmed / duo) | parental GT |
  | **PP2** | **Missense** in a gene with low benign-missense variation — gnomAD v4.1.1 missense constraint `mis.oe < 0.6` (MANE; constraint outliers excluded). Counts **independently of PP3** (both are legitimate separate ACMG lines — gene-level intolerance vs variant-level prediction), but **suppressed when BP4 fires** (a benign-predicted variant gets no gene-level pathogenic support). | gnomAD v4.1.1 constraint |
  | **PP3** | Computational damaging, graded Supporting/Moderate/Strong (see below) | AlphaMissense / REVEL |
  | **PP5** | This variant is reported pathogenic in ClinVar | ClinVar |

  **Benign**

  | Criterion | What triggers it | Source |
  |---|---|---|
  | **BA1** | gnomAD AF ≥ 5% | gnomAD v4.1 |
  | **BS1** | gnomAD AF ≥ 1% (and < 5%) | gnomAD v4.1 |
  | **BS2** | ≥ 10 homozygotes in gnomAD | gnomAD v4.1 |
  | **BP4** | Computational benign, graded (see below) | AlphaMissense / REVEL |
  | **BP6** | This variant is reported benign in ClinVar (≥1★) | ClinVar |
  | **BP7** | Synonymous with no predicted splice impact (Pangolin < 0.2) | Pangolin |

  **Not evaluated (manual curation only):** PS3/BS3 (functional), PS4 (case-control), PM1
  (hotspot/domain), PM3 (in trans), PP1/BS4 (segregation), PP2 (missense-constrained gene),
  PP4 (phenotype specificity), BP1/BP2/BP3/BP5.

  **PS1/PM5** use the ClinVar MANE-missense resource (`clinvar.MANE_missense.{PLP,BLB}.tsv`), matched on
  gene + protein residue + amino-acid change, requiring **≥1 review star**; a match also reported B/LB is
  tagged **`(conflicting)`** (still counted — flag for manual review) and detailed in the `clinvar_aa` column.
  **PP3/BP4 come from a single calibrated predictor**, graded **Supporting/Moderate/Strong**:
  **AlphaMissense** primary ([Bergquist 2025](https://doi.org/10.1016/j.gim.2025.101402): PP3
  supp ≥0.792 / mod ≥0.906 / strong ≥0.990; BP4 supp ≤0.169 / mod ≤0.099), **REVEL** fallback
  ([Pejaver 2022](https://doi.org/10.1016/j.ajhg.2022.10.013): PP3 supp ≥0.644 / mod ≥0.773 /
  strong ≥0.932; BP4 supp ≤0.290 / mod ≤0.183 / strong ≤0.016) — with a **REVEL direction-conflict
  veto**, mapped to the 2015 tiers (BP4_Moderate → supporting-benign, since 2015 has no benign-Moderate).
  **Not a final clinical call**: PM1/PP2 not assessed; PVS1 doesn't verify gene mechanism/NMD; PS1/PM5
  rely on ClinVar AA matching (no independent re-curation, and PS1 may overlap PP5 for the same variant).
- **`qc_flag`** — artifact/confidence warnings: `lowDP` (<`$QC_MIN_DP`), `lowGQ` (<`$QC_MIN_GQ`),
  `AB_het`/`AB_hom` (skewed allele balance), `homopolymer` (indel in a ≥5 homopolymer — error-prone),
  `GT_rescued` (genotype borrowed from a non-DeepVariant caller via `consensus.sh`; no VAF),
  `inh_lowqual` (carrying-parent genotype is weak), `DN_unconfirmed`, `cohort_artifact` (recurrent
  gnomAD-absent cohort artifact, present only under `--keep-cohort-artifacts` — otherwise dropped).
- **De-novo confidence [#6]:** parent VCFs here are *variant-only* (no reference depth at non-variant
  sites), so de-novo cannot be confirmed from parental coverage — `DN` rows are flagged
  `DN_unconfirmed`. Inherited rows instead get `inh_lowqual` when the parental call is low quality.
  True de-novo confirmation needs parental gVCFs/BAMs.

---

## Secondary findings (ACMG SF v3.2)

The **81 ACMG SF v3.2 genes** (`acmg_sf_v3.2.txt`) are **always** scanned, independent of the
candidate `--genes` panel, with a **stricter** gate than candidates. Findings are written into
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

### Single-variant lookup (`filtering_r.pl -v` / `-l`)

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
`<tag>` is the variant id (`chr-pos-ref-alt`) for a single `-v`, else `<first-id>+<N>`. (A
pre-annotated VCF can still be analyzed directly with `--lookup <file.germline.vep.vcf.gz>`.)

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
  `-v`/`--lookup` runs **Pangolin** on the variant inline (from the normalized annotated VCF) and
  fills `pangolin_score` + the splice rescue arm — no separate two-pass step needed. It **degrades
  gracefully**: if the `pangolin` conda env or references are missing, or Pangolin fails, it warns
  and leaves `pangolin_score` blank rather than aborting. Pass **`--no-splice`** to skip it (faster,
  and avoids the conda/GPU dependency on air-gapped hosts).
- Mechanically: build a sites-only VCF → `vep_annotate.sh` (full annotation) → Pangolin (unless
  `--no-splice`) → report-everything readable output. The annotated VCF + all splice scratch are
  removed afterward unless `--keep-vcf`.

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

Filtering thresholds (`$FREQ_AD`, `$FREQ_AR`, `$CADD_MIN`, `$REVEL_MIN`, `$AM_MIN`,
`$SPLICE_MIN`, and the cohort-artifact `$COHORT_MIN` / `$COHORT_MAX_FRAC` / `$COHORT_ART_FREQ`) are
edited directly in `filtering_r.pl`. `--keep-ar-carriers` / `KEEP_AR_CARRIERS` and
`--keep-cohort-artifacts` / `KEEP_COHORT_ARTIFACTS` toggle the two drop rules.

---

## Notes & limitations

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
- This is a **triage tool to feed manual curation**, not an automated classifier.

## Data privacy

Patient VCFs and all run outputs (`*.candidatos`, `*.pangolin*`, `*_summary.html`, logs) are PHI.

`.gitignore` is an **allow-list**: it ignores `*` and then un-ignores only pipeline code and
non-patient reference config, so patient data cannot be committed by accident. **Never remove the
leading `*` rule**; to track a new file add an explicit `!<file>` exception.

That design is what makes it safe for **this repository to be public** — only code and public
reference data are tracked. Your own paths belong in `site.env` (untracked), never in a tracked
file. If you fork this for a deployment where run outputs might land inside the working tree, verify
`git status` shows nothing patient-derived before your first push.
