# CLAUDE.md

Germline variant-filtering pipeline for **clinical candidate triage**: annotate a patient
germline VCF, then keep rare, gene-panel variants with damaging evidence (missense
pathogenicity, high CADD, splicing impact, ClinVar, LoF), label them with inheritance/recessive
context, auto-assign a triage ACMG class, and emit a curation-ready table. **This is a triage
tool that feeds manual curation — not an automated final classifier.**

See `README.md` for the full algorithm reference; this file is the working quick-start.

## ⚠️ Patient data / PHI — read first

- Patient VCFs and **all run outputs** (`*.candidatos`, `*.pangolin*`, `*_summary.html`, logs,
  `*.vcf.gz`/`*.cram`) are **PHI and must never be committed**.
- `.gitignore` is an **allow-list**: it ignores `*` and then un-ignores only code + non-patient
  reference config. **Never remove the leading `*` rule.** To track a new code/config file, add
  an explicit `!<file>` exception.
- **This repo is PUBLIC** — that is safe *because* the allow-list works, not despite it. Only code
  and public reference data are tracked. Before adding any `!<file>` exception, confirm the file
  contains no patient-derived content. Local paths go in the untracked `site.env`, never in a
  tracked file.

## Layout

| File | Purpose |
|------|---------|
| `vep_annotate.sh` | Split multiallelics **+ left-align** (`bcftools norm -m-any -f $REF_FASTA`) + Ensembl VEP (LOFTEE/REVEL/AlphaMissense/EVE/CADD + custom gnomAD v4.1 & ClinVar). → `*.germline.vep.vcf.gz`. Left-alignment is **required**: both custom sources join with `type=exact`, so a non-minimal indel misses gnomAD, reports AC=0, and collects a spurious PM2. Warns loudly if `REF_FASTA` is unset. |
| `filtering_r.pl` | The filtering algorithm (Perl, no modules). Reads annotated VCF, applies gates, writes `<proband>.<panel>.candidatos`. **All thresholds are constants at the top of this file.** Also the single-variant consult entry point (`-v`): coords/HGVS in → annotate → transposed readable view (`Lookup.<coords>.<panel>.candidatos`), gates bypassed. |
| `parse_pangolin.pl` | Reduce Pangolin output to per-variant `max(\|Δ\|)` splice score. |
| `run_filtering.sh` | End-to-end driver: emit candidates → Pangolin → final filtering → cleanup. |
| `site.sh` | Every external path (VEP/plugins/Pangolin/ClinVar-AA). Sourced by `vep_annotate.sh` + `run_filtering.sh`; exports the env vars `filtering_r.pl` reads. Override via untracked `site.env`. **No absolute personal paths in tracked code.** |
| `test/test_filtering.sh` | Regression test (synthetic, ~5s, no VEP/GPU). **Run after touching `filtering_r.pl`** — a broken gate yields a plausible table, not an error. |
| `run_wgs.sh` / `run_4probands.sh` | One-off batch drivers (WGS 2-of-4 merge; 4 DRAGEN singletons). Idempotent, log to `logs/`. |
| **Outputs** | `<proband>.<panel>.candidatos` per sample **and** `batch.<panel>.candidatos` — every proband's rows in one table, prefixed with a `sample` column, columns 2..N identical to the per-sample file. |
| `g4e-2026.txt` | Default gene panel (`gene⇥Association⇥MOI⇥GDV`). Source: Genes4Epilepsy v2026-03 (bahlolab), 1078 genes. GDV (disease + MONDO) carried over from the prior g4e-2025 for the 93 genes that had one; NO_GDV otherwise (v2026-03 has no GDV column). |
| `typevar.txt` | Consequence whitelist. |
| `mane-plus-clinical-names.txt` | MANE Select + MANE Plus Clinical transcript allow-list. |
| `acmg_sf_v3.2.txt` | 81 ACMG SF v3.2 secondary-findings genes — **always** scanned. |
| `gnomad-mis-constraint.txt` | gnomAD v4.1.1 missense constraint per MANE gene (`gene⇥mis.oe⇥mis.z⇥flags`) — drives ACMG **PP2**. |

## Common commands

```bash
# Annotate each family member — name outputs with the role suffix:
#   -P proband, -M mother, -F father  (filename drives family auto-discovery)
bash vep_annotate.sh EPIC280.raw.vcf.gz  EPIC280-P.germline.vep.vcf.gz
bash vep_annotate.sh EPIC280M.raw.vcf.gz EPIC280-M.germline.vep.vcf.gz

# Full pipeline (emit → Pangolin GPU scoring → final) over all *.germline.vep.vcf.gz
bash run_filtering.sh                 # default g4e-2026 panel → EPIC280-P.g4e.candidatos
bash run_filtering.sh my_genes.txt    # custom genes-of-interest list (forwarded to both passes)

# Run a batch from its own directory. Put every *.germline.vep.vcf.gz for the batch in
# ONE directory so the cohort artifact filter can see the whole cohort — reference files
# resolve from the repo, so the run directory needs nothing but the VCFs.
WORKDIR=/path/to/batch4-run bash run_filtering.sh

# Filtering only (no Pangolin)
perl filtering_r.pl                   # default panel
perl filtering_r.pl -l my_genes.txt   # custom panel (-l/--list is the ONLY way; no positional arg)
                                      # -l sets the PANEL; consults are -v / --lookup
perl filtering_r.pl --selftest        # built-in family-discovery self-test

# Force a specific sample as proband (overrides filename auto-discovery)
PROBAND="EPIC280-M" bash run_filtering.sh
perl filtering_r.pl --proband EPIC280-M

# Consult a SINGLE variant (coords offline; HGVS via Ensembl REST) — report everything
perl filtering_r.pl -v 'chr17-7675088-C-T'              # GRCh38 chr-pos-ref-alt
perl filtering_r.pl -v 'ENST00000269305.9:c.524G>A'     # Ensembl-transcript HGVS
perl filtering_r.pl -v 'chr17-7675088-C-T' --all-transcripts
perl filtering_r.pl -v 'chr17-7675088-C-T' -l my_genes.txt   # override the g4e panel (Association/MOI/GDV)
```

## Key conventions & gotchas

- **Role-suffix filename convention** drives family auto-discovery: `<FAMILY>-P/-M/-F`. Globs
  `*.germline.vep.vcf.gz`, groups by `<FAMILY>` prefix, analyzes each `-P` as proband. A name
  not ending in `-P/-M/-F` is ignored by discovery (still usable via `--proband`).
  **Singleton fallback:** if NO input has a `-P/-M/-F` suffix (e.g. `EPIGEN01..20`), discovery finds
  0 probands → it auto-analyzes **every** sample as a singleton proband (prints a `NOTE:`). So plainly-named
  singleton cohorts work without `--proband`.
- **Two-pass design:** if the Pangolin score map is absent, pass 1 emits the candidate list and
  stops; once scores exist, pass 2 writes the final table. `run_filtering.sh` does both.
- **Annotations resolved by name** from the CSQ header — no hard-coded column indices.
- **ClinVar exemption on the Stage-1 gates:** a panel-gene variant that is ClinVar P/LP with **≥1
  review star** passes Stage 1 regardless of the rarity ceiling and the consequence whitelist (MANE +
  panel still required). Without it, founder alleles above the AF ceiling and pathogenic non-coding /
  synonymous variants are discarded before any evidence arm runs. It earns `BS1`/`BA1` as usual, so
  the frequency tension is visible instead of silent. Same exemption the ACMG-SF path already had.
- **`start`/`end` span the REF allele** (`end = start + len(ref) - 1`), so indels report their true
  footprint and rows round-trip into the 5-field `chr-start-end-ref-alt` consult form.
- **Outputs are namespaced by panel** (`<proband>.<panel>.candidatos`), so different gene lists
  produce side-by-side results instead of overwriting.
- **Pangolin scratch is deleted after a successful pass 2** (recomputed every run — cheap, only
  the few-hundred structural-pass variants are scored). A **failed** run leaves intermediates in
  place for debugging.
- **gnomAD AC=0/AN=0 means "absent from gnomAD"** (sites-only VCF), NOT an uncallable region —
  never cite AN=0 as artifact evidence.
- **Cohort recurrent-artifact filter (internal panel-of-normals):** when a real cohort is
  auto-analyzed together (≥ `$COHORT_MIN`=**5 samples**, ≥2 probands), a candidate carried by
  ≥ `$COHORT_MIN_CARRIERS`=**3 samples** *and* ≥ `$COHORT_MAX_FRAC`=25% of the cohort **and**
  absent from gnomAD (joint **AC = 0**) is dropped as a systematic technical artifact
  (paralog/low-complexity mismapping — e.g. the recurrent SYNE1/KMT2C sites).
  **All three** conditions are required. gnomAD-absence is what keeps it founder-safe: a true
  population bottleneck allele carries a gnomAD footprint (gnomAD's large Admixed-American sample).
  The **absolute carrier floor** is what keeps small batches safe — at N=8 the 25% fraction alone
  means 2 carriers, and two unrelated probands sharing a private allele is an ordinary
  founder/relatedness event. The floor binds on small batches, the fraction binds on large cohorts.
  **Thresholds revised 2026-08:** `$COHORT_MIN` was 10 probands, which no internal batch reaches
  (they run 6–9), so the filter had never fired on a real clinical run and the KMT2C/SYNE1 artifacts
  reached every delivered table. Denominator and eligibility are both counted in **samples** now
  (previously eligibility counted probands while the denominator counted all VCFs, diluting the
  fraction ~3× on trio cohorts). Every drop is logged (`cohort_artifact drop:` lines + a per-proband
  count). **OFF** for single-variant (`-v`/`--lookup`), forced/single proband (`--proband`), and runs
  of < 5 samples. **Override:** `--keep-cohort-artifacts` (or env `KEEP_COHORT_ARTIFACTS=1`) keeps
  them tagged `flags=cohort_artifact` instead of dropping (for founder-enriched cohorts: review the
  drop log — a genuinely private founder allele would surface there).
  Self-test: `perl filtering_r.pl --selftest-cohort`.
  **The filter only sees samples in ONE run** — a batch split across per-sample directories has a
  cohort of 1. Put the whole batch's VCFs in one `$WORKDIR`.
- **One MANE row per variant:** a variant hitting >1 MANE transcript (MANE Select + MANE Plus Clinical,
  or overlapping gene models — e.g. MUTYH) is collapsed to a single row (prefer panel-primary → MANE
  Select → most evidence arms). `--lookup` still reports every annotation.
- **Dual-inheritance genes** (panel MOI has **both** AD and AR, e.g. `AD, AR`) are treated as **dominant**
  for the carrier logic: a solitary het passes through as a normal candidate (`recessive_flag` empty),
  while a genuine HOM/comp-het still gets the recessive flag. Prevents dropping a dominant-acting variant
  (LoF etc.) just because the gene also has a recessive mechanism. Only **pure** AR/XLR genes use the
  carrier path. The HOM/comp-het flag pass runs for **recessive-capable genes only**, so a purely dominant
  gene with two independent hets is never labelled `CompHet?`.
- **One `flags` column.** The recessive verdict and the QC/artifact flags are a single
  `;`-separated `flags` column (recessive verdict first), not the former `recessive_flag` + `qc_flag`
  pair. Values: `HOM` / `HEM` / `CompHet(trans)` / `CompHet?` / `carrier-only`, then `lowDP`, `lowGQ`,
  `AB_het`, `AB_hom`, `homopolymer`, `GT_rescued`, `inh_lowqual`, `DN_unconfirmed`,
  `cohort_artifact`, `clinvar_conflict`.
- **The recessive verdict is decided per gene but written per row.** The gene-level verdict governs
  which rows survive the carrier drop; the label on each row describes that row — `HOM` only on a
  homozygous row, `HEM` only on a hemizygous one, `CompHet*` on the hets that constitute it. A het in
  a gene that is hom for a *different* variant is left blank, so the flag never contradicts the row's
  own `zygosity` column.
- **Haploid genotypes (`GT=1`) are hemizygous, not no-calls.** DRAGEN emits single-allele GTs for
  non-PAR chrX/chrY in males and for chrM. They report `zygosity=hem`, are treated as biallelic-
  equivalent by the recessive logic (`HEM` flag, `AR_hem` rescue arm, `AB_hom` QC, never dropped as a
  carrier), and — critically — a **hemizygous parent now counts as a carrier**, so an inherited
  X-linked variant is no longer misreported as de novo.
- **Plain `XL` in the panel MOI is dual-inheritance** (both dominant- and recessive-capable), the same
  treatment `AD, AR` genes get. 72 of the 1078 g4e-2026 genes use plain `XL` — CDKL5, MECP2, ARX,
  IQSEC2, PCDH19, DDX3X, ATRX, SLC6A8, FLNA. They previously matched neither MOI predicate, so they
  got no HOM/CompHet flag, no `AR_hom` rescue, and the strict dominant AF ceiling.
  `XLR`/`XLD` keep their specific meanings.
- **Recessive carrier drop (DEFAULT):** a solitary het in a **pure** AR/XLR gene that is not biallelic is
  **dropped** — a single het can't explain recessive disease; carrier states are clinical noise. **True
  comp-hets are unaffected:** a gene with ≥2 gate-passing hets is a biallelic `CompHet` and kept. Same drop
  covers recessive ACMG-SF genes. **Opt-in:** `--keep-ar-carriers` / `KEEP_AR_CARRIERS=1` surfaces the
  **strong** such carriers for a targeted 2nd-hit hunt (carrier-only tier: ClinVar P/LP ≥1★, HC-LoF, or ≥2
  strong predictors AM≥0.906/CADD≥28.1/EVE-path/REVEL≥0.773, not Benign/LB; flagged
  `recessive_flag=carrier-only`). Note a common SNP (high gnomAD AF) is NOT a valid 2nd hit.
- **AR_hom rescue arm:** a **homozygous, protein-altering** (missense/inframe/stop_lost/start_lost) rare
  MANE variant in a recessive (AR/XLR) panel gene **with AB > 0.75** is kept even with no predictor/ClinVar
  support (`kept_by=AR_hom`). *General rationale* (not case-tuned): a biallelic genotype in a recessive
  disease gene is itself pathogenicity evidence under recessive inheritance, independent of missense
  predictors (calibrated mostly on dominant/het effects). **Coding-only** (intronic/splice → Pangolin arm;
  LoF → LoF arm) + **AB>0.75** (guards false-hom) keep it specific — else it floods on benign homozygous
  polypyrimidine/intron variants. `BS1/BS2/BA1` still flag benign-leaning ones.
- **Contradictions are flagged, not resolved.** `acmg_class` reaching Pathogenic/Likely pathogenic
  while a hard benign line fired (`BP6`/`BS1`/`BS2`/`BA1`) is tagged `flags=clinvar_conflict`. The
  categorical `Conflicting` verdict needs *both* sides to reach a 2-tier threshold independently, so
  a single benign criterion never blocks a pathogenic call — e.g. `PVS1,PM2,BP6` reads as
  Likely_pathogenic on a variant ClinVar calls Benign with review stars. This is a triage tool, so
  the class is left alone and the tension is made visible for the curator. `BP4` is deliberately
  excluded: a computational prediction disagreeing with PVS1/PM2 is routine, not a contradiction.
- **Splice discovery (`$SPLICE_PROBE`).** `typevar.txt` has no bare `intron_variant` or
  `synonymous_variant`, so those variants used to be dropped at Stage 1 and never scored — the
  Pangolin arm could only *upgrade* an already-whitelisted splice-region variant, never *discover*
  one. Intronic (within `$INTRON_MAX_DIST`=300 bp of an exon boundary, read from the HGVSc offset)
  and synonymous variants in panel/ACMG-SF genes are now **probed**: sent to Pangolin, and kept
  **only** if the splice arm fires (≥ `$SPLICE_MIN`). A probe that scores low simply disappears, so
  the probe set widens what can be *found* without widening the table. Probes additionally require
  **gnomAD coverage (`AN > 0`)** — see the caveat below. Disable with `--no-splice-discovery` or
  `NO_SPLICE_DISCOVERY=1`. Cost: Pangolin runs ~8 variants/s on one GPU.
- **⚠️ The gnomAD resource is MANE-restricted.** `gnomAD.joint.v4.1.mane.all.vcf.gz` covers MANE
  transcripts and flanks, *not* deep intronic sequence. An uncovered position yields `AC=""`/`AN=""`,
  which this code coerces to 0 — so `$freq` computes as 0, passes **every** rarity ceiling, and PM2
  fires on what is really an annotation gap. That is why probes require `AN > 0`: it bounds the probe
  set to where the resource can actually answer the question, and stops discovery from manufacturing
  PM2. **True deep-intronic discovery (beyond the MANE footprint, e.g. CFTR c.3718-2477C>T) needs the
  custom gnomAD VCF rebuilt with genome-wide coverage.**
- **ACMG-SF genes reach Pangolin too**, so an incidental row can carry a real `pangolin_score` and
  earn BP7 or the splice rescue. Previously only panel candidates were scored.
- **ACMG output is triage-grade**, not a final clinical call (PM1 is regional hotspot evidence, not
  a curated functional-domain assessment; PP2 is gene-level constraint only; PVS1 doesn't verify gene mechanism/NMD; PS1/PM5 rely
  on ClinVar AA matching). **BP7 requires a real Pangolin score** — an unscored synonymous variant is
  left unclassified on splicing rather than called benign. PM5 also fires for a
  single-codon in-frame deletion when a P/LP missense exists at the deleted residue (curatorial
  extension of PM5 beyond missense; tagged `(in-frame del)` in the `clinvar_aa` column).
- **De-novo can't be confirmed here** — parent VCFs are variant-only (no reference depth), so `DN`
  rows are flagged `DN_unconfirmed`. True confirmation needs parental gVCFs/BAMs.
- **Secondary findings (ACMG SF v3.2)** are always scanned independent of the panel, flagged
  `GDV=Incidental`, and carry distinct consent/reporting obligations.
- **`-l/--list <genes>` overrides the default g4e-2026 panel** — the only way to set the panel
  (no positional argument; a stray positional is a hard error). Works in normal runs and in
  `-v` consults alike; sets the Association/MOI/GDV columns and the output `<panel>` tag.
- **Lookup mode** (`filtering_r.pl -v <variant>` (repeatable), or `--lookup <annotated.vcf.gz>`)
  bypasses all gates to report variants in full; all changes are `$LOOKUP`-guarded so normal runs
  are unaffected. `-v` builds a sites-only VCF, runs `vep_annotate.sh`, **and runs Pangolin inline**
  (so `pangolin_score` + the splice rescue arm work for a single variant too) — removing the
  annotated VCF afterward unless `--keep-vcf` (splice scratch is ALWAYS removed; `--lookup` on a
  pre-annotated VCF does NOT run Pangolin at all). Consult output ALWAYS carries the `Lookup.` prefix
  (`Lookup.<tag>.<panel>.candidatos`, `<tag>` = variant id for `-v`, VCF basename for `--lookup`), so
  it is namespaced apart from cohort `<proband>.<panel>.candidatos` tables and cannot overwrite one.
  Pangolin degrades gracefully (warns,
  leaves score blank) if the env/refs are absent; **`--no-splice`** skips it. Coordinates are fully
  offline; **HGVS resolution calls the Ensembl REST API** (only the variant string, never patient
  data; needs `curl`+`jq`) and needs `ENST…` ids (the cache is Ensembl, not RefSeq).

## ACMG criteria evaluated (triage-grade)

Auto-assigned per variant by `acmg_classify`; combined per categorical ACMG 2015. Full thresholds
+ PP3/BP4 calibration are in `README.md`.

| Criterion | What triggers it |
|---|---|
| **PVS1** | LoF: LOFTEE = HC, or truncating consequence with LOFTEE ≠ LC |
| **PS1** | A **different** variant giving the same AA change is ClinVar P/LP (≥1★); the variant's own ClinVar record is excluded (the AA resource is indexed by source `chr-pos-ref-alt`), so a self-match cannot count one submission as both PS1 and PP5 |
| **PS2** / **PM6** | De novo in a **full trio** (PS2 — relatedness assumed confirmed) / assumed de novo, poor proband genotype or duo-ambiguous (PM6). PS2 is structurally unreachable without both parents: `inheritance` is only ever `DN` when both are present, singletons get `NA` and duos get `DN/IM`–`DN/IF` |
| **PM2** | Absent or singleton in gnomAD (AC ≤ 1). Counted at **Moderate** (ACMG 2015) — `$PM2_STRENGTH`. ClinGen SVI 2020 recommends Supporting, and the knob implements it, but ⚠️ **only adopt that with a points-based combiner**: ACMG 2015 has no "PVS1 + 1 supporting" rule, so under categorical combining Supporting demotes every gnomAD-absent LoF variant to VUS (measured: KCNT1, CUX2, RELN, HCN2). That is framework-mixing, not conservatism |
| **PM4** | Protein length change (in-frame indel / `stop_lost`) — **not counted when PVS1 fired**, so one protein-terminus effect cannot yield two ACMG lines via a compound consequence term |
| **PM5** | Different change — or **single-codon in-frame deletion** — at a residue with P/LP missense (≥1★) |
| **PM1** | Missense / in-frame indel inside a **PERv1** pathogenic-variant-enriched region naming **that same gene** (Pérez-Palma, *Genome Res* 2020 — PM1 is the paper's stated application). Two arms: **`PERv1_direct`** (enrichment computed on that gene) is graded by its own fold enrichment at the published Tavtigian-2018 calibration — **≥ 18.7 → `PM1_Strong`**, else Moderate; **`PERv1_paralog`** (enrichment computed across the paralog family alignment and assigned to every member — the paper's headline arm, 1,252 genes vs 215, and the one its held-out de novo test validated) is **capped at Moderate**, because transferring a family's evidence onto one member costs a tier. On overlap, direct outranks paralog. Not restricted by BP4 (unlike PP2). The winning arm and region go to `flags` as `PM1:<arm>/<PER>/<aa range>/OR=…`. Needs `$VEP_REFS/PER/PERv1.GRCh38.MANE.bed.gz`; absent → PM1 never fires. |
| **PP2** | Missense in a missense-constrained gene: gnomAD v4.1.1 `mis.oe < 0.6` (MANE, outliers excluded; `gnomad-mis-constraint.txt`; five HGNC renames resolved via `%GENE_ALIAS`, and per-panel coverage is printed at startup). Counts **independently of PP3** (separate ACMG lines), but **suppressed when BP4 fires** (no gene-level pathogenic support for a benign-predicted variant). Guard in `acmg_classify` is `!$bp4`; drop it to fire even alongside BP4. |
| **PP3** / **BP4** | Computational, graded Supp/Mod/Strong (AlphaMissense primary, REVEL fallback) |
| **PP5** / **BP6** | This variant reported P/LP (PP5) or B/LB (BP6) in ClinVar — **both require ≥1 review star** |
| **BA1** / **BS1** / **BS2** | gnomAD AF ≥ 5% / ≥ 1% / ≥ 10 homozygotes |
| **BP7** | Synonymous with no predicted splice impact (Pangolin < 0.2) |

**Not evaluated (manual curation):** PS3/BS3, PS4, PM3, PP1/BS4, PP4, BP1/BP2/BP3/BP5.

**PS1/PM5 are strictly same-gene and carry no paralog transfer** — that is PM1's job, and only via
the explicitly-labelled `PERv1_paralog` arm.

**PM1 and PS1/PM5 may both fire at one residue.** PER regional evidence is validated independently of
the individual ClinVar submissions behind it, so it is not suppressed — but ClinGen SVI cautions
against reusing one piece of evidence twice, so the co-occurrence is surfaced as
`flags=PM1_with_ps1` / `PM1_with_pm5` for the curator to rule on. **PS1/PM5 remain strictly
same-gene** (the ClinVar AA index is keyed `GENE\tAApos\tRefAA`); there is no paralog transfer
anywhere in this pipeline.

## Dependencies

- **Annotation:** Ensembl VEP (offline GRCh38 cache) + LOFTEE/REVEL/AlphaMissense/EVE/CADD
  plugins & data + custom gnomAD v4.1 and ClinVar VCFs (chr-prefixed, bgzipped, tabixed). Paths
  come from `site.sh` (`VEP`, `VEP_DATA`, `VEP_REFS`, `VEP_PLUGINS`) — override in an untracked
  `site.env`, never by editing tracked code. README "Setup" has the from-scratch install.
  Also needs `bcftools`.
- **Splice scoring:** conda env `pangolin` with PyTorch (GPU) + Pangolin, a chr-named GRCh38
  primary-assembly FASTA, and `gencode.v38.annotation.db`. Override via `CONDA_BASE`,
  `PANGOLIN_ENV`, `PANGOLIN_FASTA`, `PANGOLIN_DB`.
- **Filtering:** system Perl only (no modules). PS1/PM5 additionally read
  `clinvar.MANE_missense.{PLP,BLB}.tsv` from `$CLINVAR_AA_DIR` (**unset by default** — set it in
  `site.env`); if absent, filtering still runs and PS1/PM5 are skipped with a warning.
- **Input compatibility:** reads both single-source VCFs (e.g. DRAGEN) and the Sarek
  union-consensus output of `consensus.sh` (picks up `GT_SOURCE`/`NCALLERS`/`CONF` tags when present).
