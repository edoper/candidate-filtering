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
| `filtering_r.pl` | The filtering algorithm. Reads the annotated VCF, applies gates, writes `<proband>.<panel>.candidatos`. |
| `parse_pangolin.pl` | Convert Pangolin output into a per-variant splice-score map (`<proband>.<panel>.pangolin.tsv`). |
| `run_filtering.sh` | End-to-end driver: emit candidates → score with Pangolin → final filtering. |
| `g4e-2025.txt` | Gene panel: `gene⇥Association⇥MOI⇥GDV`. Restricts output to panel genes; supplies MOI. |
| `typevar.txt` | Consequence whitelist (atomic terms; matched per `&`-separated sub-term). |
| `mane-plus-clinical-names.txt` | MANE Select + MANE Plus Clinical transcript IDs; only these transcripts are considered. |

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
   Pangolin (GPU, de novo)  ──►  <proband>.<panel>.pangolin.csv  ──►  <proband>.<panel>.pangolin.tsv
                                   │
   Pass 2  filtering_r.pl  ──►  <proband>.<panel>.candidatos   (final, curation-ready)
```

`filtering_r.pl` auto-discovers families by **filename**: it globs
`*.germline.vep.vcf.gz`, takes the sample name as the filename minus that suffix, treats a
sample ending in `M`/`F` as a parent (when the base name also exists), and analyzes the
rest as probands (`EPIC280` = proband, `EPIC280M` = mother, `EPIC280F` = father). It is a
two-pass design — if the Pangolin score map is missing it emits the candidate list and
stops; once scores exist it produces the final table.

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

Logic runs **per transcript annotation** of each variant (so one variant can yield one
row per qualifying MANE transcript).

### Stage 1 — Structural gates (ALL required, AND)

| Gate | Source | Rule |
|------|--------|------|
| MANE transcript | `mane-plus-clinical-names.txt` | CSQ `Feature` ∈ MANE set |
| Consequence | `typevar.txt` | consequence split on `&`; kept if **any** sub-term is whitelisted |
| Gene panel | `g4e-2025.txt` (default) or a custom genes-of-interest file | CSQ `SYMBOL` ∈ panel |
| Rarity (MOI-aware) | gnomAD joint AC/AN | AF = AC/AN×100 ≤ threshold: **dominant `$FREQ_AD`=0.01%**, **recessive `$FREQ_AR`=1.0%** (MOI contains AR/XLR) |

### Stage 2 — Inclusion / rescue gate (at least ONE, OR)

A surviving variant must trip **one or more** of these. Each is independent; a `kept_by`
column records which fired.

| Arm | Threshold |
|-----|-----------|
| CADD | `$CADD_MIN` = 22 |
| AlphaMissense pathogenic | `am_class` is (likely_)pathogenic |
| EVE pathogenic | `eve_class` is Pathogenic |
| REVEL | `$REVEL_MIN` = 0.5 *(permissive; ClinGen PP3 ≈ 0.644)* |
| Pangolin (splice) | `$SPLICE_MIN` = 0.5 (max \|Δscore\|) |
| ClinVar P/LP | `ClinVar_CLNSIG` Pathogenic/Likely_pathogenic (excludes Conflicting & Benign) |
| LoF | LOFTEE `LoF=HC`, or a high-impact truncating consequence (frameshift / stop_gained / splice_donor / splice_acceptor / start_lost) unless LOFTEE downgraded it to `LC`. Covers truncating indels that CADD (SNV-only) and the missense predictors miss. |

All thresholds are single constants at the top of `filtering_r.pl`.

### Genotype-aware annotation (not gates)

- **Zygosity / GT / DP / GQ / AB** are read from the proband `FORMAT`/sample column.
- **Inheritance** uses *parental genotype* (carrier = non-ref GT, not mere site presence):
  `IB / IM / IF / DN` in a full trio; duo-ambiguous `DN/IF` (mother-only) or `DN/IM`
  (father-only); `NA` for a singleton.
- **`recessive_flag`** per gene: `HOM` (homozygous), `CompHet(trans)` (≥2 het variants
  phaseable to opposite parents — trio only), or `CompHet?` (≥2 het, unphaseable — e.g. duo).

A per-proband **run summary** prints counts (read / multiallelic-skipped / structural-pass
/ candidates) and breakdowns by `kept_by` and inheritance.

### Output columns (`<proband>.<panel>.candidatos`, TSV)

`chr, start, end, ref, alt, gene, strand, transcript, consequence, hgvs.c, hgvs.p, tpos,
revel, eve_class, eve_score, cadd, am_class, am_score, pangolin_score,
clinvar_sig, clinvar_stars, clinvar_disease, loftee, loftee_filter, loftee_flags,
gnomAD_ac, gnomAD_an, gnomAD_af, gnomAD_nhomalt, gnomAD_filter,
zygosity, GT, DP, GQ, AB, inheritance, recessive_flag, kept_by, Association, MOI, GDV`

---

## Splice scoring (Pangolin)

There is no official Pangolin VEP plugin, so splice scores are computed **de novo** by
standalone [Pangolin](https://github.com/tkzeng/Pangolin) on GPU and merged by position.
Only the proband's **structural-pass** variants are scored (a few hundred), not the whole
VCF. `parse_pangolin.pl` reduces each variant to `max(|increase|, |decrease|)`.

> Pangolin scores are cached per panel (`<proband>.<panel>.pangolin.tsv`) and regenerated
> automatically when the candidate set changes (the input-CSV checksum differs) — e.g. after
> editing the structural filters, the gene list, or re-annotating.

---

## Setup

**Annotation** requires Ensembl VEP (offline GRCh38 cache) with the LOFTEE / REVEL /
AlphaMissense / EVE / CADD plugins and data, plus the custom gnomAD and ClinVar VCFs
(chr-prefixed, bgzipped, tabixed). Paths are set near the top of `vep_annotate.sh`.

**Splice scoring** requires a Python/conda env with PyTorch (GPU) and Pangolin:

```bash
conda create -y -n pangolin -c conda-forge python=3.10 pip
conda activate pangolin
pip install torch --index-url https://download.pytorch.org/whl/cu121
pip install pyvcf3 gffutils biopython pandas pyfastx "setuptools<81"
pip install git+https://github.com/tkzeng/Pangolin.git
```

Plus a chr-named GRCh38 primary-assembly FASTA and the GENCODE annotation DB
(`gencode.v38.annotation.db`) — see `run_filtering.sh` for the expected paths.

**Filtering** needs only system Perl (no modules).

---

## Usage

```bash
# 1) Annotate each family member (proband + parents)
bash vep_annotate.sh EPIC280.germline.vcf.gz   EPIC280.germline.vep.vcf.gz
bash vep_annotate.sh EPIC280M.germline.vcf.gz  EPIC280M.germline.vep.vcf.gz

# 2) Run the full filtering pipeline (emit → Pangolin → final)
bash run_filtering.sh
#    → EPIC280.g4e-2025.candidatos
```

### Custom gene list (genes of interest)

By default the panel is `g4e-2025.txt`. To restrict to a different gene set, pass a
genes-of-interest file (one gene symbol per line; `#` comments and blanks ignored) — it is
forwarded to both passes:

```bash
bash run_filtering.sh my_genes.txt        # full pipeline with the custom list
perl filtering_r.pl  my_genes.txt         # filtering only
```

- The file may be **plain symbols** (Association/MOI/GDV columns are filled with `NA`) or the
  full 4-column g4e format (`gene⇥Association⇥MOI⇥GDV`, in which case those values are used).
- With `MOI = NA`, genes are treated as **dominant** for the rarity gate (`$FREQ_AD` = 0.01%).
  If a custom gene has recessive forms, supply its MOI (column 3 = `AR`) to get the
  recessive threshold, or relax `$FREQ_AD`.
- Outputs are **namespaced by panel** (`<proband>.<panel>.candidatos`, where `<panel>` is the
  panel-file basename, e.g. `EPIC280.g4e-2025.candidatos` vs `EPIC280.Hyperparathyroidism.candidatos`),
  so different gene lists produce **side-by-side** results instead of overwriting. Pangolin score
  caches are namespaced the same way and recomputed only when the candidate set changes (tracked
  by an input-CSV checksum).

### Forcing a proband

By default the proband is auto-detected from filenames (a sample ending in `M`/`F` is locked
in as a parent and never analyzed on its own). To analyze a specific sample — e.g. the
mother — override it:

```bash
PROBAND="EPIC280M" bash run_filtering.sh              # analyze the mother
PROBAND="EPIC280 EPIC280M" bash run_filtering.sh      # analyze both
perl filtering_r.pl --proband EPIC280M               # filtering only
```

The forced sample must have a `<name>.germline.vep.vcf.gz`. Its parents are still derived by
name (`<name>M` / `<name>F`); if they are absent (as for a mother whose own parents aren't in
the dataset) the sample is analyzed as a **singleton** — `inheritance = NA`, no compound-het
*trans* phasing (HOM and `CompHet?` flags still apply from the sample's own genotypes). Each
proband writes its own `<name>.<panel>.candidatos`, so forcing one does not overwrite another.

### Configuration (env overrides for `run_filtering.sh`)

| Variable | Default |
|----------|---------|
| `CONDA_BASE` | `$HOME/miniconda3` |
| `PANGOLIN_ENV` | `pangolin` |
| `PANGOLIN_FASTA` | `$HOME/vep_refs/pangolin/GRCh38.primary_assembly.genome.fa` |
| `PANGOLIN_DB` | `$HOME/vep_refs/pangolin/gencode.v38.annotation.db` |

Filtering thresholds (`$FREQ_AD`, `$FREQ_AR`, `$CADD_MIN`, `$REVEL_MIN`, `$SPLICE_MIN`)
are edited directly in `filtering_r.pl`.

---

## Notes & limitations

- `$FREQ_AR` = 1% is deliberately permissive (sensitivity); it can admit variants with
  many gnomAD homozygotes — the `gnomAD_nhomalt` column surfaces these for quick triage.
  Tighten if noisy.
- `REVEL ≥ 0.5` is sensitivity-tuned, below the ClinGen PP3 calibration (~0.644).
- Compound-het *trans* confirmation needs a full trio; duos report `CompHet?`.
- De-novo calls rely on parent VCF genotypes; a parental no-call (uncovered site) can
  masquerade as de novo — verify against parental depth before reporting.
- This is a **triage tool to feed manual curation**, not an automated classifier.

## Data privacy

Patient VCFs and all run outputs (`*.candidatos`, `*.pangolin*`, `*_summary.html`, logs)
are PHI and are excluded by the allow-list `.gitignore`. Keep any repository hosting this
code **private**.
