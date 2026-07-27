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
| `vep_annotate.sh` | Split multiallelics + Ensembl VEP (LOFTEE/REVEL/AlphaMissense/EVE/CADD + custom gnomAD v4.1 & ClinVar). → `*.germline.vep.vcf.gz` |
| `filtering_r.pl` | The filtering algorithm (Perl, no modules). Reads annotated VCF, applies gates, writes `<proband>.<panel>.candidatos`. **All thresholds are constants at the top of this file.** Also the single-variant consult entry point (`-v`): coords/HGVS in → annotate → transposed readable view (`Lookup.<coords>.<panel>.candidatos`), gates bypassed. |
| `parse_pangolin.pl` | Reduce Pangolin output to per-variant `max(\|Δ\|)` splice score. |
| `run_filtering.sh` | End-to-end driver: emit candidates → Pangolin → final filtering → cleanup. |
| `site.sh` | Every external path (VEP/plugins/Pangolin/ClinVar-AA). Sourced by `vep_annotate.sh` + `run_filtering.sh`; exports the env vars `filtering_r.pl` reads. Override via untracked `site.env`. **No absolute personal paths in tracked code.** |
| `test/test_filtering.sh` | Regression test (synthetic, ~5s, no VEP/GPU). **Run after touching `filtering_r.pl`** — a broken gate yields a plausible table, not an error. |
| `run_wgs.sh` / `run_4probands.sh` | One-off batch drivers (WGS 2-of-4 merge; 4 DRAGEN singletons). Idempotent, log to `logs/`. |
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

# Filtering only (no Pangolin)
perl filtering_r.pl                   # default panel
perl filtering_r.pl -l my_genes.txt   # custom panel (-l/--list is the ONLY way; no positional arg)
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
  singleton cohorts work without `--proband`. (Historically this silently produced 0 candidatos.)
- **Two-pass design:** if the Pangolin score map is absent, pass 1 emits the candidate list and
  stops; once scores exist, pass 2 writes the final table. `run_filtering.sh` does both.
- **Annotations resolved by name** from the CSQ header — no hard-coded column indices.
- **Outputs are namespaced by panel** (`<proband>.<panel>.candidatos`), so different gene lists
  produce side-by-side results instead of overwriting.
- **Pangolin scratch is deleted after a successful pass 2** (recomputed every run — cheap, only
  the few-hundred structural-pass variants are scored). A **failed** run leaves intermediates in
  place for debugging.
- **gnomAD AC=0/AN=0 means "absent from gnomAD"** (sites-only VCF), NOT an uncallable region —
  never cite AN=0 as artifact evidence.
- **Cohort recurrent-artifact filter (internal panel-of-normals):** when a real cohort is
  auto-analyzed together (≥ `$COHORT_MIN`=10 probands), a candidate carried by ≥ `$COHORT_MAX_FRAC`=25%
  of the cohort **and** absent from gnomAD (AF < `$COHORT_ART_FREQ`=0.01%) is dropped as a systematic
  technical artifact (paralog/low-complexity mismapping — e.g. the recurrent SYNE1/KMT2C sites).
  **Both** conditions are required: recurrence alone would hit founder alleles, but a true population
  bottleneck / under-represented ancestry allele carries a gnomAD footprint (gnomAD's large
  Admixed-American sample), so requiring gnomAD-absence keeps it founder-safe. Every drop is logged
  (`cohort_artifact drop:` lines + a per-proband count). **OFF** for single-variant (`-v`/`--lookup`),
  forced/single proband (`--proband`), and any run of < 10 probands — there's no cohort to compare
  against. **Override:** `--keep-cohort-artifacts` (or env `KEEP_COHORT_ARTIFACTS=1`) keeps them
  tagged `qc_flag=cohort_artifact` instead of dropping (for founder-enriched cohorts: review the drop
  log — a genuinely private founder allele would surface there). Self-test: `perl filtering_r.pl --selftest-cohort`.
- **One MANE row per variant:** a variant hitting >1 MANE transcript (MANE Select + MANE Plus Clinical,
  or overlapping gene models — e.g. MUTYH) is collapsed to a single row (prefer panel-primary → MANE
  Select → most evidence arms). `--lookup` still reports every annotation.
- **Dual-inheritance genes** (panel MOI has **both** AD and AR, e.g. `AD, AR`) are treated as **dominant**
  for the carrier logic: a solitary het passes through as a normal candidate (`recessive_flag` empty),
  while a genuine HOM/comp-het still gets the recessive flag. Prevents dropping a dominant-acting variant
  (LoF etc.) just because the gene also has a recessive mechanism. Only **pure** AR/XLR genes use the
  carrier path. The HOM/comp-het flag pass runs for **recessive-capable genes only**, so a purely dominant
  gene with two independent hets is no longer mislabeled `CompHet?`.
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
- **ACMG output is triage-grade**, not a final clinical call (PM1/PP2 not assessed; PVS1 doesn't
  verify gene mechanism/NMD; PS1/PM5 rely on ClinVar AA matching). PM5 also fires for a
  single-codon in-frame deletion when a P/LP missense exists at the deleted residue (curatorial
  extension of PM5 beyond missense; tagged `(in-frame del)` in the `clinvar_aa` column).
- **De-novo can't be confirmed here** — parent VCFs are variant-only (no reference depth), so `DN`
  rows are flagged `DN_unconfirmed`. True confirmation needs parental gVCFs/BAMs.
- **Secondary findings (ACMG SF v3.2)** are always scanned independent of the panel, flagged
  `GDV=Incidental`, and carry distinct consent/reporting obligations.
- **`-l/--list <genes>` overrides the default g4e-2026 panel** — the only way to set the panel
  (no positional argument; a stray positional is now a hard error). Works in normal runs and in
  `-v` consults alike; sets the Association/MOI/GDV columns and the output `<panel>` tag.
- **Lookup mode** (`filtering_r.pl -v <variant>` (repeatable), or `--lookup <annotated.vcf.gz>`)
  bypasses all gates to report variants in full; all changes are `$LOOKUP`-guarded so normal runs
  are unaffected. `-v` builds a sites-only VCF, runs `vep_annotate.sh`, **and runs Pangolin inline**
  (so `pangolin_score` + the splice rescue arm work for a single variant too) — removing the
  annotated VCF + splice scratch afterward unless `--keep-vcf`. Pangolin degrades gracefully (warns,
  leaves score blank) if the env/refs are absent; **`--no-splice`** skips it. Coordinates are fully
  offline; **HGVS resolution calls the Ensembl REST API** (only the variant string, never patient
  data; needs `curl`+`jq`) and needs `ENST…` ids (the cache is Ensembl, not RefSeq).

## ACMG criteria evaluated (triage-grade)

Auto-assigned per variant by `acmg_classify`; combined per categorical ACMG 2015. Full thresholds
+ PP3/BP4 calibration are in `README.md`.

| Criterion | What triggers it |
|---|---|
| **PVS1** | LoF: LOFTEE = HC, or truncating consequence with LOFTEE ≠ LC |
| **PS1** | Same AA change is ClinVar P/LP (≥1★) |
| **PS2** / **PM6** | De novo confirmed in trio (PS2) / assumed de novo, unconfirmed or duo (PM6) |
| **PM2** | Absent or singleton in gnomAD (AC ≤ 1) |
| **PM4** | Protein length change (in-frame indel / `stop_lost`) |
| **PM5** | Different change — or **single-codon in-frame deletion** — at a residue with P/LP missense (≥1★) |
| **PP2** | Missense in a missense-constrained gene: gnomAD v4.1.1 `mis.oe < 0.6` (MANE, outliers excluded; `gnomad-mis-constraint.txt`). Counts **independently of PP3** (separate ACMG lines), but **suppressed when BP4 fires** (no gene-level pathogenic support for a benign-predicted variant). Guard in `acmg_classify` is `!$bp4`; drop it to fire even alongside BP4. |
| **PP3** / **BP4** | Computational, graded Supp/Mod/Strong (AlphaMissense primary, REVEL fallback) |
| **PP5** / **BP6** | This variant reported P/LP (PP5) or B/LB ≥1★ (BP6) in ClinVar |
| **BA1** / **BS1** / **BS2** | gnomAD AF ≥ 5% / ≥ 1% / ≥ 10 homozygotes |
| **BP7** | Synonymous with no predicted splice impact (Pangolin < 0.2) |

**Not evaluated (manual curation):** PS3/BS3, PS4, PM1, PM3, PP1/BS4, PP2, PP4, BP1/BP2/BP3/BP5.

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
