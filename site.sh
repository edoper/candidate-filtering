# Shared site settings for every script in this repo. Sourced (not executed), silent.
#   . "$(dirname "$(readlink -f "$0")")/site.sh"
#
# This pipeline depends on several large external resources (VEP + cache + plugin data,
# custom gnomAD/ClinVar VCFs, Pangolin, the ClinVar amino-acid tables). Nothing about
# WHERE you put them is baked into the code — every path below can be overridden by
# exporting it beforehand, or by creating an untracked `site.env` next to this file:
#
#   # site.env
#   VEP_REFS=/data/vep_refs
#   VEP_DATA=/data/vep_cache
#   CLINVAR_AA_DIR=/data/clinvar-aa
#   PANGOLIN_ENV=pangolin
#
# See README section 0 for how to obtain each resource. The defaults are the layout
# this pipeline was developed against — a starting point, not a requirement.

CF_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export CF_REPO

# Untracked per-site overrides. Precedence: explicit export > site.env > defaults.
if [ -f "$CF_REPO/site.env" ]; then
    _keys=(); _vals=()
    for _v in VEP VEP_DATA VEP_REFS VEP_PLUGINS PERL5LIB_EXTRA HTSLIB_SO \
              CLINVAR_AA_DIR CONDA_BASE PANGOLIN_ENV PANGOLIN_FASTA PANGOLIN_DB REF_FASTA VEP_FORKS; do
        [ -n "${!_v:-}" ] && { _keys+=("$_v"); _vals+=("${!_v}"); }
    done
    . "$CF_REPO/site.env"
    for _i in "${!_keys[@]}"; do printf -v "${_keys[$_i]}" '%s' "${_vals[$_i]}"; done
    unset _v _i _keys _vals
fi

# --- Ensembl VEP -----------------------------------------------------------
: "${VEP:=$HOME/ensembl-vep/vep}"          # the vep executable
: "${VEP_DATA:=$HOME/vep_data}"            # --dir_cache (offline GRCh38 cache)
: "${VEP_REFS:=$HOME/vep_refs}"            # plugin data + custom VCFs live under here
: "${VEP_PLUGINS:=$HOME/.vep/Plugins}"     # --dir_plugins (incl. loftee/)
: "${VEP_FORKS:=4}"                        # --fork
export VEP VEP_DATA VEP_REFS VEP_PLUGINS VEP_FORKS

# Extra Perl libs VEP/LOFTEE need, and the htslib that LOFTEE's Bio::DB::HTS links.
: "${PERL5LIB_EXTRA:=$VEP_PLUGINS/loftee:$HOME/perl5/lib/perl5}"
: "${HTSLIB_SO:=$HOME/htslib/libhts.so}"
export PERL5LIB_EXTRA HTSLIB_SO

# --- Splice scoring (Pangolin) ---------------------------------------------
: "${CONDA_BASE:=$HOME/miniconda3}"
: "${PANGOLIN_ENV:=pangolin}"
: "${PANGOLIN_FASTA:=$VEP_REFS/pangolin/GRCh38.primary_assembly.genome.fa}"
: "${PANGOLIN_DB:=$VEP_REFS/pangolin/gencode.v38.annotation.db}"
export CONDA_BASE PANGOLIN_ENV PANGOLIN_FASTA PANGOLIN_DB

# filtering_r.pl uses an indexed FASTA for the homopolymer QC flag; the Pangolin
# reference is already chr-named GRCh38, so reuse it unless told otherwise.
: "${REF_FASTA:=$PANGOLIN_FASTA}"
export REF_FASTA

# --- ClinVar amino-acid tables for ACMG PS1/PM5 (optional) ------------------
# Directory holding clinvar.MANE_missense.{PLP,BLB}.tsv. Empty by default: without it
# filtering still runs and PS1/PM5 are simply skipped (with a warning). See README §0.6.
: "${CLINVAR_AA_DIR:=}"
export CLINVAR_AA_DIR
