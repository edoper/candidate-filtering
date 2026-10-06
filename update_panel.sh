#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────
#  Rebuild the default gene panel (g4e.txt) from its two public sources.
#  Run it after every Genes4Epilepsy release (March and September).
#
#    Genes4Epilepsy (bahlolab)  -> gene list, Association (phenotype), MOI
#    ClinGen gene-validity      -> GDV column + removal of weakly supported genes
#    panel_overrides.tsv        -> curated MOI fixes / KEEP / REMOVE (tracked, with reasons)
#
#  Usage:
#    update_panel.sh                       # newest G4E release + today's ClinGen export
#    update_panel.sh --g4e FILE.tsv --clingen FILE.csv [--out g4e.txt] [--overrides FILE]
#
#  GDV column = CLASS|disease|MONDO|MOI|date of ONE ClinGen assertion per gene, chosen as:
#    epilepsy/neurodevelopmental-relevant > MOI compatible with the panel > classification rank
#    > most recent; falls back to the gene's highest-ranked assertion of any disease.
#    NOT_CURATED when ClinGen has no curation for the gene.
#
#  Removal rule: a gene is dropped when it has a RELEVANT, MOI-compatible ClinGen assertion
#  and every MOI-compatible assertion it has (any disease) is Limited / Disputed / Refuted /
#  No Known Disease Relationship. Curations of a different, narrower entity (e.g. "Leigh
#  syndrome" for POLG) are not treated as relevant, so they never remove a gene on their own.
#  Removed genes are listed in the output header with the assertion that removed them.
# ──────────────────────────────────────────────────────────────────────────
set -euo pipefail
REPO="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
G4E="" CLINGEN="" OUT="$REPO/g4e.txt" OVR="$REPO/panel_overrides.tsv"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --g4e) G4E="$2"; shift 2;;
        --clingen) CLINGEN="$2"; shift 2;;
        --out) OUT="$2"; shift 2;;
        --overrides) OVR="$2"; shift 2;;
        -h|--help) sed -n '2,25p' "$0"; exit 0;;
        *) echo "ERROR: unknown argument $1" >&2; exit 1;;
    esac
done
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

G4E_SRC="$G4E"
if [[ -z "$G4E" ]]; then
    REL=$(curl -fsS https://api.github.com/repos/bahlolab/Genes4Epilepsy/contents/ \
          | grep -oE 'EpilepsyGenes_v[0-9]{4}-[0-9]{2}\.tsv' | sort -u | tail -1)
    [[ -n "$REL" ]] || { echo "ERROR: could not list Genes4Epilepsy releases" >&2; exit 1; }
    G4E_SRC="https://github.com/bahlolab/Genes4Epilepsy/blob/main/$REL"
    G4E="$TMP/$REL"
    curl -fsS -o "$G4E" "https://raw.githubusercontent.com/bahlolab/Genes4Epilepsy/main/$REL"
fi
CLINGEN_SRC="$CLINGEN"
if [[ -z "$CLINGEN" ]]; then
    CLINGEN_SRC="https://search.clinicalgenome.org/kb/gene-validity/download"
    CLINGEN="$TMP/clingen.csv"
    curl -fsS -o "$CLINGEN" "$CLINGEN_SRC"
fi
[[ -s "$G4E" && -s "$CLINGEN" ]] || { echo "ERROR: missing input" >&2; exit 1; }

python3 - "$G4E" "$CLINGEN" "$OVR" "$TMP/panel.txt" "$G4E_SRC" "$CLINGEN_SRC" <<'PY'
import csv, re, sys, os, datetime, collections
g4e_f, cg_f, ovr_f, out_f, g4e_src, cg_src = sys.argv[1:7]

RANK = {'Definitive':6,'Strong':5,'Moderate':4,'Limited':3,'Disputed':2,'Refuted':1,
        'No Known Disease Relationship':0}
# Disease labels relevant to an epilepsy panel. "Leigh syndrome" is deliberately absent:
# ClinGen curates it as its own entity, and a Limited Leigh curation says nothing about the
# gene's epilepsy phenotype (POLG, CLPB, VPS13D would otherwise be removed).
REL = re.compile(r'epilep|seizure|encephalopath|neurodevelopment|intellectual|developmental delay|'
                 r'dravet|spasm|rett|angelman|ceroid|myoclon|lissenceph|polymicrogyr|heterotop|'
                 r'cortical|megalenceph|microceph|schizenceph|pachygyr|tuberous|interferonopathy|'
                 r'aicardi|leukoenceph|leukodystroph|hypomyelin|neurodegenerat|mitochondrial disease|'
                 r'glycosylation|glycosylphosph|hypermotor|febrile|lafora|pitt-hopkins', re.I)

def panel_moi(s):
    return s.replace('AD/AR', 'AD, AR').strip()

def compatible(pmoi, cmoi):
    p = {x.strip() for x in pmoi.split(',')}
    if cmoi in ('', 'Other', 'Undetermined'): return True
    if cmoi == 'SD': return bool(p & {'AD', 'AR'})
    if cmoi in ('XL', 'XLR', 'XLD'): return bool(p & {'XL', 'XLR', 'XLD'})
    return cmoi in p

# ── inputs ──
g4e = list(csv.DictReader(open(g4e_f, encoding='utf-8'), delimiter='\t'))
for col in ('Gene', 'Inheritance', 'Phenotype(s)'):
    if col not in g4e[0]: sys.exit(f"ERROR: Genes4Epilepsy file lacks column {col}")
m = re.search(r'v(\d{4}-\d{2})', os.path.basename(g4e_f))
g4e_ver = 'v' + m.group(1) if m else os.path.basename(g4e_f)

cg = collections.defaultdict(list); hdr = None; cg_date = 'unknown'
for r in csv.reader(open(cg_f, encoding='utf-8')):
    if not r: continue
    if r[0].startswith('FILE CREATED:'): cg_date = r[0].split(':', 1)[1].strip()
    elif r[0] == 'GENE SYMBOL': hdr = r
    elif hdr and not r[0].startswith('+'):
        d = dict(zip(hdr, r)); cg[d['GENE SYMBOL']].append(d)
if hdr is None: sys.exit("ERROR: ClinGen export has no header row")

ovr = collections.defaultdict(dict)
if os.path.exists(ovr_f):
    for l in open(ovr_f, encoding='utf-8'):
        if l.startswith('#') or not l.strip(): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 4 or f[1] not in ('MOI', 'KEEP', 'REMOVE'):
            sys.exit(f"ERROR: bad override line: {l!r}")
        ovr[f[0]][f[1]] = (f[2], f[3])

def pick(asserts, moi):
    def key(a):
        return (bool(REL.search(a['DISEASE LABEL'])), compatible(moi, a['MOI']),
                RANK.get(a['CLASSIFICATION'], -1), a['CLASSIFICATION DATE'])
    return max(asserts, key=key) if asserts else None

def fmt(a):
    return '|'.join([a['CLASSIFICATION'], a['DISEASE LABEL'], a['DISEASE ID (MONDO)'],
                     a['MOI'], a['CLASSIFICATION DATE'][:10]])

rows, removed, applied = [], [], []
for r in g4e:
    g = r['Gene'].strip(); moi = panel_moi(r['Inheritance'])
    if 'MOI' in ovr[g]:
        new = ovr[g]['MOI'][0]
        if new != moi: applied.append(f"{g} MOI {moi} -> {new}")
        moi = new
    a = cg.get(g, [])
    relc = [x for x in a if REL.search(x['DISEASE LABEL']) and compatible(moi, x['MOI'])]
    anyc = [x for x in a if compatible(moi, x['MOI'])]
    best_rel = max((RANK.get(x['CLASSIFICATION'], -1) for x in relc), default=None)
    best_any = max((RANK.get(x['CLASSIFICATION'], -1) for x in anyc), default=-1)
    why = None
    if 'REMOVE' in ovr[g]:
        why = 'override: ' + ovr[g]['REMOVE'][1]
    elif best_rel is not None and best_rel <= RANK['Limited'] and best_any <= RANK['Limited']:
        if 'KEEP' in ovr[g]: applied.append(f"{g} KEEP despite ClinGen {fmt(pick(relc, moi))}")
        else: why = fmt(pick(relc, moi))
    if why:
        removed.append((g, why)); continue
    best = pick(a, moi)
    rows.append((g, r['Phenotype(s)'].strip(), moi, fmt(best) if best else 'NOT_CURATED'))

rows.sort()
today = datetime.date.today().isoformat()
with open(out_f, 'w', encoding='utf-8') as o:
    o.write(f"# Default epilepsy gene panel for filtering_r.pl, generated by update_panel.sh on {today}.\n")
    o.write(f"## panel_version=Genes4Epilepsy {g4e_ver}\n")
    o.write(f"## clingen_fileDate={cg_date}\n")
    o.write(f"# Genes4Epilepsy source: {g4e_src} (Oliver KL et al., Epilepsia 2023).\n")
    o.write(f"# ClinGen gene-disease validity source: {cg_src}\n")
    o.write("# Columns: gene <tab> Association(Phenotype) <tab> MOI(Inheritance) <tab> GDV\n")
    o.write("#   GDV = CLASS|disease|MONDO|MOI|date of the ClinGen assertion chosen by update_panel.sh\n")
    o.write("#   (epilepsy/NDD-relevant > MOI-compatible > rank > recency); NOT_CURATED = no ClinGen curation.\n")
    o.write(f"# Genes: {len(rows)} kept of {len(g4e)} in Genes4Epilepsy {g4e_ver}; {len(removed)} removed.\n")
    for x in applied: o.write(f"# OVERRIDE {x}\n")
    for g, why in removed: o.write(f"# REMOVED {g}\t{why}\n")
    for row in rows: o.write('\t'.join(row) + '\n')
c = collections.Counter(r[3].split('|')[0] for r in rows)
print(f"[panel] Genes4Epilepsy {g4e_ver}: {len(g4e)} genes; ClinGen file {cg_date}")
print(f"[panel] kept {len(rows)}, removed {len(removed)}: {', '.join(g for g, _ in removed)}")
print("[panel] GDV classes: " + ', '.join(f"{k}={v}" for k, v in c.most_common()))
for x in applied: print(f"[panel] override: {x}")
PY
mv -f "$TMP/panel.txt" "$OUT"
echo "[panel] wrote $OUT"
