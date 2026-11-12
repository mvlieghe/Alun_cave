#!/bin/bash
# ============================================================
# mn_reads.sh — oxydases du manganèse cherchées dans les READS (exploratoire)
# ============================================================
# Pourquoi : la recherche dans le co-assemblage (mn_genes.sh) ne voit pas
# les organismes rares, qui ne sont pas assemblés (cf. Gallionellaceae).
# Même principe que cyc2_reads.sh : DIAMOND blastx des reads contre des
# protéines de référence, normalisation en COPIES PAR GÉNOME (SingleM).
#
# Différence avec Cyc2 : il n'existe pas de HMM spécifique publié pour les
# oxydases du Mn, donc pas de grande base « compétitive » comme pour Cyc2.
# Pour limiter les faux positifs, on ajoute des LEURRES : des oxydases
# multicuivre qui N'oxydent PAS le Mn (CueO, CotA, CopA). Un read dont le
# meilleur partenaire est un leurre est écarté. Le seuil d'identité est
# aussi plus strict que pour Cyc2 (70 %).
# GroEL (CONTROLE_groEL) sert de TÉMOIN de normalisation : on attend
# environ 1 à 2 copies par génome.
#
# Entrées : /db/mn_oxidases/*.faa                (références, un fichier par protéine)
#           /db/mn_oxidases/leurres/*.faa        (leurres, voir MODIFICATIONS_lot4.md)
#           clean_reads/, singlem/microbial_fraction.tsv
# Sorties : mn_reads/<s>.hits.tsv, mn_reads/mn_long.tsv,
#           mn_reads/mn_copies_par_genome.tsv    (protéine x échantillon)
#
# Lancement :
#   TEST=1 bash /data/Martin/mn_reads.sh
#   screen -S alun_mnreads -L -Logfile /data/Martin/logs/mn_reads.screen.log \
#       bash /data/Martin/mn_reads.sh
#
# HISTORIQUE
#   2026-10-15 création
# ============================================================

source /data/Martin/config.sh
conda_activate bakta          # contient DIAMOND
read_samples

REF="/db/mn_oxidases"
OUT="${OUT_ROOT}/mn_reads"
mkdir -p "${OUT}"

# 1. Base : références + leurres ; en-tête = >NOM_DU_FICHIER|accession
{
    for f in "${REF}"/*.faa; do
        sed "s/^>/>$(basename "${f}" .faa)|/" "${f}"
    done
    for f in "${REF}"/leurres/*.faa; do
        sed "s/^>/>LEURRE_$(basename "${f}" .faa)|/" "${f}"
    done
} > "${OUT}/mn_refs.faa"
diamond makedb --in "${OUT}/mn_refs.faa" -d "${OUT}/mn_refs" --quiet
echo "Références : $(grep -c '>' "${OUT}/mn_refs.faa") protéines (dont $(grep -c '>LEURRE_' "${OUT}/mn_refs.faa") leurres)"

# 2. Recherche dans les reads (R1 et R2)
for s in "${SAMPLES[@]}"; do
    H="${OUT}/${s}.hits.tsv"
    if [[ -s "${H}" ]]; then echo "${s} : déjà fait -> ignoré"; continue; fi
    : > "${H}.tmp"
    for mate in R1 R2; do
        diamond blastx \
            -q "${READS_CLEAN}/${s}/${s}_${mate}.clean.fastq.gz" \
            -d "${OUT}/mn_refs.dmnd" \
            --id 70 --query-cover 80 --evalue 1e-5 \
            --max-target-seqs 1 --threads "${THREADS}" --quiet \
            --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
            >> "${H}.tmp"
    done
    mv "${H}.tmp" "${H}"
    echo "${s} : $(wc -l < "${H}") reads"
done

# 3. Agrégation et normalisation (python de base)
conda deactivate 2>/dev/null || true
python3 "${BASE}/mn_summary.py" "${OUT}" "${OUT_ROOT}/singlem/microbial_fraction.tsv" "${SAMPLES[@]}"
