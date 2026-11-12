#!/bin/bash
# ============================================================
# mn_genes.sh — gènes de l'oxydation du manganèse (EXPLORATOIRE)
# ============================================================
# Pourquoi : les dépôts noirs des grottes sont souvent des oxydes de Mn
# (Northup et al. 2003). Les gènes du fer (FeGenie) ne couvrent pas le Mn.
# Les oxydases de Mn connues sont des oxydases multicuivre (MnxG, MofA,
# MoxA, McoA) ou des peroxydases (MopA). Pas de HMM spécifique publié :
# on compare nos protéines aux protéines EXPÉRIMENTALEMENT caractérisées
# (/db/mn_oxidases/*.faa) avec DIAMOND et des seuils stricts.
# Limite : les oxydases multicuivre forment une immense famille (laccases,
# CueO…) ; une ressemblance ne prouve pas l'activité sur le Mn. Résultat
# EXPLORATOIRE, à confirmer par le contexte génomique (MAG, gènes voisins).
#
# Seuils : identité >= 35 %, couverture >= 70 % (requête ET référence),
#          E-value <= 1e-10. À 35 % d'identité sur toute la longueur, la
#          fonction est souvent conservée, mais sans garantie pour les oxydases multicuivre.
#
# Entrées : /db/mn_oxidases/*.faa, bakta/cross/cross.faa
# Sorties : mn_genes/diamond_hits.tsv, gene_hits.tsv, gene_copies_per_genome.tsv
#
# Lancement :
#   screen -S alun_mn -L -Logfile /data/Martin/logs/mn.screen.log \
#       bash /data/Martin/mn_genes.sh
#
# HISTORIQUE
#   2026-10  création
# ============================================================

source /data/Martin/config.sh
conda_activate bakta          # contient DIAMOND
read_samples

OUT="${OUT_ROOT}/mn_genes"
REFDIR="/db/mn_oxidases"
FAA="${BASE}/bakta/cross/cross.faa"
mkdir -p "${OUT}"

# 1. Base de référence : en-têtes préfixés par le nom du gène (nom du fichier)
for f in "${REFDIR}"/*.faa; do
    gene=$(basename "${f}" .faa)
    sed "s/^>/>${gene}|/" "${f}"
done > "${OUT}/mn_refs.faa"
diamond makedb --in "${OUT}/mn_refs.faa" -d "${OUT}/mn_refs" --quiet

# 2. Recherche (nos protéines contre les références)
diamond blastp -q "${FAA}" -d "${OUT}/mn_refs" \
    --id 35 --query-cover 70 --subject-cover 70 --evalue 1e-10 \
    --max-target-seqs 1 --threads "${THREADS}" --sensitive \
    --outfmt 6 qseqid sseqid pident length qcovhsp scovhsp evalue bitscore \
    -o "${OUT}/diamond_hits.tsv"

# 3. Tableau des correspondances (locus, gène, score) puis normalisation
awk 'BEGIN { OFS = "\t"; print "locus", "gene", "score" }
     { split($2, a, "|"); print $1, a[1], $8 }' "${OUT}/diamond_hits.tsv" > "${OUT}/hits_raw.tsv"
python3 "${BASE}/gene_copies.py" "${OUT}/hits_raw.tsv" "${OUT}" "${SAMPLES[@]}"
