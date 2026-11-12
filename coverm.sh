#!/bin/bash

# ============================================================
# Co-assemblage → mapping → couverture → binning
# ============================================================
#
# Entrée :
#   min100-q30-t3.contigs.fa
#   sample01_R1.fastq.gz / sample01_R2.fastq.gz
#   ...
#
# Sorties :
#   01_bam/       BAM triés et indexés
#   02_depth/     matrices de couverture
#   03_bins/      bins MetaBAT2
#
# Logiciels nécessaires :
#   bwa-mem2
#   samtools
#   coverm
#   metabat2
#
# HISTORIQUE
#   2026-10  - Ajout de la méthode "count" (nécessaire à group.py).
#            - Paramètres et échantillons centralisés (config.sh, samples.tsv).
#
# ============================================================

source /data/Martin/config.sh
conda_activate coverm
mkdir -p "${OUT_ROOT}/cross-assembly"
cd "${OUT_ROOT}/cross-assembly"


# ------------------------------------------------------------
# 1. PARAMÈTRES À MODIFIER
# ------------------------------------------------------------

# Fichier de co-assemblage
ASSEMBLY="/home/martin/martin/megahit/min100-q30-t3/cross/cross_min100-q30-t3.contigs.fa"


# Nombre de threads : THREADS, défini dans config.sh

# Nombre de threads pour samtools sort
SORT_THREADS=8

# Minimum de longueur des contigs
MIN_CONTIG_LEN=1000

# Préfixe des fichiers de sortie
PREFIX="min100-q30-t3"

# ------------------------------------------------------------
# Liste des échantillons : lue dans /data/Martin/samples.tsv
# (fonction read_samples de config.sh -> tableau SAMPLES)
# ------------------------------------------------------------

read_samples


echo "=============================================="
echo "Pipeline co-assemblage → binning"
echo "=============================================="
echo "Assembly : ${ASSEMBLY}"
echo "Threads  : ${THREADS}"
echo ""


# ------------------------------------------------------------
# 3. VÉRIFICATION DES PROGRAMMES
# ------------------------------------------------------------

for CMD in coverm; do
    if ! command -v "${CMD}" &> /dev/null; then
        echo "ERREUR : ${CMD} n'est pas installé ou n'est pas dans le PATH."
        exit 1
    fi
done


# ------------------------------------------------------------
# 8. MATRICE DE COUVERTURE AVEC COVERM
# ------------------------------------------------------------

echo ""
echo "=============================================="
echo "Étape 4 : matrice de couverture"
echo "=============================================="

BAM_LIST=()

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_LIST+=("01_bam/${SAMPLE}.sorted.bam")
done


# "count" ajouté : nombre de reads par contig, utilisé par final/group.py
# pour calculer des abondances normalisées par la profondeur de chaque échantillon.
coverm contig \
    --bam-files "${BAM_LIST[@]}" \
    --methods mean covered_fraction count \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --threads "${THREADS}" \
    --output-file "02_depth/${PREFIX}_coverm.tsv"


echo ""
echo "Matrice de couverture :"
echo "02_depth/${PREFIX}_coverm.tsv"


