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
#   2026-10  - BAM recalculés à partir des reads sans ADN humain.
#            - Paramètres et échantillons centralisés (config.sh, samples.tsv).
#            - Les bins seront évalués (CheckM2) et affinés au lot 2.
#
# ============================================================

source /data/Martin/config.sh
conda_activate metabat2
mkdir -p "${OUT_ROOT}/cross-assembly"
cd "${OUT_ROOT}/cross-assembly"
# Attention : MetaBAT2 ne supprime pas les anciens bins. Si 03_bins/ existe déjà,
# le déplacer avant de relancer, sinon anciens et nouveaux bins se mélangent.
mkdir -p 03_bins


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

for CMD in metabat2; do
    if ! command -v "${CMD}" &> /dev/null; then
        echo "ERREUR : ${CMD} n'est pas installé ou n'est pas dans le PATH."
        exit 1
    fi
done


BAM_LIST=()

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_LIST+=("01_bam/${SAMPLE}.sorted.bam")
done


# ------------------------------------------------------------
# 9. PROFILS DE COUVERTURE POUR METABAT2
# ------------------------------------------------------------
#
# MetaBAT2 peut utiliser les BAM directement pour calculer
# la profondeur.
#
# On utilise jgi_summarize_bam_contig_depths pour produire
# la matrice de profondeur.
# ------------------------------------------------------------

echo ""
echo "=============================================="
echo "Étape 5 : matrice de profondeur MetaBAT2"
echo "=============================================="


if ! command -v jgi_summarize_bam_contig_depths &> /dev/null; then
    echo "ERREUR : jgi_summarize_bam_contig_depths n'est pas disponible."
    echo "Il est normalement installé avec MetaBAT2."
    exit 1
fi


jgi_summarize_bam_contig_depths \
    --outputDepth "02_depth/${PREFIX}_metabat_depth.txt" \
    "${BAM_LIST[@]}"


echo ""
echo "Matrice MetaBAT2 :"
echo "02_depth/${PREFIX}_metabat_depth.txt"


# ------------------------------------------------------------
# 10. BINNING AVEC METABAT2
# ------------------------------------------------------------

echo ""
echo "=============================================="
echo "Étape 6 : binning MetaBAT2"
echo "=============================================="


metabat2 \
    -i "${ASSEMBLY}" \
    -a "02_depth/${PREFIX}_metabat_depth.txt" \
    -o "03_bins/${PREFIX}_bin" \
    -m 1500 \
    -t "${THREADS}"


# ------------------------------------------------------------
# 11. RÉSUMÉ
# ------------------------------------------------------------

echo ""
echo "=============================================="
echo "PIPELINE TERMINÉ"
echo "=============================================="

echo ""
echo "Co-assemblage :"
echo "  ${ASSEMBLY}"

echo ""
echo "BAM :"
echo "  01_bam/"

echo ""
echo "Matrice CoverM :"
echo "  02_depth/${PREFIX}_coverm.tsv"

echo ""
echo "Matrice MetaBAT2 :"
echo "  02_depth/${PREFIX}_metabat_depth.txt"

echo ""
echo "Bins :"
echo "  03_bins/"

echo ""
echo "Nombre de bins :"
ls -1 03_bins/*.fa 2>/dev/null | wc -l

echo ""
echo "=============================================="
echo "Fin"
echo "=============================================="
