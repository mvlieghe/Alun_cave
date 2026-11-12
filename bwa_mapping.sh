#!/bin/bash

# ============================================================
# Co-assemblage → mapping → couverture → binning
# ============================================================
#
# Entrée :
#   cross_min100-q30-t3.contigs.fa (co-assemblage MEGAHIT des 11 échantillons)
#   clean_reads/<s>/<s>_R1.clean.fastq.gz / <s>_R2.clean.fastq.gz
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
#   2026-10  - Reads d'entrée : clean_reads/ (sans ADN humain).
#            - Paramètres et échantillons centralisés (config.sh, samples.tsv),
#              auparavant recopiés dans les 3 scripts du dossier.
#            - Statistiques de mapping sur les reads primaires uniquement.
#
# ============================================================

# Paramètres communs (chemins, THREADS, samples.tsv, LC_ALL, set -euo pipefail)
source /data/Martin/config.sh

# Les chemins de sortie (01_bam/, 02_depth/, logs/) sont relatifs à ce dossier
mkdir -p "${OUT_ROOT}/cross-assembly"
cd "${OUT_ROOT}/cross-assembly"
mkdir -p 01_bam 02_depth logs


# ------------------------------------------------------------
# 1. PARAMÈTRES À MODIFIER
# ------------------------------------------------------------

# Fichier de co-assemblage
ASSEMBLY="/home/martin/martin/megahit/min100-q30-t3/cross/cross_min100-q30-t3.contigs.fa"

# Répertoire des FASTQ sans ADN humain (host_removal.sh)
READS_DIR="${READS_CLEAN}"

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

for CMD in bwa-mem2 samtools; do
    if ! command -v "${CMD}" &> /dev/null; then
        echo "ERREUR : ${CMD} n'est pas installé ou n'est pas dans le PATH."
        exit 1
    fi
done


# ------------------------------------------------------------
# 4. VÉRIFICATION DES FICHIERS FASTQ
# ------------------------------------------------------------

for SAMPLE in "${SAMPLES[@]}"; do

    R1="${READS_DIR}/${SAMPLE}/${SAMPLE}_R1.clean.fastq.gz"
    R2="${READS_DIR}/${SAMPLE}/${SAMPLE}_R2.clean.fastq.gz"

    if [[ ! -f "${R1}" ]]; then
        echo "ERREUR : fichier absent : ${R1}"
        exit 1
    fi

    if [[ ! -f "${R2}" ]]; then
        echo "ERREUR : fichier absent : ${R2}"
        exit 1
    fi

done


# ------------------------------------------------------------
# 5. INDEXATION DU CO-ASSEMBLAGE POUR BWA-MEM2
# ------------------------------------------------------------

echo ""
echo "=============================================="
echo "Étape 1 : indexation du co-assemblage"
echo "=============================================="

if [[ ! -f "${ASSEMBLY}.0123" ]]; then
    bwa-mem2 index "${ASSEMBLY}"
else
    echo "Index BWA-MEM2 déjà présent → étape ignorée."
fi


# ------------------------------------------------------------
# 6. MAPPING DES 11 ÉCHANTILLONS
# ------------------------------------------------------------

echo ""
echo "=============================================="
echo "Étape 2 : mapping des reads"
echo "=============================================="

for SAMPLE in "${SAMPLES[@]}"; do

    R1="${READS_DIR}/${SAMPLE}/${SAMPLE}_R1.clean.fastq.gz"
    R2="${READS_DIR}/${SAMPLE}/${SAMPLE}_R2.clean.fastq.gz"


    BAM="01_bam/${SAMPLE}.sorted.bam"

    echo ""
    echo "----------------------------------------------"
    echo "Échantillon : ${SAMPLE}"
    echo "----------------------------------------------"

    if [[ -f "${BAM}" ]]; then
        echo "BAM déjà présent → mapping ignoré."
        continue
    fi

    bwa-mem2 mem \
        -t "${THREADS}" \
        "${ASSEMBLY}" \
        "${R1}" \
        "${R2}" \
        2> "logs/${SAMPLE}.bwa.log" \
    | samtools sort \
        -@ "${SORT_THREADS}" \
        -o "${BAM}" \
        -

    samtools index -@ "${SORT_THREADS}" "${BAM}"

done


# ------------------------------------------------------------
# 7. STATISTIQUES DE MAPPING
# ------------------------------------------------------------

echo ""
echo "=============================================="
echo "Étape 3 : statistiques de mapping"
echo "=============================================="

STATS="02_depth/mapping_stats.tsv"

echo -e "Sample\tTotal_reads\tMapped_reads\tMapped_percent" > "${STATS}"

for SAMPLE in "${SAMPLES[@]}"; do

    BAM="01_bam/${SAMPLE}.sorted.bam"

    # Reads PRIMAIRES uniquement (-F 0x900 exclut secondaires et supplémentaires ;
    # 0x904 exclut en plus les non mappés). Avant, "view -c" comptait aussi les
    # alignements secondaires/supplémentaires, ce qui gonflait le total.
    TOTAL=$(samtools view -c -F 0x900 "${BAM}")
    MAPPED=$(samtools view -c -F 0x904 "${BAM}")

    if [[ "${TOTAL}" -gt 0 ]]; then
        PERCENT=$(awk -v m="${MAPPED}" -v t="${TOTAL}" \
            'BEGIN {printf "%.2f",100*m/t}')
    else
        PERCENT="0"
    fi

    echo -e "${SAMPLE}\t${TOTAL}\t${MAPPED}\t${PERCENT}" >> "${STATS}"

done

echo ""
echo "Statistiques disponibles dans : ${STATS}"
