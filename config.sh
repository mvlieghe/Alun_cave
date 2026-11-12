#!/bin/bash
# ============================================================
# config.sh — paramètres communs à tous les scripts du projet Alun
# ============================================================
# À charger en tête de chaque script :   source /data/Martin/config.sh
# Contient : chemins, bases de données, threads, liste des
# échantillons (samples.tsv) et fonctions utilitaires.
#
# HISTORIQUE
#   2026-10  création : centralise les chemins et la liste des
#            échantillons, auparavant recopiés dans chaque script.
#   2026-10  ajout des bases du lot 2 (SingleM, CheckM2, GTDB-Tk, FeGenie)
#   2026-10-16 LC_ALL=C remplacé par C.UTF-8 (accents illisibles dans les figures R)
# ============================================================

set -euo pipefail

# Locale : la locale allemande du serveur utilise la virgule décimale,
# ce qui fausse "sort -g" et certains calculs. On force le format C.
export LC_ALL=C.UTF-8     # point décimal (comme C) ET lecture correcte des accents (UTF-8)

# --- Chemins -------------------------------------------------
BASE="/data/Martin"                          # = /home/martin/martin (lien symbolique)
SAMPLES_TSV="${BASE}/samples.tsv"
READS_FASTP="${BASE}/fastp/min100-q30-t3"    # reads nettoyés par fastp
OUT_ROOT="${BASE}"                           # racine des sorties
READS_CLEAN="${OUT_ROOT}/clean_reads"        # reads fastp SANS ADN humain (host_removal.sh)
LOGS="${BASE}/logs"

# --- Bases de données ----------------------------------------
DB_KRAKEN="/db/krakenDB/kraken_std_old"
DB_SYLPH="/db/sylph/gtdb-r232-c200-dbv2.syl2db"
DB_HOSTILE="/db/hostile"
DB_SINGLEM="/db/singlem/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb"
DB_CHECKM2="/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
DB_GTDBTK="/db/gtdbtk_r232"
FEGENIE="/db/FeGenie"

# --- Ressources ----------------------------------------------
THREADS=20          # serveur partagé (80 cœurs) : rester raisonnable
MIN_FREE_GB=500     # espace minimum sur /data avant une étape lourde

# Fichiers temporaires sur /data plutôt que /tmp
export TMPDIR="${BASE}/tmp"
mkdir -p "${TMPDIR}" "${LOGS}"

# --- Mode test -----------------------------------------------
# TEST=1 bash script.sh  -> 1 seul échantillon, sorties dans test_run/
TEST="${TEST:-0}"
if [[ "${TEST}" == "1" ]]; then
    OUT_ROOT="${BASE}/test_run"
    READS_CLEAN="${OUT_ROOT}/clean_reads"
    mkdir -p "${OUT_ROOT}"
    echo ">>> MODE TEST : 1 échantillon, sorties dans ${OUT_ROOT}"
fi

# --- Fonctions -----------------------------------------------

# read_samples : remplit le tableau SAMPLES depuis samples.tsv
# (colonne 1 ; en-tête et lignes commençant par # ignorés)
read_samples() {
    mapfile -t SAMPLES < <(awk -F'\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $1}' "${SAMPLES_TSV}")
    if [[ "${TEST}" == "1" ]]; then
        SAMPLES=("${SAMPLES[0]}")
    fi
    echo "Échantillons (${#SAMPLES[@]}) : ${SAMPLES[*]}"
}

# conda_activate <env> : active /mnt/conda_ext4/envs/<env>.
# "set +u" pendant l'activation : les scripts de conda utilisent des
# variables non définies, ce qui ferait planter "set -u".
source /home/martin/miniconda3/etc/profile.d/conda.sh
conda_activate() {
    set +u
    conda activate "/mnt/conda_ext4/envs/$1"
    set -u
}

# check_disk : arrête le script s'il reste moins de MIN_FREE_GB Go sur /data
check_disk() {
    local free_gb
    free_gb=$(df --output=avail -BG "${BASE}" | tail -1 | tr -dc '0-9')
    if (( free_gb < MIN_FREE_GB )); then
        echo "ERREUR : ${free_gb} Go libres sur ${BASE} (minimum ${MIN_FREE_GB} Go)."
        exit 1
    fi
    echo "Espace libre : ${free_gb} Go"
}
