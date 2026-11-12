#!/bin/bash
# ============================================================
# host_removal.sh — retrait des reads humains avec hostile
# ============================================================
# Pourquoi : Kraken2 attribue DIRECTEMENT des reads à Homo sapiens
# (ex. 6 661 reads dans 1B2, sans redistribution Bracken) : c'est de
# l'ADN humain apporté à la manipulation. On le retire avant toute
# analyse (taxonomie, mapping, statistiques).
#
# Outil : hostile (Constantinides et al. 2023, Bioinformatics) ;
# bowtie2 contre T2T-CHM13 + HLA (index "human-t2t-hla"), conçu
# pour retirer l'humain en perdant le moins possible de reads microbiens.
#
# Entrée  : fastp/min100-q30-t3/<s>/<s>_min100-q30-t3_R{1,2}.fastq.gz
# Sorties : clean_reads/<s>/<s>_R{1,2}.clean.fastq.gz
#           clean_reads/<s>/<s>_hostile.json      (rapport hostile)
#           clean_reads/host_removal_stats.tsv    (résumé)
#
# Lancement :
#   TEST=1 bash /data/Martin/host_removal.sh       (100 000 paires de 1B2)
#   screen -S alun_hostile -L -Logfile /data/Martin/logs/hostile.screen.log \
#       bash /data/Martin/host_removal.sh
#
# HISTORIQUE
#   2026-10  création
# ============================================================

source /data/Martin/config.sh
conda_activate hostile
export HOSTILE_CACHE_DIR="${DB_HOSTILE}"
read_samples
check_disk

for s in "${SAMPLES[@]}"; do
    R1="${READS_FASTP}/${s}/${s}_min100-q30-t3_R1.fastq.gz"
    R2="${READS_FASTP}/${s}/${s}_min100-q30-t3_R2.fastq.gz"
    OUT="${READS_CLEAN}/${s}"
    mkdir -p "${OUT}"

    if [[ -f "${OUT}/${s}_R1.clean.fastq.gz" ]]; then
        echo "${s} : déjà traité -> ignoré"
        continue
    fi

    # Mode test : 100 000 premières paires seulement
    if [[ "${TEST}" == "1" ]]; then
        seqkit head -n 100000 "${R1}" -o "${OUT}/test_R1.fastq.gz"
        seqkit head -n 100000 "${R2}" -o "${OUT}/test_R2.fastq.gz"
        R1="${OUT}/test_R1.fastq.gz"
        R2="${OUT}/test_R2.fastq.gz"
    fi

    echo "== ${s} : hostile"
    # Si la version installée refuse une option, comparer avec "hostile clean --help" (étape A2)
    hostile clean \
        --fastq1 "${R1}" --fastq2 "${R2}" \
        --aligner bowtie2 --index human-t2t-hla \
        --threads "${THREADS}" \
        --output "${OUT}" \
        > "${OUT}/${s}_hostile.json"

    # hostile nomme ses sorties <entrée>.clean_1/_2.fastq.gz : on renomme
    mv "${OUT}"/*.clean_1.fastq.gz "${OUT}/${s}_R1.clean.fastq.gz"
    mv "${OUT}"/*.clean_2.fastq.gz "${OUT}/${s}_R2.clean.fastq.gz"
    rm -f "${OUT}"/test_R?.fastq.gz
done

# Résumé (une ligne par échantillon)
python3 - "${READS_CLEAN}" "${SAMPLES[@]}" <<'EOF'
import json, sys
base, samples = sys.argv[1], sys.argv[2:]
with open(f"{base}/host_removal_stats.tsv", "w") as out:
    out.write("sample\treads_in\treads_out\treads_removed\tprop_removed\n")
    for s in samples:
        rep = json.load(open(f"{base}/{s}/{s}_hostile.json"))
        rep = rep[0] if isinstance(rep, list) else rep
        out.write(f"{s}\t{rep['reads_in']}\t{rep['reads_out']}\t"
                  f"{rep['reads_removed']}\t{rep['reads_removed_proportion']}\n")
EOF
echo "Résumé : ${READS_CLEAN}/host_removal_stats.tsv"
