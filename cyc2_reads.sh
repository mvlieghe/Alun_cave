#!/bin/bash
# ============================================================
# cyc2_reads.sh — gène Cyc2 compté directement dans les READS
# ============================================================
# Pourquoi : les ferro-oxydants (~0,3 %) sont trop rares pour être
# assemblés ; iron_genes.sh (co-assemblage) ne trouve aucun de leurs Cyc2.
# On aligne donc les reads (DIAMOND blastx) sur la base Cyc2 annotée
# (cyc2_refdb.sh) : chaque read Cyc2 est rattaché au taxon du génome de
# référence le plus proche (attribution compétitive).
#
# Seuils : identité >= 60 % (protéine), >= 80 % du read aligné, E <= 1e-5.
#   Assez souple pour des Gallionellaceae de grotte non séquencées,
#   assez strict pour la famille (pas l'espèce).
#
# Abondance en COPIES PAR GÉNOME :
#   couverture Cyc2 = Σ (longueur du read / longueur du gène en nt)
#   équivalents-génomes = bases bactériennes+archéennes / taille moyenne
#                         de génome  (singlem/microbial_fraction.tsv)
#   copies par génome = couverture Cyc2 / équivalents-génomes
#   -> 0,003 = ~0,3 % des génomes portent un Cyc2 (de ce taxon).
#
# Entrées : clean_reads/, /db/cyc2_gtdb_r220/, singlem/microbial_fraction.tsv
# Sorties : cyc2_reads/<s>.hits.tsv          (sorties DIAMOND brutes)
#           cyc2_reads/cyc2_long.tsv         (échantillon x cluster x famille x genre)
#           cyc2_reads/cyc2_summary.tsv      (une ligne par échantillon ; utilisé par les statistiques)
#
# Lancement :
#   TEST=1 bash /data/Martin/cyc2_reads.sh
#   screen -S alun_cyc2 -L -Logfile /data/Martin/logs/cyc2_reads.screen.log \
#       bash /data/Martin/cyc2_reads.sh
#
# HISTORIQUE
#   2026-10  création
# ============================================================

source /data/Martin/config.sh
conda_activate bakta          # contient DIAMOND
read_samples

DB="/db/cyc2_gtdb_r220"
OUT="${OUT_ROOT}/cyc2_reads"
mkdir -p "${OUT}"

for s in "${SAMPLES[@]}"; do
    H="${OUT}/${s}.hits.tsv"
    if [[ -s "${H}" ]]; then
        echo "${s} : déjà fait -> ignoré"; continue
    fi
    : > "${H}.tmp"
    for mate in R1 R2; do
        diamond blastx \
            -q "${READS_CLEAN}/${s}/${s}_${mate}.clean.fastq.gz" \
            -d "${DB}/cyc2_refs.dmnd" \
            --id 60 --query-cover 80 --evalue 1e-5 \
            --max-target-seqs 1 --threads "${THREADS}" --quiet \
            --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
            >> "${H}.tmp"
    done
    mv "${H}.tmp" "${H}"
    echo "${s} : $(wc -l < "${H}") reads Cyc2"
done

# Agrégation et normalisation (python de base : pandas)
conda deactivate 2>/dev/null || true
python3 "${BASE}/cyc2_summary.py" "${OUT}" "${DB}/cyc2_taxonomy.tsv" \
    "${OUT_ROOT}/singlem/microbial_fraction.tsv" "${BASE}/iron_oxidizers.tsv" "${SAMPLES[@]}"
