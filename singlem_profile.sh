#!/bin/bash
# ============================================================
# singlem_profile.sh — PROFIL TAXONOMIQUE PRINCIPAL (SingleM)
# ============================================================
# Pourquoi SingleM (Woodcroft et al. 2025, Nat Biotechnol) :
#  sylph et Kraken2 ne décrivent que < 8 % des reads de nos échantillons,
#  car la communauté de la grotte est très éloignée des génomes de
#  référence. SingleM profile à partir de 59 gènes marqueurs universels
#  placés dans la taxonomie GTDB : un organisme sans génome proche est
#  compté quand même, et classé au rang le plus précis possible
#  (ex. une Gallionellaceae inconnue -> « f__Gallionellaceae »).
#
# Entrée  : clean_reads/<s>/<s>_R{1,2}.clean.fastq.gz
# Sorties : singlem/profiles/<s>.profile.tsv     (profil par échantillon)
#           singlem/microbial_fraction.tsv       (part des reads microbiens)
#           singlem/table-*.tsv                  (tables échantillon x taxon par rang)
#
# Lancement :
#   TEST=1 bash /data/Martin/singlem_profile.sh
#   screen -S alun_singlem -L -Logfile /data/Martin/logs/singlem.screen.log \
#       bash /data/Martin/singlem_profile.sh
#
# HISTORIQUE
#   2026-10  création
#   2026-10  microbial_fraction calculée sur le profil brut (avant renommage) ; reprises
#            indépendantes pour le profil et pour la fraction.
# ============================================================

source /data/Martin/config.sh
conda_activate singlem
read_samples
check_disk

OUT="${OUT_ROOT}/singlem"
mkdir -p "${OUT}/profiles" "${OUT}/otu" "${OUT}/microbial_fraction"

for s in "${SAMPLES[@]}"; do
    R1="${READS_CLEAN}/${s}/${s}_R1.clean.fastq.gz"
    R2="${READS_CLEAN}/${s}/${s}_R2.clean.fastq.gz"
    P_RAW="${OUT}/profiles/${s}.profile.raw.tsv"   # profil brut, noms SingleM (« 1B2_R1.clean »)
    P="${OUT}/profiles/${s}.profile.tsv"           # même profil, échantillon renommé « 1B2 »
    F="${OUT}/microbial_fraction/${s}.tsv"

    # 1. Profil (sauté s'il existe déjà)
    if [[ ! -s "${P_RAW}" ]]; then
        echo "== ${s} : singlem pipe"
        singlem pipe -1 "${R1}" -2 "${R2}" \
            --metapackage "${DB_SINGLEM}" \
            --threads "${THREADS}" \
            --otu-table "${OUT}/otu/${s}.otu.tsv" \
            -p "${P_RAW}"
    fi

    # 2. Fraction microbienne : à calculer sur le profil BRUT, car SingleM vérifie
    #    que le nom d'échantillon du profil est celui déduit des fichiers de reads.
    if [[ ! -s "${F}" ]] || [[ $(wc -l < "${F}") -lt 2 ]]; then
        singlem microbial_fraction -1 "${R1}" -2 "${R2}" -p "${P_RAW}" \
            --metapackage "${DB_SINGLEM}" \
            | awk -F'\t' -v OFS='\t' -v s="${s}" 'NR > 1 { $1 = s } { print }' > "${F}"
    fi

    # 3. Copie renommée pour les tables (colonne 1 = nom du projet)
    awk -F'\t' -v OFS='\t' -v s="${s}" 'NR > 1 { $1 = s } { print }' "${P_RAW}" > "${P}"
done

# Résumé des fractions microbiennes (une ligne par échantillon)
awk 'FNR == 1 && NR != 1 { next } { print }' "${OUT}"/microbial_fraction/*.tsv > "${OUT}/microbial_fraction.tsv"

# Tables échantillon x taxon, un fichier par rang
# (si l'option n'existe pas dans la version installée : voir « singlem summarise --help »)
singlem summarise \
    --input-taxonomic-profiles "${OUT}"/profiles/*.profile.tsv \
    --output-species-by-site-relative-abundance-prefix "${OUT}/table"

echo "Tables : ${OUT}/table-*.tsv"
