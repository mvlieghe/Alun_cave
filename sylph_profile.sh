#!/bin/bash
# ============================================================
# sylph_profile.sh — profil taxonomique PRINCIPAL (sylph + GTDB r232)
# ============================================================
# Pourquoi sylph plutôt que Kraken2/Bracken comme profil principal :
#  - sylph n'annonce une espèce que si le génome de référence est
#    retrouvé à >= ~95 % d'ANI sur l'ensemble de ses k-mers. Ce n'est
#    pas le cas de Kraken2, qui classe un read sur quelques k-mers, ni de
#    Bracken, qui redistribue les reads vers les espèces de la base
#    (ex. B. mallei dans 2UY1 : 10 reads Kraken -> 28 181 après Bracken) ;
#  - taxonomie GTDB, plus cohérente pour les bactéries environnementales.
# Limite : seuls les organismes proches d'un génome GTDB sont profilés.
#   L'option -u estime la part NON profilée (colonne Sequence_abundance).
#
# Entrée  : clean_reads/<s>/<s>_R{1,2}.clean.fastq.gz
# Sorties : sylph/sylph_profile.tsv              (sortie brute, par génome)
#           sylph/taxprof/<s>.sylphmpa            (lignées GTDB complètes)
#           sylph/sylph_relative_abundance.tsv    (tous rangs, d__ -> s__)
#           sylph/sylph_sequence_abundance.tsv    (idem, en % des reads)
# Les tables par rang seront extraites dans le script R (lot 2).
#
# Lancement :
#   TEST=1 bash /data/Martin/sylph_profile.sh
#   screen -S alun_sylph -L -Logfile /data/Martin/logs/sylph.screen.log \
#       bash /data/Martin/sylph_profile.sh
#
# HISTORIQUE
#   2026-10  création
# ============================================================

source /data/Martin/config.sh
read_samples
check_disk

OUT="${OUT_ROOT}/sylph"
mkdir -p "${OUT}/sketches" "${OUT}/taxprof"

R1S=(); R2S=()
for s in "${SAMPLES[@]}"; do
    R1S+=("${READS_CLEAN}/${s}/${s}_R1.clean.fastq.gz")
    R2S+=("${READS_CLEAN}/${s}/${s}_R2.clean.fastq.gz")
done

# 1. Esquisses (sketches) des reads ; -S = nom de chaque échantillon
conda_activate sylph
sylph sketch -1 "${R1S[@]}" -2 "${R2S[@]}" -S "${SAMPLES[@]}" \
    -d "${OUT}/sketches" -t "${THREADS}"

# 2. Profil contre GTDB r232 ; -u : estime la fraction non profilée
sylph profile "${DB_SYLPH}" "${OUT}"/sketches/*.sylsp \
    -u -t "${THREADS}" -o "${OUT}/sylph_profile.tsv"

# 3. Taxonomie GTDB (métadonnées dans /db/sylph, cf. config sylph-tax)
conda_activate sylph-tax
sylph-tax taxprof "${OUT}/sylph_profile.tsv" -t GTDB_r232 -o "${OUT}/taxprof/"
sylph-tax merge "${OUT}"/taxprof/*.sylphmpa --column relative_abundance \
    -o "${OUT}/sylph_relative_abundance.tsv"
sylph-tax merge "${OUT}"/taxprof/*.sylphmpa --column sequence_abundance \
    -o "${OUT}/sylph_sequence_abundance.tsv"

echo "Profils : ${OUT}/sylph_relative_abundance.tsv"
