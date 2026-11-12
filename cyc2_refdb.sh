#!/bin/bash
# ============================================================
# cyc2_refdb.sh — base de référence Cyc2 ANNOTÉE TAXONOMIQUEMENT (une seule fois)
# ============================================================
# Pourquoi : dans le co-assemblage, aucun Cyc2 de ferro-oxydant n'a été
# trouvé. Les Gallionellaceae sont trop rares (~0,3 % selon SingleM) pour
# être assemblées : leurs gènes restent dans les reads. On cherchera donc
# Cyc2 directement dans les READS (cyc2_reads.sh), contre une base de
# protéines Cyc2 dont on connaît le génome et donc le taxon.
#
# Construction : les 3 HMM Cyc2 de FeGenie sont recherchés dans TOUS les
# protéomes représentatifs de GTDB r220 (déjà sur le disque, ~113 000
# génomes), avec les seuils de bitscore de FeGenie. Prendre toute la base,
# et pas seulement les ferro-oxydants, permet une attribution
# COMPÉTITIVE : un read issu d'un Cyc2 de non-ferro-oxydant sera rattaché
# à son vrai taxon au lieu d'être attribué à tort aux Gallionellaceae.
# (r220 ici vs r232 pour SingleM : les noms de famille sont stables entre
#  les deux versions pour les taxons qui nous intéressent.)
#
# Sorties : /db/cyc2_gtdb_r220/cyc2_hits.tsv    (génome, protéine, HMM, bitscore)
#           /db/cyc2_gtdb_r220/cyc2_refs.faa    (en-tête : >génome|protéine|cluster)
#           /db/cyc2_gtdb_r220/cyc2_refs.dmnd   (base DIAMOND)
#           /db/cyc2_gtdb_r220/cyc2_taxonomy.tsv (génome -> taxonomie GTDB)
# Durée : plusieurs heures (hmmsearch sur ~113 000 protéomes, 20 en parallèle).
#
# Lancement :
#   screen -S alun_cyc2db -L -Logfile /data/Martin/logs/cyc2_refdb.screen.log \
#       bash /data/Martin/cyc2_refdb.sh
#
# HISTORIQUE
#   2026-10  création
#   2026-10-13 hmmsearch lit un fichier temporaire et non plus un tube (seul le cluster 1 était recherché) ;
#              contrôle des 3 clusters ajouté
# ============================================================

source /data/Martin/config.sh

GTDB_PROT="/db/GTDB/R220/protein_faa_reps"
OUT="/db/cyc2_gtdb_r220"
HMMDIR="${FEGENIE}/hmms/iron/iron_oxidation"
CUTOFFS="${FEGENIE}/hmms/iron/HMM-bitcutoffs.txt"
mkdir -p "${OUT}"

# 1. Un seul fichier HMM contenant les 3 clusters Cyc2
cat "${HMMDIR}"/Cyc2_repCluster{1,2,3}.hmm > "${OUT}/cyc2.hmm"

# 2. hmmsearch sur chaque protéome (20 en parallèle, 1 cœur chacun).
#    Seuil large (E <= 1e-5) ici ; le seuil FeGenie est appliqué à l'étape 3.
#    Chaque ligne de sortie est préfixée par l'identifiant du génome.
search_one() {
    f="$1"
    g=$(basename "${f}" _protein.faa.gz)
    # Protéome décompressé dans un fichier temporaire : hmmsearch relit la base
    # de séquences pour CHAQUE modèle du fichier .hmm (3 ici). Lue depuis un
    # tube (« - »), elle n'est disponible qu'une fois -> seul le 1er modèle
    # (cluster 1) était recherché.
    tmp=$(mktemp -p "${TMPDIR}" cyc2.XXXXXX.faa)
    zcat "${f}" > "${tmp}"
    hmmsearch --cpu 1 -E 1e-5 --noali -o /dev/null --tblout /dev/stdout \
        "${OUT}/cyc2.hmm" "${tmp}" 2>/dev/null \
        | awk -v g="${g}" '!/^#/ { print g "\t" $1 "\t" $3 "\t" $6 }'
    rm -f "${tmp}"
}
export -f search_one
export OUT

if [[ ! -s "${OUT}/cyc2_hits_raw.tsv" ]]; then
    find "${GTDB_PROT}/bacteria" "${GTDB_PROT}/archaea" -name '*_protein.faa.gz' \
        | xargs -P "${THREADS}" -I{} bash -c 'search_one "$@"' _ {} \
        > "${OUT}/cyc2_hits_raw.tsv.tmp"
    mv "${OUT}/cyc2_hits_raw.tsv.tmp" "${OUT}/cyc2_hits_raw.tsv"
fi


# Contrôle : les 3 clusters Cyc2 doivent avoir des correspondances
n_cl=$(cut -f3 "${OUT}/cyc2_hits_raw.tsv" | sort -u | wc -l)
if (( n_cl < 3 )); then
    echo "ERREUR : ${n_cl} cluster(s) Cyc2 seulement dans cyc2_hits_raw.tsv (3 attendus)."
    echo "Supprimer ${OUT}/cyc2_hits_raw.tsv et vérifier la recherche hmmsearch."
    exit 1
fi

# 3. Seuils FeGenie, puis meilleur cluster par protéine
awk -F'\t' -v OFS='\t' '
    NR == FNR { cut[$1] = $2; next }                  # fichier des seuils
    ($3 in cut) && $4 >= cut[$3] { print }
' "${CUTOFFS}" "${OUT}/cyc2_hits_raw.tsv" \
    | sort -t$'\t' -k1,1 -k2,2 -k4,4gr \
    | awk -F'\t' '!seen[$1 FS $2]++' > "${OUT}/cyc2_hits.tsv"
echo "Protéines Cyc2 retenues : $(wc -l < "${OUT}/cyc2_hits.tsv")"

# 4. Extraction des séquences (en-tête : >génome|protéine|cluster)
: > "${OUT}/cyc2_refs.faa"
cut -f1 "${OUT}/cyc2_hits.tsv" | sort -u | while read -r g; do
    dom=$([[ -f "${GTDB_PROT}/bacteria/${g}_protein.faa.gz" ]] && echo bacteria || echo archaea)
    awk -F'\t' -v g="${g}" '$1 == g { print $2 "\t" $3 }' "${OUT}/cyc2_hits.tsv" > "${TMPDIR}/ids.tsv"
    seqkit grep -f <(cut -f1 "${TMPDIR}/ids.tsv") "${GTDB_PROT}/${dom}/${g}_protein.faa.gz" \
        | seqkit replace -p '^(\S+).*' -r "${g}|\${1}|{kv}" -k "${TMPDIR}/ids.tsv" \
        >> "${OUT}/cyc2_refs.faa"
done

# 5. Taxonomie GTDB des génomes retenus
zcat /db/GTDB/R220/bac120_taxonomy_r220.tsv.gz /db/GTDB/R220/ar53_taxonomy_r220.tsv.gz \
    | awk -F'\t' 'NR == FNR { keep[$1]; next } ($1 in keep)' \
        <(cut -f1 "${OUT}/cyc2_hits.tsv" | sort -u) - > "${OUT}/cyc2_taxonomy.tsv"

# 6. Base DIAMOND
conda_activate bakta          # contient DIAMOND
diamond makedb --in "${OUT}/cyc2_refs.faa" -d "${OUT}/cyc2_refs" --quiet

echo "Base Cyc2 : ${OUT}/cyc2_refs.dmnd"
echo "Répartition par famille (top 15) :"
awk -F'\t' '{ split($2, t, ";"); print t[5] }' "${OUT}/cyc2_taxonomy.tsv" | sort | uniq -c | sort -nr | head -15
