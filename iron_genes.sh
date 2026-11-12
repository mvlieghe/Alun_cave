#!/bin/bash
# ============================================================
# iron_genes.sh — gènes de la ferro-oxydation dans le co-assemblage
# ============================================================
# Pourquoi : appartenir à un genre « ferro-oxydant » ne garantit pas le
# trait. On cherche donc directement les gènes, avec les HMM de FeGenie
# (Garber et al. 2020, Front Microbiol) et leurs seuils de bitscore.
# Marqueur principal : cyc2, cytochrome de la membrane externe des
# ferro-oxydants neutrophiles. Attention : parmi ses 3 clusters, le
# cluster 3 est le moins clairement lié à l'oxydation du fer.
#
# Abondance exprimée en COPIES PAR GÉNOME :
#   Σ couverture des contigs portant le gène
#   / médiane des couvertures de 16 gènes de protéines ribosomiques
#     monocopie (Hug et al. 2016), présents une fois par génome.
# Cette normalisation corrige la profondeur de séquençage et se
# lit directement : 0,05 = 5 % des génomes portent le gène.
#
# Entrées : bakta/cross/cross.faa, cross.tsv  (annotation du co-assemblage)
#           cross-assembly/02_depth/min100-q30-t3_coverm.tsv  (couverture des contigs)
# Sorties : iron_genes/hmm/*.tbl                     (résultats hmmsearch bruts)
#           iron_genes/iron_gene_hits.tsv            (protéines retenues + contig + taxon Kraken)
#           iron_genes/iron_gene_copies_per_genome.tsv (gène x échantillon)
#
# Lancement :
#   screen -S alun_iron -L -Logfile /data/Martin/logs/iron.screen.log \
#       bash /data/Martin/iron_genes.sh
#
# HISTORIQUE
#   2026-10  création
#   2026-10  normalisation déplacée dans gene_copies.py (partagée avec mn_genes.sh) ;
#            seuil TIGRFAM « trusted » pour les HMM sans seuil FeGenie (MtrB_TIGR03509).
#   Interprétation : marqueur principal = Cyc2 clusters 1-2 ; MtoA/MtrB sont ambigus
#   (homologues de MtrA/MtrB des RÉDUCTEURS de fer), cluster 3 de Cyc2 peu fiable.
# ============================================================

source /data/Martin/config.sh
read_samples

OUT="${OUT_ROOT}/iron_genes"
mkdir -p "${OUT}/hmm"
FAA="${BASE}/bakta/cross/cross.faa"
BAKTA_TSV="${BASE}/bakta/cross/cross.tsv"
COVERM="${BASE}/cross-assembly/02_depth/min100-q30-t3_coverm.tsv"
K2NAMES="${BASE}/cross-assembly/final/k2_name.tsv"
HMMDIR="${FEGENIE}/hmms/iron/iron_oxidation"
CUTOFFS="${FEGENIE}/hmms/iron/HMM-bitcutoffs.txt"   # 2 colonnes : HMM, seuil de bitscore

# 1. hmmsearch : un HMM à la fois, avec le seuil de bitscore de FeGenie
for h in "${HMMDIR}"/*.hmm; do
    name=$(basename "${h}" .hmm)
    cut=$(awk -v n="${name}" '$1 == n { print $2 }' "${CUTOFFS}")
    if [[ -n "${cut}" ]]; then
        thr=(-T "${cut}")                  # seuil calibré par FeGenie
    elif grep -q '^TC ' "${h}"; then
        thr=(--cut_tc)                     # seuil « trusted » du modèle (ex. MtrB_TIGR03509)
    else
        echo "ATTENTION : pas de seuil pour ${name} -> E-value 1e-10 utilisée"
        thr=(-E 1e-10)
    fi
    hmmsearch --cpu "${THREADS}" "${thr[@]}" --noali \
        --tblout "${OUT}/hmm/${name}.tbl" "${h}" "${FAA}" > /dev/null
done

# 2. Tableau des correspondances (protéine, HMM, bitscore) à partir des sorties hmmsearch
awk 'BEGIN { OFS = "\t"; print "locus", "gene", "score" }
     !/^#/ { print $1, $3, $6 }' "${OUT}"/hmm/*.tbl > "${OUT}/hits_raw.tsv"

# 3. Normalisation en copies par génome (script partagé avec mn_genes.sh)
python3 "${BASE}/gene_copies.py" "${OUT}/hits_raw.tsv" "${OUT}" "${SAMPLES[@]}"

