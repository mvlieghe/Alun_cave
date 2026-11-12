#!/bin/bash

# Classification Kraken2 des READS (et non des contigs) de chaque échantillon.
# À lancer après host_removal.sh ; Bracken ensuite (bracken_reads.sh).
#
# HISTORIQUE
#   2026-10  - Entrée : clean_reads/ (sans ADN humain) au lieu de fastp/.
#            - Ajout de --confidence 0.1 : un read n'est classé que si au moins
#              10 % de ses k-mers soutiennent l'attribution. Avant (0 par défaut),
#              un seul k-mer suffisait, alors que 85 à 94 % des k-mers des reads
#              n'ont aucune correspondance dans la base (communauté de grotte
#              mal représentée).
#            - Échantillons lus dans samples.tsv ; chemins dans config.sh.
#            - Sortie par read compressée (gzip).

source /data/Martin/config.sh
read_samples

echo "start kraken2"

# Reads sans ADN humain (un sous-dossier par échantillon, cf. host_removal.sh)

READSPATH="$READS_CLEAN"

# Export kraken2 database path (to be changed)

DBPATH="$DB_KRAKEN"

# Export main kraken2 output path (to be changed)

OUTPUTPATH="${OUT_ROOT}/kraken2"

# Extensions des reads (sorties de host_removal.sh)

fextension="_R1.clean.fastq.gz"
rextension="_R2.clean.fastq.gz"

# Parameters (to be changed)

threads=$THREADS       # défini dans config.sh

# (change nothing after this line)

# Nouveau dossier de sortie : les anciens résultats (reads_min100-q30-t3/)
# sont conservés pour comparaison

outputdir="reads_clean_conf0.1"


# Create output dir

mkdir -p "$OUTPUTPATH/$outputdir"

# Loop to run kraken2

for sampleid in "${SAMPLES[@]}"; do
        mkdir -p "$OUTPUTPATH/$outputdir/$sampleid"
        cd "$OUTPUTPATH/$outputdir/$sampleid"
        # --confidence 0.1 : voir HISTORIQUE. Sortie par read compressée (~500 Mo -> ~100 Mo).
        kraken2 --db "$DBPATH" --paired --gzip-compressed --threads "$threads" \
                --confidence 0.1 \
                --report "$sampleid"_report.tsv --use-names \
                "$READSPATH/$sampleid/${sampleid}${fextension}" \
                "$READSPATH/$sampleid/${sampleid}${rextension}" \
                | gzip > "$sampleid"_k2.tsv.gz
        echo "finish kraken2 for sample $sampleid"
done



echo "finish kraken2"
