#!/bin/bash

# Estimation Bracken des abondances à partir des rapports Kraken2 sur les
# reads (kraken2_reads.sh). Rangs calculés : P, C, O, F, G.
#
# HISTORIQUE
#   2026-10  - Rang S (espèce) supprimé : Bracken redistribue vers des espèces
#              très représentées dans la base les reads que Kraken2 n'a classés
#              qu'au rang genre ou famille (ex. B. mallei dans 2UY1 : 10 reads
#              Kraken -> 28 181 ; Salmonella dans 1B2 : 123 -> 10 547).
#              Les abondances au rang espèce n'ont pas de sens ici.
#            - Fusion des tables corrigée (un seul en-tête, script relançable).
#            - Boucle sur samples.tsv (avant : dossiers megahit, y compris
#              all/ et cross/).
#            - Entrée : kraken2/reads_clean_conf0.1/ (reads sans ADN humain).

source /data/Martin/config.sh
read_samples

echo "start kraken2"

# Export kraken2 database path (to be changed)

DBPATH="$DB_KRAKEN"

# Export main kraken2 output path (to be changed)

OUTPUTPATH="${OUT_ROOT}/kraken2"

# Parameters (to be changed)

threads=20
# Longueur de read pour Bracken : 100 = seule distribution de k-mers
# disponible dans la base (database100mers.kmer_distrib). Valide pour nos
# reads (~145 pb après fastp) : Bracken demande une valeur <= longueur réelle.
readlen=100
threshold=5

# (change nothing after this line)

# Dossier produit par kraken2_reads.sh

output="reads_clean_conf0.1"

# Loop to run kraken2

for sampleid in "${SAMPLES[@]}"; do
        cd "$OUTPUTPATH/$output/$sampleid"
        # Rang S (espèce) retiré : voir HISTORIQUE en tête du fichier.
	level="G"
	bracken -d $DBPATH -i "$sampleid"_report.tsv -o "$sampleid"_bracken_$level.tsv -r $readlen -l $level -t $threshold
	level="F"
	bracken -d $DBPATH -i "$sampleid"_report.tsv -o "$sampleid"_bracken_$level.tsv -r $readlen -l $level -t 1
	level="O"
        bracken -d $DBPATH -i "$sampleid"_report.tsv -o "$sampleid"_bracken_$level.tsv -r $readlen -l $level -t 1
	level="C"
        bracken -d $DBPATH -i "$sampleid"_report.tsv -o "$sampleid"_bracken_$level.tsv -r $readlen -l $level -t 1
	level="P"
        bracken -d $DBPATH -i "$sampleid"_report.tsv -o "$sampleid"_bracken_$level.tsv -r $readlen -l $level -t 1
        echo "finish bracken for sample $sampleid"
done

cd "$OUTPUTPATH/$output"
# Fusion des tables par rang : UN SEUL en-tête + colonne "sample".
# (Avant : "sed -i" ajoutait l'échantillon à chaque ligne, en-tête compris,
#  et "cat" empilait 11 en-têtes ; relancer le script ajoutait une 2e colonne.)
for level in P C O F G; do
        awk -F'\t' -v OFS='\t' '
                FNR == 1 { split(FILENAME, p, "/"); s = p[1]   # dossier = échantillon
                           if (NR == 1) print $0, "sample"
                           next }
                { print $0, s }' */*_bracken_"$level".tsv > bracken-"$level"-full.tsv
done

echo "finish kraken2"
