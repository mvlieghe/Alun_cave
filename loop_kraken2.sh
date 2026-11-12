#!/bin/bash

# ============================================================
# OBSOLÈTE (2026-10) — ne plus utiliser. Kraken2 sur les contigs ; remplacé par kraken2_reads.sh (reads).
# Voir /data/Martin/README_ANALYSES.md
# ============================================================

# Script for executing kraken2 classification based on megahit assembly. Intended to be used after the script loop_megahit.sh

echo "start kraken2"

# Export assembly directory, directory should only contain one folder per sample and will be used as output prefix (to be changed)

CONTIGSPATH="/home/martin/martin/megahit/min100-q30-t3"

# Export kraken2 database path (to be changed)

DBPATH="/db/krakenDB/kraken_std_old"

# Export main kraken2 output path (to be changed)

OUTPUTPATH="/home/martin/martin/kraken2"

# Export contigs extension (to be changed)

extension=".contigs.fa"

# Parameters (to be changed)

threads=10

# (change nothing after this line)

# Export output name

cd $CONTIGSPATH
output=${PWD##*/}

# Create output dir

mkdir $OUTPUTPATH/$output

# Loop to run kraken2

for folder in $CONTIGSPATH/*; do
        cd $folder
        sampleid=${PWD##*/}
	mkdir $OUTPUTPATH/$output/$sampleid
	cd $OUTPUTPATH/$output/$sampleid
        k2 classify $CONTIGSPATH/$sampleid/"$sampleid"_$output$extension --db $DBPATH --threads $threads --unclassified-out "$sampleid"_"$output"_unclass.fna --classified-out "$sampleid"_"$output"_class.fna --output "$sampleid"_"$output"_k2.txt --report "$sampleid"_"$output"_k2-report.txt --log "$sampleid"_"$output"_k2-log.txt
        echo "finish kraken2 for sample $sampleid"
done


echo "finish kraken2"
