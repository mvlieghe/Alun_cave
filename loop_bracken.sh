#!/bin/bash

# ============================================================
# OBSOLÈTE (2026-10) — ne plus utiliser. Bracken sur les contigs ; remplacé par bracken_reads.sh.
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

threads=20
readlen=100
level="S"
threshold=5

# (change nothing after this line)

# Export output name

cd $CONTIGSPATH
output=${PWD##*/}

# Loop to run kraken2

for folder in $CONTIGSPATH/*; do
        cd $folder
        sampleid=${PWD##*/}
	cd $OUTPUTPATH/$output/$sampleid
 #      kraken2 $CONTIGSPATH/$sampleid/"$sampleid"_$output$extension --db=$DBPATH --threads $threads --report "$sampleid"_"$output".kreport2 > "$sampleid"_"$output".kraken2
	level="G"
	bracken -d $DBPATH -i "$sampleid"_"$output".kreport2 -o "$sampleid"_"$output"-thresh5-$level.bracken -r $readlen -l $level -t $threshold
	level="F"
	bracken -d $DBPATH -i "$sampleid"_"$output".kreport2 -o "$sampleid"_"$output"-thresh5-$level.bracken -r $readlen -l $level -t 1
	level="O"
        bracken -d $DBPATH -i "$sampleid"_"$output".kreport2 -o "$sampleid"_"$output"-thresh5-$level.bracken -r $readlen -l $level -t 1
	level="C"
        bracken -d $DBPATH -i "$sampleid"_"$output".kreport2 -o "$sampleid"_"$output"-thresh5-$level.bracken -r $readlen -l $level -t 1
	level="P"
        bracken -d $DBPATH -i "$sampleid"_"$output".kreport2 -o "$sampleid"_"$output"-thresh5-$level.bracken -r $readlen -l $level -t 1
	sed -i "s/$/\t$sampleid/" *-thresh5-*
        echo "finish bracken for sample $sampleid"
done

cd $OUTPUTPATH/$output
cat */*-thresh5-G.bracken > bracken-G-full.tsv
cat */*-thresh5-F.bracken > bracken-F-full.tsv
cat */*-thresh5-O.bracken > bracken-O-full.tsv
cat */*-thresh5-C.bracken > bracken-C-full.tsv
cat */*-thresh5-P.bracken > bracken-P-full.tsv
echo "finish kraken2"
