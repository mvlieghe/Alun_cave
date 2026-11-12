#!/bin/bash
# ============================================================
# check_suspect_taxa.sh — ces « espèces » Bracken sont-elles réelles ?
# ============================================================
# Pourquoi : Bracken attribue des abondances élevées à des pathogènes
# (Salmonella, E. coli, B. mallei, B. pertussis…) alors que Kraken2
# n'avait classé que quelques reads au rang espèce (voir README, §A).
# Test : on aligne les reads sur UN génome de référence par espèce
# (identité >= 95 %).
#  - Espèce réellement présente : reads répartis sur tout le génome,
#    donc fraction couverte proche de l'attendu 1 - exp(-profondeur).
#  - Artefact : reads concentrés sur quelques gènes conservés, donc
#    fraction couverte très inférieure à l'attendu.
#
# Sorties : suspect_taxa/refs/<espèce>.fna
#           suspect_taxa/coverm_genome.tsv
#           suspect_taxa/verdict.tsv
#
# Lancement :
#   screen -S alun_suspect -L -Logfile /data/Martin/logs/suspect.screen.log \
#       bash /data/Martin/check_suspect_taxa.sh
#
# HISTORIQUE
#   2026-10  création
# ============================================================

source /data/Martin/config.sh
read_samples

OUT="${OUT_ROOT}/suspect_taxa"
mkdir -p "${OUT}/refs"

# Un génome de référence RefSeq par espèce suspecte (souches types / de référence)
declare -A REFS=(
    [Salmonella_enterica]=GCF_000006945.2     # Typhimurium LT2
    [Escherichia_coli]=GCF_000005845.2        # K-12 MG1655
    [Klebsiella_pneumoniae]=GCF_000240185.1   # HS11286
    [Pseudomonas_aeruginosa]=GCF_000006765.1  # PAO1
    [Burkholderia_mallei]=GCF_000011705.1     # ATCC 23344
    [Burkholderia_ambifaria]=GCF_000203915.1  # AMMD
    [Bordetella_pertussis]=GCF_000195715.1    # Tohama I
)

# Téléchargement ciblé : 7 génomes (~35 Mo au total), pas de téléchargement en masse
conda_activate ncbi-datasets
for name in "${!REFS[@]}"; do
    acc="${REFS[$name]}"
    [[ -s "${OUT}/refs/${name}.fna" ]] && continue
    datasets download genome accession "${acc}" --include genome \
        --filename "${OUT}/refs/${name}.zip"
    unzip -p "${OUT}/refs/${name}.zip" "ncbi_dataset/data/${acc}/*.fna" > "${OUT}/refs/${name}.fna"
    rm "${OUT}/refs/${name}.zip"
done

# Alignement de tous les échantillons (aucun BAM conservé)
R1S=(); R2S=()
for s in "${SAMPLES[@]}"; do
    R1S+=("${READS_CLEAN}/${s}/${s}_R1.clean.fastq.gz")
    R2S+=("${READS_CLEAN}/${s}/${s}_R2.clean.fastq.gz")
done

conda_activate coverm
coverm genome \
    -1 "${R1S[@]}" -2 "${R2S[@]}" \
    --genome-fasta-files "${OUT}"/refs/*.fna \
    --methods relative_abundance mean covered_fraction \
    --min-read-percent-identity 95 --min-read-aligned-percent 75\
    --min-covered-fraction 0 \
    -t "${THREADS}" -o "${OUT}/coverm_genome.tsv"

# Verdict : fraction couverte observée / attendue (1 - exp(-profondeur))
python3 - "${OUT}/coverm_genome.tsv" "${OUT}/verdict.tsv" <<'EOF'
import sys, math, pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t")
rows = []
for col in [c for c in df.columns if c.endswith(" Mean")]:
    sample = col[: -len(" Mean")]
    for _, r in df.iterrows():
        if r["Genome"] == "unmapped":        # ligne ajoutée par coverm, sans objet ici
            continue
        mean, cf = r[col], r[f"{sample} Covered Fraction"]
        attendu = 1 - math.exp(-mean)
        if mean < 0.01:
            verdict = "trop peu de reads : absent ou indétectable"
        elif cf / attendu < 0.3:
            verdict = "ARTEFACT probable : reads concentrés sur peu de régions"
        elif cf / attendu > 0.7:
            verdict = "présence plausible : couverture homogène"
        else:
            verdict = "intermédiaire : à examiner"
        rows.append([sample, r["Genome"], mean, cf, round(attendu, 4),
                     round(cf / attendu, 3) if attendu > 0 else None, verdict])
pd.DataFrame(rows, columns=["sample", "genome", "mean_depth", "covered_fraction",
                            "covered_expected", "ratio", "verdict"]
             ).to_csv(sys.argv[2], sep="\t", index=False)
EOF
echo "Verdicts : ${OUT}/verdict.tsv"
