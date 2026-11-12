#!/usr/bin/env python3
# ============================================================
# cyc2_summary.py — agrège les reads Cyc2 (cyc2_reads.sh) en copies par génome
# ============================================================
# Usage : python3 cyc2_summary.py <dossier> <cyc2_taxonomy.tsv>
#                 <microbial_fraction.tsv> <iron_oxidizers.tsv> <échantillons...>
# Sorties : <dossier>/cyc2_long.tsv, <dossier>/cyc2_summary.tsv
#
# copies par génome = Σ(longueur read / (3 x longueur protéine))
#                     / (bases bactéries+archées / taille moyenne de génome)
#
# HISTORIQUE
#   2026-10  création
# ============================================================
import re
import sys
import pandas as pd

out, taxo_file, mf_file, fe_file, samples = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
COLS = ["read", "ref", "pident", "length", "qlen", "slen", "evalue", "bitscore"]

# Taxonomie GTDB des génomes de référence
taxo = pd.read_csv(taxo_file, sep="\t", names=["genome", "lineage"]).set_index("genome")["lineage"]

def rank(lineage, prefix):
    """Extrait un rang (ex. 'f__') d'une lignée GTDB ; '' si absent."""
    for part in str(lineage).split(";"):
        if part.startswith(prefix):
            return part[len(prefix):]
    return ""

def clean(name):
    """Retire le suffixe GTDB (Pseudomonas_E -> Pseudomonas)."""
    return re.sub(r"_[A-Z]+$", "", name)

# Liste des ferro-oxydants de niveau 1 (familles ou genres)
fe = pd.read_csv(fe_file, sep="\t", comment="#")
fe1 = set(fe.loc[fe["niveau"] == 1, "nom"])

# Équivalents-génomes par échantillon (SingleM)
mf = pd.read_csv(mf_file, sep="\t").set_index("sample")
genome_eq = mf["bacterial_archaeal_bases"] / mf["average_bacterial_archaeal_genome_size"]

rows = []
for s in samples:
    h = pd.read_csv(f"{out}/{s}.hits.tsv", sep="\t", names=COLS)
    if h.empty:
        continue
    ref = h["ref"].str.split("|", expand=True)
    h["genome"], h["cluster"] = ref[0], ref[2]
    h["lineage"] = h["genome"].map(taxo)
    h["family"] = h["lineage"].map(lambda x: rank(x, "f__"))
    h["genus"] = h["lineage"].map(lambda x: rank(x, "g__"))
    h["cov"] = h["qlen"] / (3 * h["slen"])          # contribution du read à la couverture du gène
    g = (h.groupby(["cluster", "family", "genus"])
           .agg(n_reads=("read", "size"), cov=("cov", "sum"))
           .reset_index())
    g["copies_per_genome"] = g["cov"] / genome_eq[s]
    g.insert(0, "sample", s)
    rows.append(g)

long = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
    columns=["sample", "cluster", "family", "genus", "n_reads", "cov", "copies_per_genome"])
long.to_csv(f"{out}/cyc2_long.tsv", sep="\t", index=False)

# Résumé par échantillon (utilisé par stats_col_vs_unc.R)
cl12 = long["cluster"].isin(["Cyc2_repCluster1", "Cyc2_repCluster2"])
is_fe1 = long["family"].map(clean).isin(fe1) | long["genus"].map(clean).isin(fe1)
summary = pd.DataFrame(index=samples)
summary.index.name = "sample"
summary["cyc2_cl12_total"] = long[cl12].groupby("sample")["copies_per_genome"].sum()
summary["cyc2_cl12_fe_niveau1"] = long[cl12 & is_fe1].groupby("sample")["copies_per_genome"].sum()
summary["cyc2_cl12_gallionellaceae"] = long[cl12 & (long["family"] == "Gallionellaceae")].groupby("sample")["copies_per_genome"].sum()
summary["cyc2_cl3_total"] = long[~cl12].groupby("sample")["copies_per_genome"].sum()
summary["reads_cyc2_cl12"] = long[cl12].groupby("sample")["n_reads"].sum()
summary = summary.fillna(0)
summary.to_csv(f"{out}/cyc2_summary.tsv", sep="\t")
pd.set_option("display.width", 200)
print(summary.round(5).to_string())
