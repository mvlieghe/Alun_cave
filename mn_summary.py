#!/usr/bin/env python3
# ============================================================
# mn_summary.py — reads d'oxydases du Mn (mn_reads.sh) -> copies par génome
# ============================================================
# copies par génome = Σ(longueur read / (3 x longueur protéine))
#                     / (bases bactéries+archées / taille moyenne de génome)
# Les reads dont le meilleur partenaire est un LEURRE sont comptés à part
# (colonne « leurres ») puis écartés.
# Usage : python3 mn_summary.py <dossier> <microbial_fraction.tsv> <échantillons...>
#
# HISTORIQUE
#   2026-10-15 création
# ============================================================
import sys
import pandas as pd

out, mf_file, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
COLS = ["read", "ref", "pident", "length", "qlen", "slen", "evalue", "bitscore"]
mf = pd.read_csv(mf_file, sep="\t").set_index("sample")
genome_eq = mf["bacterial_archaeal_bases"] / mf["average_bacterial_archaeal_genome_size"]

rows = []
for s in samples:
    h = pd.read_csv(f"{out}/{s}.hits.tsv", sep="\t", names=COLS)
    if h.empty:
        continue
    h["proteine"] = h["ref"].str.split("|").str[0]
    h["cov"] = h["qlen"] / (3 * h["slen"])
    g = h.groupby("proteine").agg(n_reads=("read", "size"), cov=("cov", "sum"),
                                  identite_moy=("pident", "mean")).reset_index()
    g["copies_per_genome"] = g["cov"] / genome_eq[s]
    g.insert(0, "sample", s)
    rows.append(g)

long = pd.concat(rows, ignore_index=True)
long.to_csv(f"{out}/mn_long.tsv", sep="\t", index=False)
wide = (long.pivot_table(index="proteine", columns="sample", values="copies_per_genome", fill_value=0)
            .reindex(columns=samples, fill_value=0))
wide.index.name = "gene"
wide.to_csv(f"{out}/mn_copies_par_genome.tsv", sep="\t")
pd.set_option("display.width", 250)
print(wide.round(4).to_string())
print("\nReads par protéine :")
print(long.groupby("proteine")["n_reads"].sum().to_string())
