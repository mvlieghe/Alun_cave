#!/usr/bin/env python3
# ============================================================
# gene_copies.py — abondance de gènes en COPIES PAR GÉNOME (co-assemblage)
# ============================================================
# Utilisé par iron_genes.sh et mn_genes.sh.
#
# Principe :
#   copies(gène, échantillon) =
#       Σ couverture moyenne des contigs portant le gène
#     / médiane, sur 16 gènes de protéines ribosomiques monocopie
#       (Hug et al. 2016), de la somme des couvertures de leurs contigs
#   -> 0,05 = environ 5 % des génomes portent le gène ; corrige la
#      profondeur de séquençage (3UBr2 a ~2x plus de reads).
#
# Usage : python3 gene_copies.py <hits.tsv> <dossier_sortie> <échantillon1> ...
#   hits.tsv : 3 colonnes avec en-tête  locus  gene  score
#              (une ligne par correspondance ; la meilleure par protéine est gardée)
# Sorties  : <dossier>/gene_hits.tsv              (protéines retenues + contig + taxon Kraken)
#            <dossier>/gene_copies_per_genome.tsv (gène x échantillon)
#
# HISTORIQUE
#   2026-10  création (extrait d'iron_genes.sh pour être partagé)
# ============================================================
import sys
import pandas as pd

BASE = "/data/Martin"
BAKTA_TSV = f"{BASE}/bakta/cross/cross.tsv"
COVERM = f"{BASE}/cross-assembly/02_depth/min100-q30-t3_coverm.tsv"
K2NAMES = f"{BASE}/cross-assembly/final/k2_name.tsv"
RP16 = {"rplB", "rplC", "rplD", "rplE", "rplF", "rplN", "rplO", "rplP",
        "rplR", "rplV", "rplX", "rpsC", "rpsH", "rpsJ", "rpsQ", "rpsS"}

hits_file, out, samples = sys.argv[1], sys.argv[2], sys.argv[3:]

# Meilleure correspondance par protéine
hits = (pd.read_csv(hits_file, sep="\t")
          .sort_values("score", ascending=False)
          .drop_duplicates("locus"))

# Annotation Bakta : protéine -> contig ; repérage des protéines ribosomiques
cols = ["contig", "type", "start", "stop", "strand", "locus", "gene_bakta", "product", "dbxrefs"]
bakta = pd.read_csv(BAKTA_TSV, sep="\t", comment="#", names=cols, usecols=range(9))

# Couverture moyenne de chaque contig, par échantillon
cov = pd.read_csv(COVERM, sep="\t").set_index("Contig")
cov = cov[[f"{s}.sorted Mean" for s in samples]]
cov.columns = samples

# Profondeur d'un « génome moyen » : médiane sur les 16 gènes RP
rp = bakta[bakta["gene_bakta"].isin(RP16)][["gene_bakta", "contig"]]
rp_depth = rp.join(cov, on="contig").groupby("gene_bakta")[samples].sum().median()
print("Profondeur génome moyen (médiane RP16) :")
print(rp_depth.round(2).to_string())

# Gènes trouvés : contig, produit Bakta, taxon Kraken du contig, couvertures
h = hits.merge(bakta[["locus", "contig", "product"]], on="locus", how="left")
k2 = pd.read_csv(K2NAMES, sep="\t")
h = h.merge(k2, left_on="contig", right_on="Contig", how="left").drop(columns="Contig")
h = h.join(cov, on="contig")
h.to_csv(f"{out}/gene_hits.tsv", sep="\t", index=False)

copies = h.groupby("gene")[samples].sum() / rp_depth
copies.to_csv(f"{out}/gene_copies_per_genome.tsv", sep="\t")
print("\nCopies par génome :")
print(copies.round(4).to_string())
