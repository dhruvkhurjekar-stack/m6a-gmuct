import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parents[1] / "data"

flower_path = DATA_DIR / "13059_2022_2814_MOESM3_ESM.xlsx"
supp_path = DATA_DIR / "Supplementary Tables.xlsx"

# ------------------------------------------------------------
# Load flower m6A table
# ------------------------------------------------------------
# Real headers are on Excel row 3, so use header=2
flower = pd.read_excel(flower_path, header=2)
flower.columns = flower.columns.astype(str).str.strip()

print("Flower table columns:")
print(flower.columns.tolist())
print(flower.head())

# Column A = gene_id
flower = flower.rename(columns={"gene_id": "GeneID"})

# Column H = flower control m6A peaks
flower_col = flower.columns[7]
print("Using Column H as flower m6A column:", flower_col)

flower["GeneID"] = flower["GeneID"].astype(str).str.strip()

# Treat any nonzero value in Column H as flower m6A present
flower["Flower_m6A_present"] = flower[flower_col].notna() & (flower[flower_col] != 0)

methylated_genes = set(flower.loc[flower["Flower_m6A_present"], "GeneID"])

print("Total genes in flower table:", flower["GeneID"].nunique())
print("Genes with flower m6A peak:", len(methylated_genes))

# ------------------------------------------------------------
# Load NAR Table 2 cleavage sites
# ------------------------------------------------------------
header_row_idx = 3
df2_raw = pd.read_excel(supp_path, sheet_name="Sheet2", header=header_row_idx)

new_cols = df2_raw.iloc[0].astype(str).str.strip().tolist()
cleavage = df2_raw.iloc[1:].copy()
cleavage.columns = new_cols
cleavage = cleavage.dropna(how="all")

print("\nCleavage columns:")
print(cleavage.columns.tolist())
print(cleavage.head())

# Clean gene IDs and fraction bias
cleavage = cleavage.rename(columns={"ID": "GeneID"})
cleavage = cleavage.dropna(subset=["GeneID"])
cleavage["GeneID"] = cleavage["GeneID"].astype(str).str.strip()

cleavage["log2 Fold-change"] = pd.to_numeric(cleavage["log2 Fold-change"], errors="coerce")
cleavage = cleavage.dropna(subset=["log2 Fold-change"])

cleavage["FractionBias"] = np.where(
    cleavage["log2 Fold-change"] > 0,
    "Nucleoplasm",
    "Cytoplasm"
)

print("\nFractionBias counts in all cleavage rows:")
print(cleavage["FractionBias"].value_counts())

print("\nFractionBias counts in mRNA cleavage rows:")
print(cleavage[cleavage["Type"] == "mRNA"]["FractionBias"].value_counts())

print("\nUnique mRNA genes by FractionBias:")
print(
    cleavage[cleavage["Type"] == "mRNA"]
    [["GeneID", "FractionBias"]]
    .drop_duplicates()
    ["FractionBias"]
    .value_counts()
)

# Optional but likely most relevant: mRNA only
cleavage_mrna = cleavage[cleavage["Type"] == "mRNA"].copy()

print("\nTotal cleavage rows with GeneID:", len(cleavage))
print("Total mRNA cleavage rows:", len(cleavage_mrna))
print("Unique mRNA cleavage genes:", cleavage_mrna["GeneID"].nunique())

# ------------------------------------------------------------
# Gene-level overlap with flower m6A genes
# ------------------------------------------------------------
# Separate unique mRNA cleavage genes by compartment
nuc_genes = set(
    cleavage_mrna.loc[cleavage_mrna["FractionBias"] == "Nucleoplasm", "GeneID"]
)

cyto_genes = set(
    cleavage_mrna.loc[cleavage_mrna["FractionBias"] == "Cytoplasm", "GeneID"]
)

print("\nUnique nucleoplasmic mRNA cleavage genes:", len(nuc_genes))
print("Unique cytoplasmic mRNA cleavage genes:", len(cyto_genes))

# Count flower m6A overlap
nuc_m6a_plus = len(nuc_genes & methylated_genes)
nuc_m6a_minus = len(nuc_genes - methylated_genes)

cyto_m6a_plus = len(cyto_genes & methylated_genes)
cyto_m6a_minus = len(cyto_genes - methylated_genes)

# Build main table
table = pd.DataFrame(
    {
        "Flower_m6A_plus": [nuc_m6a_plus, cyto_m6a_plus],
        "Flower_m6A_minus": [nuc_m6a_minus, cyto_m6a_minus],
    },
    index=["Nucleoplasm", "Cytoplasm"]
)

print("\n--- MAIN CONTINGENCY TABLE: mRNA cleavage genes ---")
print(table)

pct = table.div(table.sum(axis=1), axis=0).mul(100).round(2)

print("\n--- PERCENT TABLE ---")
print(pct)

# ------------------------------------------------------------
# Fisher exact test
# ------------------------------------------------------------
from scipy.stats import fisher_exact

fisher_table = [
    [nuc_m6a_plus, nuc_m6a_minus],
    [cyto_m6a_plus, cyto_m6a_minus],
]

odds_ratio, p_value = fisher_exact(fisher_table, alternative="less")

print("\n--- FISHER EXACT TEST ---")
print("Table [[nuc m6A+, nuc m6A-], [cyto m6A+, cyto m6A-]]:")
print(fisher_table)
print("Odds ratio:", odds_ratio)
print("p-value:", p_value)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
gene_level = pd.DataFrame({
    "GeneID": list(nuc_genes) + list(cyto_genes),
    "FractionBias": ["Nucleoplasm"] * len(nuc_genes) + ["Cytoplasm"] * len(cyto_genes)
})
gene_level["Flower_m6A_present"] = gene_level["GeneID"].isin(methylated_genes)

gene_out = DATA_DIR / "flower_m6a_mrna_cleavage_gene_overlap.csv"
table_out = DATA_DIR / "flower_m6a_nucleus_vs_cytoplasm_table.csv"
pct_out = DATA_DIR / "flower_m6a_nucleus_vs_cytoplasm_percent.csv"

gene_level.to_csv(gene_out, index=False)
table.to_csv(table_out)
pct.to_csv(pct_out)

print("\nWrote:", gene_out)
print("Wrote:", table_out)
print("Wrote:", pct_out)