import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ------------------------------------------------------------
# ANALYSIS 1: Coordinate-level Col-0 m6A overlap
# ------------------------------------------------------------

coordinate_summary = pd.DataFrame({
    "Analysis": ["All Table 2 cleavage sites", "mRNA cleavage sites only"],
    "Total_sites": [6504, 2243],
    "Observed_overlaps": [31, 31],
    "Expected_overlaps": [35.72, 12.32],
    "Observed_expected_ratio": [0.868, 2.517],
    "Percent_overlap": [0.4766, 1.3821],
    "Binomial_p_value": [0.8077, 5.18e-06],
})

coordinate_summary.to_csv(
    DATA_DIR / "publication_table_coordinate_col0_m6a_overlap.csv",
    index=False
)

print("\nCoordinate-level summary:")
print(coordinate_summary)

# Figure 1A: observed vs expected overlaps
fig, ax = plt.subplots(figsize=(6, 5))

x = range(len(coordinate_summary))
bar_width = 0.35

ax.bar(
    [i - bar_width / 2 for i in x],
    coordinate_summary["Observed_overlaps"],
    width=bar_width,
    label="Observed"
)

ax.bar(
    [i + bar_width / 2 for i in x],
    coordinate_summary["Expected_overlaps"],
    width=bar_width,
    label="Expected"
)

ax.set_xticks(list(x))
ax.set_xticklabels(coordinate_summary["Analysis"], rotation=20, ha="right")
ax.set_ylabel("Number of cleavage sites overlapping Col-0 m6A peaks")
ax.set_title("Observed vs expected coordinate overlap")
ax.legend()

for i, row in coordinate_summary.iterrows():
    ax.text(i - bar_width / 2, row["Observed_overlaps"] + 1, str(row["Observed_overlaps"]),
            ha="center", fontsize=9)
    ax.text(i + bar_width / 2, row["Expected_overlaps"] + 1, f"{row['Expected_overlaps']:.1f}",
            ha="center", fontsize=9)

plt.tight_layout()
plt.savefig(FIG_DIR / "figure1_coordinate_observed_vs_expected.png", dpi=300)
plt.close()

# Figure 1B: observed / expected ratio
fig, ax = plt.subplots(figsize=(6, 5))

ax.bar(
    coordinate_summary["Analysis"],
    coordinate_summary["Observed_expected_ratio"]
)

ax.axhline(1, linestyle="--", linewidth=1)

ax.set_ylabel("Observed / expected overlap")
ax.set_title("mRNA cleavage sites show stronger Col-0 m6A overlap")
ax.set_xticklabels(coordinate_summary["Analysis"], rotation=20, ha="right")

for i, v in enumerate(coordinate_summary["Observed_expected_ratio"]):
    ax.text(i, v + 0.05, f"{v:.2f}x", ha="center", fontsize=9)

plt.tight_layout()
plt.savefig(FIG_DIR / "figure2_coordinate_observed_expected_ratio.png", dpi=300)
plt.close()

# ------------------------------------------------------------
# ANALYSIS 1B: Nucleoplasm/cytoplasm split among coordinate overlaps
# ------------------------------------------------------------

coordinate_fraction_table = pd.DataFrame({
    "Fraction": ["Nucleoplasm", "Cytoplasm"],
    "mRNA_coordinate_overlaps": [30, 1],
})

coordinate_fraction_table["Percent"] = (
        coordinate_fraction_table["mRNA_coordinate_overlaps"]
        / coordinate_fraction_table["mRNA_coordinate_overlaps"].sum()
        * 100
).round(2)

coordinate_fraction_table.to_csv(
    DATA_DIR / "publication_table_coordinate_overlap_fractionbias.csv",
    index=False
)

print("\nCoordinate overlap fraction table:")
print(coordinate_fraction_table)

fig, ax = plt.subplots(figsize=(5, 5))

ax.bar(
    coordinate_fraction_table["Fraction"],
    coordinate_fraction_table["mRNA_coordinate_overlaps"]
)

ax.set_ylabel("Number of overlapping mRNA cleavage sites")
ax.set_title("Fraction bias among Col-0 m6A-overlapping mRNA cleavage sites")

for i, row in coordinate_fraction_table.iterrows():
    ax.text(
        i,
        row["mRNA_coordinate_overlaps"] + 0.5,
        f"{int(row['mRNA_coordinate_overlaps'])}\n({row['Percent']}%)",
        ha="center",
        fontsize=9
    )

plt.tight_layout()
plt.savefig(FIG_DIR / "figure3_coordinate_overlap_fractionbias.png", dpi=300)
plt.close()

# ------------------------------------------------------------
# ANALYSIS 2: Gene-level flower CS m6A overlap
# ------------------------------------------------------------

flower_table = pd.DataFrame({
    "Cleavage_gene_group": [
        "Nucleoplasmic mRNA cleavage genes",
        "Cytoplasmic mRNA cleavage genes"
    ],
    "Flower_m6A_plus": [477, 70],
    "Flower_m6A_minus": [299, 10],
})

flower_table["Total_genes"] = (
        flower_table["Flower_m6A_plus"] + flower_table["Flower_m6A_minus"]
)

flower_table["Percent_flower_m6A_plus"] = (
        flower_table["Flower_m6A_plus"] / flower_table["Total_genes"] * 100
).round(2)

flower_table["Fisher_exact_p_value"] = ["7.43e-07", "7.43e-07"]

flower_table.to_csv(
    DATA_DIR / "publication_table_flower_gene_level_m6a_overlap.csv",
    index=False
)

print("\nFlower gene-level summary:")
print(flower_table)

# Figure 4: percent flower m6A+
fig, ax = plt.subplots(figsize=(6, 5))

ax.bar(
    flower_table["Cleavage_gene_group"],
    flower_table["Percent_flower_m6A_plus"]
)

ax.set_ylabel("Percent of genes with flower CS m6A peak")
ax.set_title("Flower CS m6A peaks are more common in cytoplasmic mRNA cleavage genes")
ax.set_ylim(0, 100)
ax.set_xticklabels(flower_table["Cleavage_gene_group"], rotation=20, ha="right")

for i, row in flower_table.iterrows():
    ax.text(
        i,
        row["Percent_flower_m6A_plus"] + 2,
        f"{row['Percent_flower_m6A_plus']}%",
        ha="center",
        fontsize=9
    )

plt.tight_layout()
plt.savefig(FIG_DIR / "figure4_flower_percent_m6a_plus.png", dpi=300)
plt.close()

# Figure 5: stacked count bar
fig, ax = plt.subplots(figsize=(6, 5))

ax.bar(
    flower_table["Cleavage_gene_group"],
    flower_table["Flower_m6A_plus"],
    label="Flower m6A+"
)

ax.bar(
    flower_table["Cleavage_gene_group"],
    flower_table["Flower_m6A_minus"],
    bottom=flower_table["Flower_m6A_plus"],
    label="Flower m6A-"
)

ax.set_ylabel("Number of mRNA cleavage genes")
ax.set_title("Flower CS m6A status by cleavage compartment")
ax.set_xticklabels(flower_table["Cleavage_gene_group"], rotation=20, ha="right")
ax.legend()

plt.tight_layout()
plt.savefig(FIG_DIR / "figure5_flower_stacked_counts.png", dpi=300)
plt.close()

# ------------------------------------------------------------
# Combined table for quick discussion
# ------------------------------------------------------------

combined_summary = pd.DataFrame({
    "Analysis": [
        "Coordinate-level Col-0 m6A overlap: all Table 2 sites",
        "Coordinate-level Col-0 m6A overlap: mRNA sites only",
        "Gene-level flower CS m6A overlap: nucleoplasmic mRNA cleavage genes",
        "Gene-level flower CS m6A overlap: cytoplasmic mRNA cleavage genes",
    ],
    "Main_result": [
        "31/6504 sites overlap Col-0 m6A peaks",
        "31/2243 mRNA sites overlap Col-0 m6A peaks",
        "477/776 genes have flower CS m6A peaks",
        "70/80 genes have flower CS m6A peaks",
    ],
    "Percent": [
        "0.48%",
        "1.38%",
        "61.47%",
        "87.50%",
    ],
    "Statistic": [
        "Observed/expected = 0.87; binomial p = 0.81",
        "Observed/expected = 2.52; binomial p = 5.18e-06",
        "Fisher comparison vs cytoplasm p = 7.43e-07",
        "Fisher comparison vs nucleoplasm p = 7.43e-07",
    ],
    "Interpretation": [
        "No global enrichment across all cleavage sites",
        "Significant enrichment in mRNA cleavage sites",
        "Lower flower m6A frequency than cytoplasmic group",
        "Higher flower m6A frequency than nucleoplasmic group",
    ]
})

combined_summary.to_csv(
    DATA_DIR / "combined_publication_summary_table.csv",
    index=False
)

print("\nCombined summary:")
print(combined_summary)

print("\nDone. Figures saved in:")
print(FIG_DIR)
print("\nTables saved in:")
print(DATA_DIR)