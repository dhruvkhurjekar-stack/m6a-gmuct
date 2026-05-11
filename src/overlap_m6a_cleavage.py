import pandas as pd
import pyranges as pr
import numpy as np
from pathlib import Path
from scipy.stats import binomtest, fisher_exact

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parents[1] / "data"
DATA_DIR.mkdir(exist_ok=True)

# ------------------------------------------------------------
# Arabidopsis TAIR10 nuclear chromosome sizes
# ------------------------------------------------------------
GENOME_SIZES = {
    "Chr1": 30427671,
    "Chr2": 19698289,
    "Chr3": 23459830,
    "Chr4": 18585056,
    "Chr5": 26975502,
}
TOTAL_GENOME_BP = sum(GENOME_SIZES.values())

# ------------------------------------------------------------
# Load high-confidence Col-0 m6A peaks
# File format:
# row 1 = title
# row 2 = blank
# row 3+ = data
# ------------------------------------------------------------
m6a_path = DATA_DIR / "col0m6asites.xlsx"

m6a_df = pd.read_excel(
    m6a_path,
    sheet_name=0,
    header=None,
    skiprows=2,
    usecols=[0, 1, 2, 3]
)

m6a_df.columns = ["Chromosome", "Start", "End", "Strand"]
m6a_df = m6a_df.dropna(subset=["Chromosome", "Start", "End"])
m6a_df["Chromosome"] = m6a_df["Chromosome"].astype(str).str.strip()
m6a_df["Start"] = m6a_df["Start"].astype(int)
m6a_df["End"] = m6a_df["End"].astype(int)
m6a_df["End"] = m6a_df["End"].where(m6a_df["End"] > m6a_df["Start"], m6a_df["Start"] + 1)

m6a_gr = pr.PyRanges(m6a_df[["Chromosome", "Start", "End"]])

print("Loaded Col-0 m6A peaks:", len(m6a_gr))
print(m6a_df.head())

# ------------------------------------------------------------
# Load NAR Supplementary Table 2 cleavage sites
# ------------------------------------------------------------
supp_path = DATA_DIR / "Supplementary Tables.xlsx"
xls = pd.ExcelFile(supp_path)
print("\nAvailable sheets:", xls.sheet_names)

header_row_idx = 3
df2_raw = pd.read_excel(supp_path, sheet_name="Sheet2", header=header_row_idx)

new_cols = df2_raw.iloc[0].astype(str).str.strip().tolist()
df2 = df2_raw.iloc[1:].copy()
df2.columns = new_cols
df2 = df2.dropna(how="all")

print("\nCleaned Sheet2 columns:", list(df2.columns))
print(df2.head(5).to_string(index=False))

# ------------------------------------------------------------
# Clean cleavage coordinates
# ------------------------------------------------------------
df2["Chromosome"] = df2["Chromosome"].astype(str).str.strip()
df2["Chromosome"] = df2["Chromosome"].apply(lambda x: f"Chr{x}" if x.isdigit() else x)

df2["Start"] = df2["Start"].astype(int)
df2["End"] = df2["End"].astype(int)
df2["End"] = df2["End"].where(df2["End"] > df2["Start"], df2["Start"] + 1)

if "log2 Fold-change" not in df2.columns:
    raise ValueError("Expected 'log2 Fold-change' column was not found.")

df2["FractionBias"] = np.where(df2["log2 Fold-change"] > 0, "Nucleoplasm", "Cytoplasm")

keep_cols = [
    "Chromosome", "Start", "End",
    "FractionBias", "log2 Fold-change",
    "Contrast", "Coverage", "Type", "ID", "Name"
]
keep_cols = [c for c in keep_cols if c in df2.columns]

cleavage_gr = pr.PyRanges(df2[keep_cols])

print("\nLoaded Table 2 cleavage sites:", len(cleavage_gr))
print("Cleavage chroms (sample):", cleavage_gr.df["Chromosome"].unique()[:5])
print("Cleavage Start range:", cleavage_gr.df["Start"].min(), cleavage_gr.df["Start"].max())
print("m6A Start range:", m6a_gr.df["Start"].min(), m6a_gr.df["Start"].max())
print("m6A End range:", m6a_gr.df["End"].min(), m6a_gr.df["End"].max())

# ------------------------------------------------------------
# Strict overlap
# ------------------------------------------------------------
overlap = cleavage_gr.join(m6a_gr)

print("\nOverlap columns:", overlap.df.columns.tolist())
print(overlap.df.head(10))

overlap_sites = overlap.df[["Chromosome", "Start", "End"]].drop_duplicates()

print("\n--- STRICT OVERLAP SUMMARY ---")
print("Total Table 2 cleavage sites:", len(df2))
print("Unique cleavage sites inside m6A peak coordinates:", len(overlap_sites))
print("Percent of Table 2 cleavage sites inside m6A peaks:",
      round(100 * len(overlap_sites) / len(df2), 4), "%")

counts = (
    overlap.df.groupby(["Chromosome", "Start", "End"])
    .size()
    .sort_values(ascending=False)
)
if len(counts) > 0:
    print("Max m6A peaks hit by a single cleavage site:", counts.iloc[0])
    print("Sites hitting >1 m6A peak:", (counts > 1).sum())
else:
    print("No overlaps found.")

# ------------------------------------------------------------
# FractionBias summary for overlapping sites
# ------------------------------------------------------------
unique_sites_by_bias = (
    overlap.df[["Chromosome", "Start", "End", "FractionBias"]]
    .drop_duplicates()
    .groupby("FractionBias")
    .size()
    .sort_values(ascending=False)
)

print("\nUnique overlapping cleavage sites by FractionBias:")
print(unique_sites_by_bias.to_string())

pct_bias = (
    overlap.df[["Chromosome", "Start", "End", "FractionBias"]]
    .drop_duplicates()["FractionBias"]
    .value_counts(normalize=True)
    .mul(100)
    .round(2)
)

print("\nPercent of unique overlapping cleavage sites by FractionBias:")
print(pct_bias.to_string())

# ------------------------------------------------------------
# Distance to nearest m6A peak
# ------------------------------------------------------------
nearest = cleavage_gr.nearest(m6a_gr)
distances = nearest.df["Distance"]

print("\n--- DISTANCE TO NEAREST m6A PEAK ---")
print(distances.describe())

print("\nNumber of cleavage sites within 50 bp of an m6A peak:",
      (distances.abs() <= 50).sum())
print("Number within 100 bp:",
      (distances.abs() <= 100).sum())

# ------------------------------------------------------------
# Genome coverage / expected-by-chance analysis
# ------------------------------------------------------------
m6a_merged = m6a_gr.merge()
m6a_bp_covered = (m6a_merged.df["End"] - m6a_merged.df["Start"]).sum()
p_null = m6a_bp_covered / TOTAL_GENOME_BP

observed_k = len(overlap_sites)
n_sites = len(df2)
expected_overlaps = n_sites * p_null

binom_result = binomtest(observed_k, n_sites, p=p_null, alternative="greater")

print("\n--- GENOME COVERAGE / BINOMIAL TEST ---")
print("Total genome bp:", TOTAL_GENOME_BP)
print("m6A bp covered (merged peaks):", m6a_bp_covered)
print("Fraction of genome covered by m6A peaks:", round(p_null, 6))
print("Observed overlaps:", observed_k)
print("Expected overlaps by chance:", round(expected_overlaps, 2))
print("Observed / expected:", round(observed_k / expected_overlaps, 3) if expected_overlaps > 0 else np.nan)
print("Binomial enrichment p-value:", binom_result.pvalue)

# ------------------------------------------------------------
# Fisher exact test for nucleoplasm vs cytoplasm enrichment
# ------------------------------------------------------------
total_bias_counts = df2["FractionBias"].value_counts()
overlap_bias_counts = (
    overlap.df[["Chromosome", "Start", "End", "FractionBias"]]
    .drop_duplicates()["FractionBias"]
    .value_counts()
)

nuc_overlap = int(overlap_bias_counts.get("Nucleoplasm", 0))
cyto_overlap = int(overlap_bias_counts.get("Cytoplasm", 0))

nuc_nonoverlap = int(total_bias_counts.get("Nucleoplasm", 0) - nuc_overlap)
cyto_nonoverlap = int(total_bias_counts.get("Cytoplasm", 0) - cyto_overlap)

contingency = [
    [nuc_overlap, nuc_nonoverlap],
    [cyto_overlap, cyto_nonoverlap],
]

oddsratio, fisher_p = fisher_exact(contingency, alternative="greater")

print("\n--- FISHER TEST: NUCLEOPLASM ENRICHMENT AMONG OVERLAPS ---")
print("Contingency table [[nuc_overlap, nuc_nonoverlap], [cyto_overlap, cyto_nonoverlap]]:")
print(contingency)
print("Odds ratio:", oddsratio)
print("Fisher exact p-value:", fisher_p)

# ------------------------------------------------------------
# mRNA-only analysis
# ------------------------------------------------------------
if "Type" in df2.columns:
    df2_mrna = df2[df2["Type"] == "mRNA"].copy()
    cleavage_gr_mrna = pr.PyRanges(df2_mrna[keep_cols])

    overlap_mrna = cleavage_gr_mrna.join(m6a_gr)
    overlap_sites_mrna = overlap_mrna.df[["Chromosome", "Start", "End"]].drop_duplicates()

    nearest_mrna = cleavage_gr_mrna.nearest(m6a_gr)
    dist_mrna = nearest_mrna.df["Distance"]

    observed_k_mrna = len(overlap_sites_mrna)
    n_sites_mrna = len(df2_mrna)
    expected_overlaps_mrna = n_sites_mrna * p_null
    binom_result_mrna = binomtest(observed_k_mrna, n_sites_mrna, p=p_null, alternative="greater")

    print("\n--- mRNA-ONLY ANALYSIS ---")
    print("Total mRNA cleavage sites:", n_sites_mrna)
    print("mRNA cleavage sites inside m6A peaks:", observed_k_mrna)
    print("Percent mRNA cleavage sites inside m6A peaks:",
          round(100 * observed_k_mrna / n_sites_mrna, 4) if n_sites_mrna > 0 else np.nan, "%")
    print("Expected mRNA overlaps by chance:", round(expected_overlaps_mrna, 2))
    print("Observed / expected (mRNA):",
          round(observed_k_mrna / expected_overlaps_mrna, 3) if expected_overlaps_mrna > 0 else np.nan)
    print("Binomial enrichment p-value (mRNA):", binom_result_mrna.pvalue)
    print("mRNA cleavage sites within 50 bp of an m6A peak:", (dist_mrna.abs() <= 50).sum())
    print("mRNA cleavage sites within 100 bp of an m6A peak:", (dist_mrna.abs() <= 100).sum())

    overlap_bias_mrna = (
        overlap_mrna.df[["Chromosome", "Start", "End", "FractionBias"]]
        .drop_duplicates()["FractionBias"]
        .value_counts()
    )
    total_bias_mrna = df2_mrna["FractionBias"].value_counts()

    nuc_overlap_mrna = int(overlap_bias_mrna.get("Nucleoplasm", 0))
    cyto_overlap_mrna = int(overlap_bias_mrna.get("Cytoplasm", 0))
    nuc_nonoverlap_mrna = int(total_bias_mrna.get("Nucleoplasm", 0) - nuc_overlap_mrna)
    cyto_nonoverlap_mrna = int(total_bias_mrna.get("Cytoplasm", 0) - cyto_overlap_mrna)

    contingency_mrna = [
        [nuc_overlap_mrna, nuc_nonoverlap_mrna],
        [cyto_overlap_mrna, cyto_nonoverlap_mrna],
    ]
    oddsratio_mrna, fisher_p_mrna = fisher_exact(contingency_mrna, alternative="greater")

    print("mRNA Fisher table:", contingency_mrna)
    print("mRNA Fisher odds ratio:", oddsratio_mrna)
    print("mRNA Fisher p-value:", fisher_p_mrna)

# ------------------------------------------------------------
# Permutation test by chromosome
# ------------------------------------------------------------
rng = np.random.default_rng(42)
n_perm = 1000

# Merge m6A peaks for interval lookup
m6a_merged_df = m6a_merged.df[["Chromosome", "Start", "End"]].copy()

# Count observed overlaps by simple interval inclusion after deduplication
obs_df = overlap_sites.copy()
observed_perm_stat = len(obs_df)

# Precompute chromosome-specific merged peak arrays
peak_lookup = {}
for chrom in GENOME_SIZES:
    tmp = m6a_merged_df[m6a_merged_df["Chromosome"] == chrom]
    peak_lookup[chrom] = list(zip(tmp["Start"].to_numpy(), tmp["End"].to_numpy()))

# Unique cleavage sites only for permutation
perm_df = df2[["Chromosome", "Start", "End"]].drop_duplicates().copy()
perm_df["Width"] = perm_df["End"] - perm_df["Start"]

def count_random_overlaps(random_df, peak_lookup_dict):
    count = 0
    for row in random_df.itertuples(index=False):
        chrom_peaks = peak_lookup_dict.get(row.Chromosome, [])
        for s, e in chrom_peaks:
            if row.Start < e and row.End > s:
                count += 1
                break
    return count

perm_counts = []

for _ in range(n_perm):
    rand_df = perm_df.copy()

    new_starts = []
    new_ends = []

    for row in rand_df.itertuples(index=False):
        chrom_len = GENOME_SIZES[row.Chromosome]
        width = int(row.Width)
        max_start = max(1, chrom_len - max(width, 1))
        rand_start = rng.integers(0, max_start)
        rand_end = rand_start + max(width, 1)
        new_starts.append(rand_start)
        new_ends.append(rand_end)

    rand_df["Start"] = new_starts
    rand_df["End"] = new_ends

    perm_counts.append(count_random_overlaps(rand_df, peak_lookup))

perm_counts = np.array(perm_counts)
perm_p = (np.sum(perm_counts >= observed_perm_stat) + 1) / (n_perm + 1)

print("\n--- PERMUTATION TEST ---")
print("Observed unique overlaps:", observed_perm_stat)
print("Permutation mean overlaps:", round(perm_counts.mean(), 3))
print("Permutation std overlaps:", round(perm_counts.std(), 3))
print("Permutation min/max:", int(perm_counts.min()), int(perm_counts.max()))
print("Empirical permutation p-value:", perm_p)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
out_overlap_sites = DATA_DIR / "table2_cleavage_sites_inside_col0_m6a_peaks.csv"
overlap_sites.to_csv(out_overlap_sites, index=False)

out_overlap_pairs = DATA_DIR / "table2_cleavage_x_col0_m6a_overlap_pairs.csv"
overlap.df.to_csv(out_overlap_pairs, index=False)

out_bias = DATA_DIR / "table2_col0_m6a_overlap_fractionbias_summary.csv"
unique_sites_by_bias.reset_index(name="n_unique_sites").to_csv(out_bias, index=False)

summary_rows = [
    {"metric": "n_col0_m6a_peaks", "value": len(m6a_gr)},
    {"metric": "n_table2_cleavage_sites", "value": len(df2)},
    {"metric": "n_unique_overlaps", "value": observed_k},
    {"metric": "pct_unique_overlaps", "value": round(100 * observed_k / len(df2), 4)},
    {"metric": "m6a_bp_covered", "value": int(m6a_bp_covered)},
    {"metric": "genome_bp", "value": int(TOTAL_GENOME_BP)},
    {"metric": "p_null", "value": p_null},
    {"metric": "expected_overlaps_by_chance", "value": expected_overlaps},
    {"metric": "observed_over_expected", "value": observed_k / expected_overlaps if expected_overlaps > 0 else np.nan},
    {"metric": "binomial_pvalue", "value": binom_result.pvalue},
    {"metric": "perm_mean_overlaps", "value": perm_counts.mean()},
    {"metric": "perm_pvalue", "value": perm_p},
    {"metric": "n_within_50bp", "value": int((distances.abs() <= 50).sum())},
    {"metric": "n_within_100bp", "value": int((distances.abs() <= 100).sum())},
    {"metric": "nuc_overlap", "value": nuc_overlap},
    {"metric": "cyto_overlap", "value": cyto_overlap},
    {"metric": "fisher_pvalue_nucleoplasm_enrichment", "value": fisher_p},
]

out_summary = DATA_DIR / "table2_col0_m6a_analysis_summary.csv"
pd.DataFrame(summary_rows).to_csv(out_summary, index=False)

print("\nWrote:", out_overlap_sites)
print("Wrote:", out_overlap_pairs)
print("Wrote:", out_bias)
print("Wrote:", out_summary)