# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import os
import sys
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import gzip
import re

# --- Snakemake I/O ---
raw_data_path = snakemake.input.raw_data
vcf_path = snakemake.input.vcf_path_for_lengths
output_csv = snakemake.output.data
report_csv = snakemake.output.report
plot_file = snakemake.output.plot
exclude_file = snakemake.output.exclude_scaffolds

# --- Directory and Logging Setup ---
output_dir = os.path.dirname(report_csv)
os.makedirs(output_dir, exist_ok=True)
print(f"Reading raw data from: {raw_data_path}")
print(f"Reading VCF header from: {vcf_path}")

# --- Helper Functions ---
def get_chromosome_lists(vcf_path, exclude_file, report_csv, output_csv, plot_file):
    """
    Quickly reads the VCF header, returning a dictionary of chromosomes > 5Mb (ref_chroms),
    a dictionary of all chromosomes > 1Mb (long_chroms), and a list of chromosomes < 1Mb.
    """
    long_chroms = {}
    short_chroms = []
    ref_chroms = {}
    with pysam.VariantFile(vcf_path, 'r') as vcf_reader:
        for contig in vcf_reader.header.contigs.values():
            if contig.length: # Ensure length is not None
                if contig.length >= 5_000_000:
                    ref_chroms[contig.name] = contig.length
                if contig.length >= 1_000_000:
                    long_chroms[contig.name] = contig.length
                else:
                    short_chroms.append(contig.name)

    if not long_chroms:
        print("No chromosomes > 1Mb found in VCF header.")
        print("Writing empty output files and exiting successfully.")
        # Write empty exclude file (no chromosomes to exclude since none are large enough)
        with open(exclude_file, 'w') as f:
            pass  # Empty file
        # Write empty report
        pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj', 'chr_len']).to_csv(report_csv, index=False)
        # Write empty data output
        pd.DataFrame(columns=['Sample', 'Chromosome', 'Depth', 'Heterozygosity', 'Normalized_Depth', 'Normalized_Heterozygosity']).to_csv(output_csv, index=False)
        # Write empty plot file
        with open(plot_file, 'w') as f:
            pass  # Empty file
        sys.exit(0)

    if not ref_chroms:
        print("No chromosomes > 5Mb found in VCF header — assembly is too fragmented for sex chromosome normalization.")
        print("Writing empty output files and exiting successfully.")
        with open(exclude_file, 'w') as f:
            pass  # Empty file — no sex chromosomes identified, nothing to exclude
        pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj', 'chr_len']).to_csv(report_csv, index=False)
        pd.DataFrame(columns=['Sample', 'Chromosome', 'Depth', 'Heterozygosity', 'Normalized_Depth', 'Normalized_Heterozygosity']).to_csv(output_csv, index=False)
        with open(plot_file, 'w') as f:
            pass  # Empty file
        sys.exit(0)

    print(f"Found {len(ref_chroms)} reference chromosomes (> 5Mb) for normalization.")
    print(f"Found {len(long_chroms)} total chromosomes > 1Mb for analysis.")
    print(f"Found {len(short_chroms)} chromosomes < 1Mb to be excluded.")
    return ref_chroms, long_chroms, short_chroms

# --- MODIFICATION 1: Updated the function signature ---
def process_raw_data(data_path, all_chrom_lengths, ref_chrom_keys):
    """
    Processes raw data to calculate depth and heterozygosity.
    Normalization medians are calculated ONLY from reference chromosomes (>5Mb).
    The analysis is performed on all chromosomes > 1Mb.
    """
    print(f"Analyzing {len(all_chrom_lengths)} chromosomes > 1Mb")

    # Running totals per (chromosome, sample): memory grows with chromosomes x
    # samples, not with variants. Holding every genotype of a chromosome first
    # ran out of memory at 32 GB on 1-Layia (71 samples, 14 GB VCF).
    chromosomes_to_process = set(all_chrom_lengths.keys())
    het = re.compile(r'0[/|]1|1[/|]0')
    totals = {}  # (chrom, sample) -> [depth sum, sites with depth, het sites]

    with gzip.open(data_path, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 100_000_000 == 0:
                print(f"  Processed {line_num:,} lines...")
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 4 or parts[0] not in chromosomes_to_process:
                continue
            try:
                depth = int(parts[2]) if parts[2] != '.' else None
            except ValueError:
                continue
            if depth is None:
                continue
            t = totals.get((parts[0], parts[1]))
            if t is None:
                t = totals[(parts[0], parts[1])] = [0, 0, 0]
            t[0] += depth
            t[1] += 1
            if het.search(parts[3]):
                t[2] += 1

    print(f"Combining results from {len({c for c, _ in totals})} chromosomes...")

    if not totals:
        print("Warning: No valid data found in file")
        return pd.DataFrame()

    # Chromosomes in file order, samples sorted within each, as before.
    summary = pd.DataFrame(
        [(sample, depth_sum, n, n_het, chrom) for (chrom, sample), (depth_sum, n, n_het) in totals.items()],
        columns=['Sample', 'Total_Depth', 'Variant_Count', 'Het_Count', 'Chromosome'],
    )
    order = {c: i for i, c in enumerate(dict.fromkeys(summary['Chromosome']))}
    summary = (summary.assign(_o=summary['Chromosome'].map(order))
               .sort_values(['_o', 'Sample'], kind='stable').drop(columns='_o')
               .reset_index(drop=True))

    summary['Depth'] = summary['Total_Depth'] / summary['Variant_Count']
    summary['chr_len'] = summary['Chromosome'].map(all_chrom_lengths)
    summary.dropna(subset=['chr_len'], inplace=True)
    summary['Heterozygosity'] = summary['Het_Count'] / summary['chr_len']
    
    # --- MODIFICATION 2: Calculate medians using only the reference chromosomes ---
    print("Calculating normalization medians using reference chromosomes (>= 5Mb)...")
    
    # Create a temporary DataFrame containing only the reference chromosomes for median calculation
    ref_summary = summary[summary['Chromosome'].isin(ref_chrom_keys)]

    median_stats = ref_summary.groupby('Sample', observed=True).agg(
        Depth_SampleMedian=('Depth', 'median'),
        Heterozygosity_SampleMedian=('Heterozygosity', 'median')
    ).reset_index()
    
    # Merge the reference-based medians back into the main summary DataFrame
    summary = summary.merge(median_stats, on='Sample', how='left')

    bad_samples_mask = summary['Depth_SampleMedian'].isna() | summary['Heterozygosity_SampleMedian'].isna()
    if bad_samples_mask.any():
        bad_sample_names = summary.loc[bad_samples_mask, 'Sample'].unique()
        print(f"Warning: Removing {len(bad_sample_names)} sample(s) with failed median calculation: {list(bad_sample_names)}")
        summary = summary[~summary['Sample'].isin(bad_sample_names)]

    epsilon = 1e-9
    summary['Normalized_Depth'] = summary['Depth'] / (summary['Depth_SampleMedian'] + epsilon)
    summary['Normalized_Heterozygosity'] = summary['Heterozygosity'] / (summary['Heterozygosity_SampleMedian'] + epsilon)

    # Replace any remaining inf or nan values that might have been created
    summary.replace([np.inf, -np.inf], np.nan, inplace=True)

    final_df = summary[['Sample', 'Chromosome', 'Depth', 'Heterozygosity', 'Normalized_Depth', 'Normalized_Heterozygosity']]
    final_df.to_csv(output_csv, index=False)
    return final_df

# --- Main Execution ---
ref_chroms, long_chroms, small_chromosomes = get_chromosome_lists(vcf_path, exclude_file, report_csv, output_csv, plot_file)

# --- MODIFICATION 3: Pass the ref_chroms dictionary keys to the processing function ---
df = process_raw_data(raw_data_path, long_chroms, ref_chroms.keys())

if df.empty:
    pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj']).to_csv(report_csv, index=False)
    with open(exclude_file, 'w') as f:
        for chrom in sorted(small_chromosomes):
            f.write(f"{chrom}\n")
    sys.exit("The processed dataframe is empty. Wrote empty reports and exiting.")

df['Cluster'] = np.nan

# --- Identify Potential Sex Chromosomes (Clustering and Stats) ---
print("Identifying sex chromosome candidates...")
results = []
grouped_chroms = df.groupby('Chromosome', observed=True)

for chrom, group in grouped_chroms:
    group_copy = group.copy()
    group_copy.replace([np.inf, -np.inf], np.nan, inplace=True)
    group_copy.dropna(subset=['Normalized_Depth', 'Normalized_Heterozygosity'], inplace=True)

    if len(group_copy['Sample'].unique()) < 20:
        continue

    # Ensure no NaN values remain before KMeans
    if group_copy['Normalized_Depth'].isna().any() or len(group_copy) == 0:
        continue

    kmeans = KMeans(n_clusters=2, random_state=42, n_init='auto')
    group_copy['Cluster'] = kmeans.fit_predict(group_copy[['Normalized_Depth']])
    
    df.loc[group_copy.index, 'Cluster'] = group_copy['Cluster']
    
    cluster_counts = group_copy['Cluster'].value_counts()
    sil_score = silhouette_score(group_copy[['Normalized_Depth']], group_copy['Cluster'])

    if all(count >= 10 for count in cluster_counts) and sil_score > 0.80:
        cluster0_het = group_copy[group_copy['Cluster'] == 0]['Normalized_Heterozygosity']
        cluster1_het = group_copy[group_copy['Cluster'] == 1]['Normalized_Heterozygosity']
        
        if len(cluster0_het.dropna()) > 1 and len(cluster1_het.dropna()) > 1:
            t_stat, p_val = ttest_ind(cluster0_het, cluster1_het, equal_var=False, nan_policy='omit')
            if not np.isnan(p_val):
                results.append([chrom, p_val, sil_score])

# --- Statistical Analysis and Reporting ---
if results:
    results_df = pd.DataFrame(results, columns=['Chromosome', 'p_value', "sil_score"])
    results_df['p_adj'] = multipletests(results_df['p_value'], method='bonferroni')[1]
    significant_results = results_df[results_df['p_adj'] < 0.0001].sort_values(by='p_adj')
    
    # --- ADDED LINE ---
    # Map the chromosome lengths from the long_chroms dictionary to the results
    significant_results['chr_len'] = significant_results['Chromosome'].map(long_chroms)
    
    significant_results.to_csv(report_csv, index=False)
    significant_chromosomes = set(significant_results['Chromosome'])
else:
    # Updated to include the new column in the empty file case
    pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj', 'chr_len']).to_csv(report_csv, index=False)
    significant_chromosomes = set()

print(f"Found {len(significant_chromosomes)} significant chromosome(s): {significant_chromosomes or 'None'}")
df['Cluster'] = df['Cluster'].astype('category')

# --- Write Exclude File ---
chromosomes_to_exclude = set(small_chromosomes).union(significant_chromosomes)
print(f"Writing {len(chromosomes_to_exclude)} chromosomes to exclude file: {exclude_file}")
with open(exclude_file, 'w') as f:
    for chrom in sorted(list(chromosomes_to_exclude)):
        f.write(f"{chrom}\n")

# --- Plotting ---
print("Generating plots...")
for chrom in significant_chromosomes:
    subset = df[df['Chromosome'] == chrom]
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=subset, x='Normalized_Depth', y='Normalized_Heterozygosity',
        hue='Cluster', palette='coolwarm', alpha=0.8, s=50
    )
    plt.xlabel('Normalized Depth (vs. Autosomal Median)')
    plt.ylabel('Normalized Heterozygosity Rate')
    plt.title(f'Sex Chromosome Candidate: {chrom}')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Depth Cluster')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'plot_{chrom}_details.png'), dpi=300)
    plt.close()

if not df.empty:
    plt.figure(figsize=(12, 8))
    autosomes = df[~df['Chromosome'].isin(significant_chromosomes)]
    if not autosomes.empty:
        plt.scatter(autosomes['Normalized_Depth'], autosomes['Heterozygosity'], color='grey', alpha=0.3, label='Autosomes')

    palette = sns.color_palette('bright', n_colors=len(significant_chromosomes))
    for i, chrom in enumerate(sorted(list(significant_chromosomes))):
        subset = df[df['Chromosome'] == chrom]
        if not subset.empty:
            plt.scatter(subset['Normalized_Depth'], subset['Heterozygosity'], color=palette[i], alpha=0.9, label=chrom)

    plt.xlabel('Normalized Depth (vs. Autosomal Median)')
    plt.ylabel('Heterozygosity Rate (Het sites per base)')
    plt.title('Normalized Depth vs. Heterozygosity')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Chromosome', loc='best', markerscale=1.5)
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()

print("Script finished successfully.")