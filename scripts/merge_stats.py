import pandas as pd

rec_rates = snakemake.input.rec_rates #tsv
window_stats = snakemake.input.window_stats #csv
callability_mask = snakemake.input.callability_mask #tsv
outfile = snakemake.output.combined_stats #csv


# example files
# rec_rates = 'data/rec_rates/Allenopithecus_nigroviridis_rec_rates.tsv'
# window_stats = 'data/stats/Allenopithecus_nigroviridis.window.stats.csv'
# callability_mask = 'data/callable_fraction/Allenopithecus_nigroviridis.callable_fraction.csv'
# outfile = 'test_combined_stats.csv'

# Read the dataframes from files
df1 = pd.read_csv(window_stats, sep=',') # CSV file
df2 = pd.read_csv(rec_rates, sep='\t')  # TSV file
df3 = pd.read_csv(callability_mask, sep='\t')  # TSV file

# Rename the chromosome column in df2 and df3 to match df1
df2 = df2.rename(columns={'chr': 'chrom'})
df3 = df3.rename(columns={'chr': 'chrom'})

# Convert start and end columns to integers for all dataframes
for df in [df1, df2, df3]:
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)

# Merge the dataframes
merged_df = df1.merge(df2, on=['chrom', 'start', 'end'], how='outer')
final_merged_df = merged_df.merge(df3, on=['chrom', 'start', 'end'], how='outer')

# Sort the dataframe by chromosome and start position
final_merged_df = final_merged_df.sort_values(['chrom', 'start'])

# Save to CSV if needed
final_merged_df.to_csv(outfile, index=False)