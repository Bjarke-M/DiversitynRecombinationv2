import pandas as pd

input_file = snakemake.input.fraction_callable
accepted_fraction = snakemake.params.accepted_fraction
outfile = snakemake.output.positive_mask

# read file
df = pd.read_csv(input_file, sep='\t')
# remove regions of low callability
df = df[df['freq_median']>accepted_fraction]
bed = pd.DataFrame()
bed['CHROM'] = df['chr']
bed['BEG'] = df['start']
bed['END'] = df['end']
bed.to_csv(outfile,sep='\t', index=False, header=False)