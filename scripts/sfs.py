import allel
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

vcf = snakemake.input.vcf
sfs_out = snakemake.output.data
plot_out = snakemake.output.plot

# example
# vcf = 'data/vcfs/Gorilla_gorilla.masked.biallelic.vcf.gz'
# sfs_out = 'results/stats/sfs/Gorilla_gorilla.sfs.csv'
# plot_out = 'results/stats/sfs/Gorilla_gorilla.sfs.pdf'

# Read VCF file - only do this once
callset = allel.read_vcf(vcf)
gt_raw = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']
chroms = callset['variants/CHROM']

#filter sites with a missing genotype
missing_sites = np.any(gt_raw.is_missing(),axis=1) # returns a vector with len = #sites and true if site has a missing genotype
gt = gt_raw[~missing_sites] # reverses the bool to remove sites that have missing genotype

print(len(gt_raw),len(gt))
print(gt[:5])

ac = gt.count_alleles()
folded = allel.sfs_folded(ac)


# Save SFS
df = pd.DataFrame(folded)
df.to_csv(sfs_out, index = False)


# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot as histogram instead of line
x = np.arange(len(folded))
ax.bar(x, folded, width=0.8, align='center', alpha=0.7)

# Customize the plot
ax.set_xlabel('Minor allele count')
ax.set_ylabel('Count')
ax.set_title('Folded Site Frequency Spectrum')

# Make x-axis ticks match the actual counts
ax.set_xticks(x)

# Save to PDF
plt.savefig(plot_out, format='pdf', bbox_inches='tight', dpi=300)
plt.close()
