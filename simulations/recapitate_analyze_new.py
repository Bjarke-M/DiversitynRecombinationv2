#!/usr/bin/env python3
import sys
import tskit
import msprime
import pyslim
import numpy as np
import pandas as pd
from pandas import IntervalIndex

# Parse arguments
scenario = sys.argv[1]
ne = int(sys.argv[2])
s = float(sys.argv[3]) if sys.argv[3] != 'NA' else None
c = int(sys.argv[4])
mu = float(sys.argv[5])
mu_key = str(sys.argv[6])
recomb_ends_file = sys.argv[7]
recomb_rates_file = sys.argv[8]
output_file = sys.argv[9]
input_files = sys.argv[10:]

print(input_files)


# Load recombination map
ends = pd.read_csv(recomb_ends_file, names=['ends'], index_col=False)['ends'].tolist()
rates = pd.read_csv(recomb_rates_file, names=['rates'], index_col=False)['rates'].tolist()
recomb_map = msprime.RateMap(position=ends, rate=rates)

# Create recombination rate dataframe for mapping
rates_df = pd.DataFrame({
    'left': recomb_map.position[:-1].astype(int),
    'right': recomb_map.position[1:].astype(int),
    'rate': recomb_map.rate*c #plus c
})
rate_intervals = IntervalIndex.from_arrays(rates_df['left'], rates_df['right'], closed='both')


df_list = []
for replicate in input_files:
    ts = tskit.load(replicate)
    # Recapitate
    recap = pyslim.recapitate(ts, recombination_rate=recomb_map, ancestral_Ne=ne, record_provenance=False)
    # Simplify and add mutations
    tree = recap.simplify()
    next_id = pyslim.next_slim_mutation_id(tree)
    mts = msprime.sim_mutations(tree, rate=1e-7, model=msprime.SLiMMutationModel(type=0, next_id=next_id), keep=False)

    # Calculate diversity in windows
    windows = np.linspace(0, mts.sequence_length, 1300).astype(int)
    pi = mts.diversity(mts.samples(), windows=windows)
    taj_d = mts.Tajimas_D(mts.samples(), windows=windows)
    seg_sites = mts.segregating_sites(sample_sets=mts.samples(), windows=windows)

    df = pd.DataFrame({
    'left': windows[:-1],
    'right': windows[1:],
    'pi': pi,
    'taj_d': taj_d,
    'seg_sites': seg_sites
})
    df['mid'] = df['left'] + ((df['right'] - df['left']) // 2)
    df['span'] = df['right'] - df['left']

    # Map recombination rates to windows
    weighted_rates = []
    for idx, window in df.iterrows():
        overlap_mask = rate_intervals.overlaps(pd.Interval(window['left'], window['right'], closed='both'))
        overlapping_rates = rates_df[overlap_mask].copy()

        if len(overlapping_rates) == 0:
            weighted_rates.append(np.nan)
        else:
            overlap_starts = np.maximum(overlapping_rates['left'].values, window['left'])
            overlap_ends = np.minimum(overlapping_rates['right'].values, window['right'])
            overlap_lengths = overlap_ends - overlap_starts
            weights = overlap_lengths / overlap_lengths.sum()
            weighted_rate = (overlapping_rates['rate'].values * weights).sum()
            weighted_rates.append(weighted_rate)

    df['rate'] = weighted_rates
    df['replicate'] = replicate
    df_list.append(df)

# Merge all replicates
merged_df = pd.concat(df_list, ignore_index=True)

 # Bin by recombination rate and calculate median and variance
merged_df['bin'] = pd.qcut(merged_df['rate'], q=20)
grouped_df = merged_df.groupby('bin', observed=False).agg({
    'pi': ['median', 'var'],
    'taj_d': 'median',
    'seg_sites': 'median',
    'rate': 'median'
})
grouped_df.columns = ['_'.join(col).strip('_') for col in grouped_df.columns.values]
grouped_df = grouped_df.rename(columns={'pi_median': 'pi', 'taj_d_median': 'taj_d', 'seg_sites_median': 'seg_sites', 'rate_median': 'rate'})
grouped_df['cm_mb'] = grouped_df['rate'] * 1e8 * c

# Add metadata columns
grouped_df['scenario'] = scenario
grouped_df['ne'] = ne
grouped_df['s'] = s if s is not None else 'NA'
grouped_df['mu'] = mu
grouped_df['mu_key'] = mu_key
grouped_df['c'] = c


# Save
grouped_df.to_csv(output_file)