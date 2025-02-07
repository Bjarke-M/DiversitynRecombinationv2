import allel
import numpy as np
import pandas as pd

def analyze_variants_in_windows(vcf_path, bed_path=None, window_size=None):
    """
    Analyze variants in windows, using either BED file windows or fixed-size windows.
    """
    callset = allel.read_vcf(vcf_path)
    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = callset['variants/POS']
    chroms = callset['variants/CHROM']
    ac = gt.count_alleles()
    
    results = []
    
    # Process each chromosome separately
    for chrom in np.unique(chroms):
        chrom_mask = chroms == chrom
        chrom_pos = pos[chrom_mask]
        chrom_ac = ac[chrom_mask]
        
        if bed_path:
            bed_df = pd.read_csv(bed_path, sep='\t', header=None,
                               names=['chrom', 'start', 'end'])
            # Filter BED regions for current chromosome
            chrom_bed = bed_df[bed_df['chrom'] == chrom]
            if len(chrom_bed) == 0:
                continue
                
            for _, region in chrom_bed.iterrows():
                # Get variants in this region
                region_mask = (chrom_pos >= region['start']) & (chrom_pos <= region['end'])
                region_pos = chrom_pos[region_mask]
                region_ac = chrom_ac[region_mask]
                
                if len(region_pos) == 0:
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': 0,
                        'n_singletons': 0,
                        'singleton_proportion': 0,
                        'mean_af': 0
                    })
                    continue
                
                # Calculate statistics for this region
                afs = region_ac.to_frequencies()
                mean_af = np.mean(afs[:, 1]) if afs.shape[1] > 1 else 0
                is_singleton = region_ac.is_singleton(1)
                
                results.append({
                    'chrom': chrom,
                    'start': region['start'],
                    'end': region['end'],
                    'n_variants': len(region_pos),
                    'n_singletons': np.sum(is_singleton),
                    'singleton_proportion': np.mean(is_singleton),
                    'mean_af': mean_af
                })
        else:
            # Use fixed-size windows
            counts, windows = allel.windowed_count(chrom_pos, size=window_size)
            
            # Calculate statistics for each window
            afs = chrom_ac.to_frequencies()
            mean_afs = allel.windowed_statistic(
                chrom_pos, afs[:, 1],
                statistic=lambda x: np.mean(x) if len(x) > 0 else 0,
                size=window_size
            )
            
            is_singleton = chrom_ac.is_singleton(1)
            singleton_counts = allel.windowed_count(chrom_pos[is_singleton], size=window_size)[0]
            
            # Calculate singleton proportion
            singleton_proportion = np.divide(
                singleton_counts, 
                counts, 
                out=np.zeros_like(singleton_counts, dtype=float), 
                where=counts > 0
            )
            
            for i in range(len(windows)):
                results.append({
                    'chrom': chrom,
                    'start': windows[i][0],
                    'end': windows[i][1],
                    'n_variants': counts[i],
                    'n_singletons': singleton_counts[i],
                    'singleton_proportion': singleton_proportion[i],
                    'mean_af': mean_afs[i]
                })
    
    return pd.DataFrame(results)

# Run the function
# vcf = snakemake.input.vcf
# bed = snakemake.input.bed
# df = analyze_variants_in_windows(vcf, bed)
# df.to_csv(snakemake.output.csv, index=False)

vcf= 'data/vcfs/Alouatta_puruensis_no_0429.vcf'
bed = 'data/vcfs/Alouatta_puruensis_100kb.bed'
df = analyze_variants_in_windows(vcf, bed)
df.to_csv(snakemake.output.csv, index=False)