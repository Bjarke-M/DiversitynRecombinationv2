import allel
import numpy as np
import pandas as pd
from typing import Optional, Tuple

def analyze_windows_combined(
    vcf_path: str,
    bed_path: Optional[str] = None,
    window_size: Optional[int] = None,
    min_sites: int = 3
) -> pd.DataFrame:
    """
    Analyze multiple population genetic statistics in windows, using either BED file windows 
    or fixed-size windows. Combines diversity, singleton, and Tajima's D calculations.
    
    Parameters:
    -----------
    vcf_path : str
        Path to the VCF file
    bed_path : str, optional
        Path to BED file defining windows. If not provided, fixed-size windows will be used.
    window_size : int, optional
        Size of windows in base pairs when not using BED file
    min_sites : int, optional
        Minimum number of segregating sites required to compute Tajima's D
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing combined window-based statistics
    """
    # Read VCF file - only do this once
    callset = allel.read_vcf(vcf_path)
    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = callset['variants/POS']
    chroms = callset['variants/CHROM']
    
    results = []
    
    # Process each chromosome separately
    for chrom in np.unique(chroms):
        chrom_mask = chroms == chrom
        chrom_pos = pos[chrom_mask]
        chrom_gt = gt[chrom_mask]
        
        # Calculate allele counts once for this chromosome
        ac = chrom_gt.count_alleles()
        
        if bed_path:
            # Read BED file
            bed_df = pd.read_csv(bed_path, sep='\t', header=None,
                               names=['chrom', 'start', 'end'])
            # Filter BED regions for current chromosome
            chrom_bed = bed_df[bed_df['chrom'] == chrom]
            
            if len(chrom_bed) == 0:
                continue
                
            for _, region in chrom_bed.iterrows():
                # Get variants in this region - calculate once for all metrics
                region_mask = (chrom_pos >= region['start']) & (chrom_pos <= region['end'])
                region_pos = chrom_pos[region_mask]
                region_ac = ac[region_mask]
                region_size = region['end'] - region['start']
                
                if len(region_pos) == 0:
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': 0,
                        'diversity': 0,
                        'n_singletons': 0,
                        'singleton_proportion': 0,
                        'tajima_d': np.nan,
                        'bases': region_size
                    })
                    continue
                
                try:
                    # Calculate diversity
                    diversity, _, n_bases, n_variants = allel.windowed_diversity(
                        region_pos,
                        region_ac,
                        start=region['start'],
                        stop=region['end'],
                        windows=[(region['start'], region['end'])]
                    )
                    
                    # Calculate singleton statistics
                    is_singleton = region_ac.is_singleton(1)
                    rev_singleton = region_ac.is_singleton(0)
                    n_singletons = np.sum(is_singleton+rev_singleton)
                    singleton_proportion = np.mean(is_singleton+rev_singleton)
                    
                    # Calculate Tajima's D
                    tajima_d = allel.tajima_d(region_ac, min_sites=min_sites)
                    
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': n_variants[0],
                        'diversity': diversity[0],
                        'n_singletons': n_singletons,
                        'singleton_proportion': singleton_proportion,
                        'tajima_d': tajima_d,
                        'bases': n_bases[0]
                    })
                    
                except Exception as e:
                    print(f"Warning: Analysis failed for region {chrom}:{region['start']}-{region['end']}: {str(e)}")
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': len(region_pos),
                        'diversity': np.nan,
                        'n_singletons': np.nan,
                        'singleton_proportion': np.nan,
                        'tajima_d': np.nan,
                        'bases': region_size
                    })
                    
        else:
            # Use fixed-size windows
            try:
                # Calculate diversity
                diversity, windows, n_bases, n_variants = allel.windowed_diversity(
                    chrom_pos,
                    ac,
                    size=window_size,
                    start=chrom_pos[0],
                    stop=chrom_pos[-1]
                )
                
                # Calculate singleton statistics for the same windows
                is_singleton = ac.is_singleton(1)
                singleton_counts = allel.windowed_count(chrom_pos[is_singleton], 
                                                      size=window_size,
                                                      start=chrom_pos[0],
                                                      stop=chrom_pos[-1])[0]
                
                # Calculate singleton proportion
                singleton_proportion = np.divide(
                    singleton_counts,
                    n_variants,
                    out=np.zeros_like(singleton_counts, dtype=float),
                    where=n_variants > 0
                )
                
                
                # Calculate Tajima's D
                tajima_d, windows_td, counts_td = allel.windowed_tajima_d(
                    chrom_pos,
                    ac,
                    size=window_size,
                    start=chrom_pos[0],
                    stop=chrom_pos[-1],
                    min_sites=min_sites
                )
                
                # Collect results
                for i in range(len(windows)):
                    results.append({
                        'chrom': chrom,
                        'start': windows[i][0],
                        'end': windows[i][1],
                        'n_variants': n_variants[i],
                        'diversity': diversity[i],
                        'n_singletons': singleton_counts[i],
                        'singleton_proportion': singleton_proportion[i],
                        'tajima_d': tajima_d[i],
                        'bases': n_bases[i]
                    })
                    
            except Exception as e:
                print(f"Warning: Analysis failed for chromosome {chrom}: {str(e)}")
    
    return pd.DataFrame(results)

# Example usage

#Run the function
vcf = snakemake.input.vcf
bed = snakemake.input.bed
outfile = snakemake.output.csv
df = analyze_windows_combined(vcf, bed)
df.to_csv(outfile, index=False)

#tester
# vcf = 'data/vcfs/Alouatta_puruensis.vcf.gz'
# bed = 'data/beds/Alouatta_puruensis.bed'
# df = analyze_windows_combined(vcf, bed)
# df.to_csv('data/stats/Alouatta_puruensis.window.stats.csv', index=False)