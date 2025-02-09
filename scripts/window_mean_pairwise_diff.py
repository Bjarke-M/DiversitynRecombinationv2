import allel
import numpy as np
import pandas as pd

def analyze_diversity_in_windows(vcf_path, bed_path=None, window_size=None):
    """
    Analyze nucleotide diversity in windows, using either BED file windows or fixed-size windows.
    
    Parameters:
    -----------
    vcf_path : str
        Path to the VCF file
    bed_path : str, optional
        Path to BED file defining windows. If not provided, fixed-size windows will be used.
    window_size : int, optional
        Size of windows in base pairs when not using BED file
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing window-based diversity statistics
    """
    # Read VCF file
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
        
        # Calculate allele counts for the chromosome
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
                # Get variants in this region
                region_mask = (chrom_pos >= region['start']) & (chrom_pos <= region['end'])
                region_pos = chrom_pos[region_mask]
                region_ac = ac[region_mask]
                
                if len(region_pos) == 0:
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': 0,
                        'diversity': 0,
                        'bases': region['end'] - region['start']
                    })
                    continue
                
                try:
                    # Calculate diversity using windowed_diversity
                    diversity, _, n_bases, n_variants = allel.windowed_diversity(
                        region_pos, 
                        region_ac,
                        start=region['start'],
                        stop=region['end'],
                        windows=[(region['start'], region['end'])]
                    )
                    
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': n_variants[0],
                        'diversity': diversity[0],
                        'bases': n_bases[0]
                    })
                    
                except Exception as e:
                    print(f"Warning: Diversity calculation failed for region {chrom}:{region['start']}-{region['end']}: {str(e)}")
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': len(region_pos),
                        'diversity': np.nan,
                        'bases': region['end'] - region['start']
                    })
        else:
            # Use fixed-size windows
            try:
                # Calculate diversity using windowed_diversity
                diversity, windows, n_bases, n_variants = allel.windowed_diversity(
                    chrom_pos,
                    ac,
                    size=window_size,
                    start=chrom_pos[0],
                    stop=chrom_pos[-1]
                )
                
                # Collect results
                for i in range(len(windows)):
                    results.append({
                        'chrom': chrom,
                        'start': windows[i][0],
                        'end': windows[i][1],
                        'n_variants': n_variants[i],
                        'diversity': diversity[i],
                        'bases': n_bases[i]
                    })
                    
            except Exception as e:
                print(f"Warning: Diversity calculation failed for chromosome {chrom}: {str(e)}")
    
    return pd.DataFrame(results)

# Example usage

# Run the function
# vcf = snakemake.input.vcf
# bed = snakemake.input.bed
# df = analyze_variants_in_windows(vcf, bed)
vcf = 'data/vcfs/Alouatta_puruensis.vcf.gz'
bed = 'data/beds/Alouatta_puruensis.bed'
df = analyze_diversity_in_windows(vcf, bed)
df.to_csv('data/stats/Alouatta_puruensis.pi2.csv', index=False)