import allel
import numpy as np
import pandas as pd

def analyze_pairwise_diff_in_windows(vcf_path, bed_path=None, window_size=None):
    """
    Analyze mean pairwise differences in windows, using either BED file windows or fixed-size windows.
    
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
        DataFrame containing window-based statistics
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
                region_gt = chrom_gt[region_mask]
                
                if len(region_pos) == 0:
                    results.append({
                        'chrom': chrom,
                        'start': region['start'],
                        'end': region['end'],
                        'n_variants': 0,
                        'mean_pairwise_diff': 0,
                        'variance_pairwise_diff': 0
                    })
                    continue
                
                # Calculate mean pairwise difference
                ac = region_gt.count_alleles()
                mpd = allel.mean_pairwise_difference(ac)
                
                # Calculate variance of pairwise differences
                # First get all pairwise differences
                n_samples = region_gt.n_samples
                pairwise_diffs = []
                
                for i in range(n_samples):
                    for j in range(i + 1, n_samples):
                        # Get genotypes for each sample pair
                        gt_i = region_gt[:, i]
                        gt_j = region_gt[:, j]
                        # Count differences
                        diff = np.sum(gt_i != gt_j, axis=1)
                        pairwise_diffs.extend(diff)
                
                var_pd = np.var(pairwise_diffs) if pairwise_diffs else 0
                
                results.append({
                    'chrom': chrom,
                    'start': region['start'],
                    'end': region['end'],
                    'n_variants': len(region_pos),
                    'mean_pairwise_diff': mpd,
                    'variance_pairwise_diff': var_pd
                })
        else:
            # Use fixed-size windows
            windows = allel.position_windows(chrom_pos, size=window_size)
            
            for start, stop in windows:
                # Get variants in this window
                window_mask = (chrom_pos >= start) & (chrom_pos < stop)
                window_gt = chrom_gt[window_mask]
                
                if window_gt.shape[0] == 0:
                    results.append({
                        'chrom': chrom,
                        'start': start,
                        'end': stop,
                        'n_variants': 0,
                        'mean_pairwise_diff': 0,
                        'variance_pairwise_diff': 0
                    })
                    continue
                
                # Calculate mean pairwise difference
                ac = window_gt.count_alleles()
                mpd = allel.mean_pairwise_difference(ac)
                
                # Calculate variance of pairwise differences
                n_samples = window_gt.n_samples
                pairwise_diffs = []
                
                for i in range(n_samples):
                    for j in range(i + 1, n_samples):
                        gt_i = window_gt[:, i]
                        gt_j = window_gt[:, j]
                        diff = np.sum(gt_i != gt_j, axis=1)
                        pairwise_diffs.extend(diff)
                
                var_pd = np.var(pairwise_diffs) if pairwise_diffs else 0
                
                results.append({
                    'chrom': chrom,
                    'start': start,
                    'end': stop,
                    'n_variants': window_gt.shape[0],
                    'mean_pairwise_diff': mpd,
                    'variance_pairwise_diff': var_pd
                })
    
    return pd.DataFrame(results)

# Run the function
vcf = snakemake.input.vcf
bed = snakemake.input.bed
df = analyze_variants_in_windows(vcf, bed)
df.to_csv(snakemake.output.csv, index=False)