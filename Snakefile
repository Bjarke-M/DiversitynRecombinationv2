import yaml

# Load the configuration file
configfile: "config.yaml"

# Access configuration values
species = config["species_mapping"].keys()
# Species with more than 5 samples
n=5
n_samples_species = [species for species in config["species_mapping"] if len(config["species_mapping"][species]['samples']) >= n]


######### Command to run the pipeline for the new data #########

# snakemake --executor slurm -j 200 --sdm conda --retries 3 --keep-going --default-resources slurm_account=primatediversity 

######### import the rules from the rules file #########

# TODO: Add rule for subsetting and preparing the gVCFs
### Rule for subsetting and preparing the gVCFs ###
#   Concatenate the different batches before liftover.
#   Prepare gVCFs for callability mask. 
#   Filter monomorphic sites for downstream analysis. 
#include: "rules/combine_batches_gvcfs.smk"

# TODO: Add rule for callability mask
### Rule for callability mask ###
#   Prepare the callability mask for the gVCFs.  
#include: "rules/callability_mask.smk"
# TODO:
### Rule for liftover ###
#   Liftover the gVCFs to the hg38 reference genome.
#   Liftover the gVCFs to the hg38 reference genome using Crossmap.
#include: "rules/liftover.smk"
# TODO:
### Rule for degeneracy annotation ###
#   Annotate the degeneracy of the gVCFs.
#include: "rules/annotate_degeneracy.smk"

# TODO: GET RECOMB RATES AND MERGE PER SPECIES AND LATER MERGE ALL
### Rule for removing bad samples ###
#   Remove samples with low coverage.
include: "rules/remove_bad_samples.smk"

### Rule for getting windowbased summary statistics ###
#   Get summary statistics for the VCFs.
#   Estimate windowbased Pi and Tajima's D.
include: "rules/window_stats.smk"

### Rule for getting SFS ###
include: "rules/get_sfs.smk"

### Rule for Merging stat files ###
#   Merge the windowbased summary statistics.
include: "rules/merge_stats.smk"

# Define the pipeline
rule all:
    input:
        expand("data/vcfs/{species}.masked.biallelic.vcf.gz", species=species),
        expand("data/vcfs/{species}.masked.biallelic.vcf.gz.tbi", species=species),
        #expand('data/stats/{species}.window.stats.csv',species=species),
        #expand('data/stats/{species}.combined.stats.csv',species=species),
        #"data/stats/all_species.combined.stats.csv",
        #expand('data/callable_fraction/{species}_positive_mask.bed',species='Gorilla_gorilla'),
        #expand('data/vcfs/{species}.masked.biallelic.vcf.gz',species='Gorilla_gorilla'),
        expand('results/stats/sfs/{species}.sfs.csv',species=n_samples_species),
        #'results/stats/sfs/Aotus_vociferans.sfs.csv'
ruleorder:
    filter_bcf_samples >
    filter_multiallelic_sites_nonvariant_sites_and_missing_genotypes >
    window_stats >
    positive_mask >
    get_sfs >
    merge_per_species >
    merge_all >
    index 


rule index:
    input:
        "{path}.vcf.gz"
    output:
        "{path}.vcf.gz.tbi"
    resources:
        mem_mb= 1000*8,
        runtime= 60*4
    conda:
        "env/tabix.yaml"
    shell:
        "tabix -p vcf {input}"