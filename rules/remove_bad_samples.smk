def get_species_samples(wildcards):
    species = wildcards.species
    return config["species_mapping"][species]["samples"]



rule filter_bcf_samples:
    input:
        bcf="data/bcfs/{species}.bcf.gz"
    output:
        vcf="data/vcfs/{species}.vcf.gz"
    params:
        # Join sample IDs with newlines for direct use in shell command
        samples=lambda w: "\\n".join(get_species_samples(w))
    resources:
        mem_mb=1000*8,
        runtime=60*4
    conda:
        "../env/bcftools.yaml"
    shell:
        # Use echo and pipe to provide sample names directly to bcftools
        """
        echo -e '{params.samples}' | bcftools view \
            --samples-file - \
            --force-samples \
            {input.bcf} \
            -Oz \
            -o {output.vcf} \
        """


rule positive_mask:
    input:
        fraction_callable = 'data/callable_fraction/{species}.callable_fraction.csv'
    params:
        accepted_fraction = 0.5,
    output:
        positive_mask= 'data/callable_fraction/{species}_positive_mask.bed'
    conda:
        '../env/base.yaml'
    resources:
        mem_mb= 1000*4,
        runtime = 30
    script:
        '../scripts/positive_mask.py'


# rule mask_regions:
#     input:
#         vcf = 'data/vcfs/{species}.vcf.gz',
#         mask_bed = 'data/callable_fraction/{species}_positive_mask.bed'
#     output:
#         masked_vcf = 'data/vcfs/{species}.masked.vcf.gz'
#     resources:
#         mem_mb = 1000*5,
#         runtime = 60*4
#     conda:
#         '../env/bcftools.yaml'
#     shell:
#         """
#         bcftools view -R {input.mask_bed} {input.vcf} -Oz -o {output.masked_vcf}
#         """

rule filter_multiallelic_sites_nonvariant_sites_and_missing_genotypes:
    input:
        vcf = 'data/vcfs/{species}.vcf.gz',
        tbi = 'data/vcfs/{species}.vcf.gz.tbi',
        mask_bed = 'data/callable_fraction/{species}_positive_mask.bed'
    output:
        biallelic = 'data/vcfs/{species}.masked.biallelic.vcf.gz'
    resources:
        mem_mb = 1000*5,
        runtime = 60*6
    conda:
        '../env/bcftools.yaml'
    shell:
        """
        bcftools view -R {input.mask_bed} {input.vcf} -Ou | bcftools view -m2 -M2 -e 'GT[*] = "mis" || MAC == 0' -Ou | bcftools sort -Oz -o {output.biallelic}
        """

#bcftools view -m2 -M2 {input.masked_vcf} -Ou | bcftools filter -e "MAC == 0" -Ou | bcftools view -e 'GT[*] = "mis"' Ou |