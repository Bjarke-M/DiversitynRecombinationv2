def get_species_samples(wildcards):
    species = wildcards.species
    print(config["species_mapping"][species]["samples"])
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
        runtime=60
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