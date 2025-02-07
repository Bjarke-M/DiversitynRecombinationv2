rule mean_pairwise_difference:
    input:
        vcf = "data/vcfs/{species}.vcf.gz",
        tbi = "data/vcfs/{species}.vcf.gz.tbi",
        bed = "data/beds/{species}.bed"
    output:
        "data/stats/{species}.pi.csv"
    conda:
        "../env/base.yaml"
    script:
        "../scripts/window_mean_pairwise_diff.py"