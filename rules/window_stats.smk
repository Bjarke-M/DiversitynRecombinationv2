def get_mem_mb(wildcards, attempt):
    return attempt * 1000 * 32

def get_time(wildcards, attempt):
    return attempt * 60 *2


rule window_stats:
    input:
        vcf = "data/vcfs/{species}.masked.biallelic.vcf.gz",
        tbi = "data/vcfs/{species}.masked.biallelic.vcf.gz.tbi",
        bed = "data/beds/{species}.bed"
    output:
        csv ="data/stats/{species}.window.stats.csv"
    resources:
        mem_mb=get_mem_mb,
        runtime= get_time
    conda:
        "../env/base.yaml"
    script:
        "../scripts/window_stats.py"