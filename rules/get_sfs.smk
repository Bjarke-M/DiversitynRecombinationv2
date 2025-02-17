rule get_sfs:
    input:
        vcf = 'data/vcfs/{species}.masked.biallelic.vcf.gz',
        tbi = 'data/vcfs/{species}.masked.biallelic.vcf.gz.tbi'
    output:
        data = 'results/stats/sfs/{species}.sfs.csv',
        plot = 'results/stats/sfs/{species}.sfs.pdf'
    conda:
        '../env/base.yaml'
    resources:
        mem_mb = 1000*86,
        runtime = 60*2,
    script:
        '../scripts/sfs.py'


