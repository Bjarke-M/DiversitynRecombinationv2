def get_time(wildcards, attempt):
    return attempt * 30


rule merge_per_species:
    input:
        rec_rates = "data/rec_rates/{species}_rec_rates.tsv",
        window_stats = "data/stats/{species}.window.stats.csv",
        callability_mask = "data/callable_fraction/{species}.callable_fraction.csv", # actually a tsv.. ups
    output:
        combined_stats = "data/stats/{species}.combined.stats.csv",
    resources:
        mem_mb = 1000*4,
        runtime=get_time,
    conda:
        "../env/base.yaml"
    script:
        "../scripts/merge_stats.py"

rule merge_all:
    input:
        file_list = expand("data/stats/{species}.combined.stats.csv", species=species),
    output:
        all_stats= "data/stats/all_species.combined.stats.csv",
    resources:
        mem_mb = 1000*12,
        runtime=60,
    conda:
        "../env/base.yaml"
    script:
        "../scripts/merge_all.py"

# example files
# rec_rates = 'data/rec_rates/Allenopithecus_nigroviridis_rec_rates.tsv'
# window_stats = 'data/stats/Allenopithecus_nigroviridis.window.stats.csv'
# callability_mask = 'data/callable_fraction/Allenopithecus_nigroviridis.callable_fraction.csv'
# outfile = 'test_combined_stats.csv'