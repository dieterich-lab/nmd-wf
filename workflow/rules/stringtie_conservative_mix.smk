
strand = {"fr-firststrand": "--fr", "fr-secondstrand": "--rf"}

rule stringtie_conservative_mix:
    input:
        short_bam="{group}.bam",
        long_bam=lambda x: ...
    output:
        "stringtie_conservative_mix/{group}.gtf",
    params:
        strandness=strand.get(config.get("strandness", ""), ""),
        min_junct_coverage=config.get("min_junct_coverage", 3),
        min_isoform_proportion=config.get("min_isoform_proportion", 0.001),
        minimum_read_per_bp_coverage=config.get("minimum_read_per_bp_coverage", 3),
    conda:
    	"envs/stringtie_2.15.yaml",
    log:
        "logs/stringtie_conservative_mix/{group}.log",
    shell:
        "stringtie {input.short_bam} -o {output} "
        "-p {threads} "
        "--mix {input.short_bam} "
        "{params.strandness} "
        "-c {params.minimum_read_per_bp_coverage} "
        "-j {params.min_junct_coverage} "
        "-f {params.min_isoform_proportion} "
        "2> {log} "
