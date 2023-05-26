
rule stringtie_merge_bam:
    input:
        lambda wc: ["mappings/{}_{}.bam".format(*x) for x in d[wc.group]],
    output:
        bam="stringtie/merged_bam/{group}.bam",
    threads: 10
    log:
        "logs/stringtie_merge_bam/{group}.log",
    shell:
        "samtools merge {output.bam} {input} "
        "--threads {threads} | "
        "samtools index {output.bam} {output.bai} "
        "2> {log} "
