import pandas as pd

mode, ps, sample = glob_wildcards("stringtie_{mode}/{ps}/{sample}.gtf")
ps2, sample2 = glob_wildcards("stringtie/{ps}/{sample}.gtf")

configfile: "config/config.yaml"
sample = pd.Series(sample)
remove = sample.str.contains('combined') + sample.str.contains('annotated') 
mode = pd.Series(mode)[~remove]
ps = pd.Series(ps)[~remove]
sample = sample[~remove]


rule all:
    input:
        expand(
          "stringtie_{mode}/{ps}/{sample}.stats",
          zip,
          mode=mode,
          ps=ps,
          sample=sample),
        expand(
            "stringtie/{ps}/{sample}.stats",
            zip,
            ps=ps2,
            sample=sample2)


rule gffcompare:
    input:
        "stringtie_{mode}/{ps}/{sample}.gtf",
    output:
        "stringtie_{mode}/{ps}/{sample}.stats"
    conda:
        "envs/gffcompare.yaml"
    params:
        ref=config["ref"],
        prefix=lambda wc, output: output[0][:-6],
    shell:
        "gffcompare -Q -r {params.ref} -o {params.prefix} {input}"

rule gffcompare2:
    input:
        "stringtie/{ps}/{sample}.gtf",
    output:
        "stringtie/{ps}/{sample}.stats"
    conda:
        "envs/gffcompare.yaml"
    params:
        ref=config["ref"],
        prefix=lambda wc, output: output[0][:-6],
    shell:
        "gffcompare -Q -r {params.ref} -o {params.prefix} {input}"
