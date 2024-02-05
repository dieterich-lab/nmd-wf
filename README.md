# NMD workflow

[![License](https://img.shields.io/github/license/dieterich-lab/nmd-wf)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/dieterich-lab/nmd-wf)](https://github.com/dieterich-lab/nmd-wf/issues)

This repository exists for reproducibility purpose. The data generated on this workflow powers the [NMDtxDB](https://github.com/dieterich-lab/nmdtxdb). Raw data is available at the SRA [PRJNA1054031](https://www.ncbi.nlm.nih.gov/sra/PRJNA1054031). RNA-seq reads need to be pre-processed and aligned. 

## Workflow description

The workflow comprises two parts. The first part comprises a Snakemake workflow (`workflow`). The second part enables the CDS detection and integration. 

## Usage

### Part 1 

This refers to the workflow to generate the de novo transcriptome, and compute DGE and DTE.

```{bash}
snakemake --jobs 10 --cores 10 --profile slurm --printshellcmds --reason --use-singularity --use-conda --use-envmodule
```

To produce the DAG:
```{bash}
snakemake --rulegraph | dot -Tsvg > rulegraph.sv
```

### Part 2
This refers to the workflow for CDS detection. Here an example using sequences trimmed by the Ensembl start codon:


```{bash}

awk '{ print $1 "\t" $7-1 "\t" $8 "\t" $4 "\t" 1 "\t" $6; }' GRCh38.102.gtf > ref_cds.bed

Rscript cds/StartATG_to_cDNA.R ref_cds.bed

perl longorf2_fwd_v2.pl --input GRCh38.102.fa --startcodon ref_cds_cDNA.bed > ensembl_longorf2.fa 
```

See [longorf_integration_bed12](cds/longorf_integration_bed12.R) script, which details how the multiple source integration is done. 

## License

This project is licensed under the [MIT](LICENSE). 

## Funding
This work was supported by the DFG Research Infrastructure West German Genome Center, project 407493903, as part of the Next-Generation Sequencing Competence Network, project 423957469.

