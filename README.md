# NMD workflow


[![License](https://img.shields.io/github/license/dieterich-lab/nmd-wf)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/dieterich-lab/nmd-wf)](https://github.com/dieterich-lab/nmd-wf/issues)

This repository exists for reproducibility porpusese. The data generated on this workflow powers the [NMDtxDB](https://github.com/dieterich-lab/nmdtxdb). Raw data is avaiable at the SRA [PRJNA1054031](https://www.ncbi.nlm.nih.gov/sra/PRJNA1054031). RNA-seq reads need to be pre-processed and aligment before. 

## Workflow description

The workflow comprises of two parts. The first part comprises a snakemake workflow (`workflow`). The second part enables the CDS detection and integration. 

## Usage

```
snakemake --jobs 10 --cores 10 --profile slurm --printshellcmds --reason --use-singularity --use-conda --use-envmodule

```

To produce the DAG:
```
snakemake --rulegraph | dot -Tsvg > rulegraph.sv
```

## License

This project is licensed under the [MIT](LICENSE). 

## Funding
This work was supported by the DFG Research Infrastructure West German Genome Center, project 407493903, as part of the Next-Generation Sequencing Competence Network, project 423957469.