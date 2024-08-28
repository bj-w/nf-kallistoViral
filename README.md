# nf-kallistoViral
This workflow is a wrapper around the method described by Luebbert et al., 2024<sup>[1](https://www.biorxiv.org/content/10.1101/2023.12.11.571168v2)</sup> in which Kallisto is used to quantify reads of viral origin in RNA-seq data (bulk + since-cell). 
This wrapper is currently **bulk only**.

---

## Requirements
- [Nextflow](https://github.com/nextflow-io/nextflow) (tested with v24.04.4.5917)
- [Singularity](https://github.com/sylabs/singularity) (tested with v3.11.5-1.el7)

This workflow uses a containerised version of [kallisto | bustools](https://github.com/pachterlab/kallistobustools/) v0.28.2<sup>[2](https://www.nature.com/articles/s41587-021-00870-2)</sup> running kallisto v0.50.1<sup>[3](https://www.nature.com/articles/nbt.3519)</sup> and bustools v0.43.2<sup>[2](https://www.nature.com/articles/s41587-021-00870-2)</sup>.

---
## Set up
This workflow requires several reference files providing information on both viral and 'host' species to be downloaded prior to running. Here, 'host' species refers to the species of origin of the sequenced samples.

Firstly, the reference genome and transcriptome (cDNA) of the host species is required in FASTA format, along with a GTF annotation. 
An example for the GENCODE v46 release for Humans is shown below.
```bash
# Genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
# cDNA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
# GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
```

Also required are viral reference files (modified PalmDB files) which can be found [here](https://github.com/pachterlab/LSCHWCP_2023/tree/main/PalmDB) and downloaded as shown below
```bash
# fasta
wget https://raw.githubusercontent.com/pachterlab/LSCHWCP_2023/main/PalmDB/palmdb_rdrp_seqs.fa
# T2G
wget https://raw.githubusercontent.com/pachterlab/LSCHWCP_2023/main/PalmDB/palmdb_clustered_t2g.txt
```

---
## Instructions
This workflow is written in Nextflow - see [here](https://www.nextflow.io/docs/latest/install.html) for installation info.
Alternatively, Nextflow can be installed into a conda/mamba environments (my preferred method).
```bash
micromamba create -n nextflow nextflow
```

Nextflow configurations are important for optimal functioning - more info [here](https://www.nextflow.io/docs/latest/config.html).
An example config file for use with Newcastle University's Rocket HPC is provided - `rocket.config`


### Workflow parameters
This workflow requires essential parameters described in the table below.
- `${projectDir}` describes the directory where kallistoViral.nf is located
- FASTQ files should be concatenated across lanes if necessary
	- In the example above, the FASTQ files are in the format:
		- SampleA_R1.fastq.gz, SampleA_R2.fastq.gz, 
		- SampleB_R1.fastq.gz, SampleB_R2.fastq.gz

| Parameter    | Description                                                                                     | Example                                                       |
|--------------|-------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| fastq        | Path to FASTQ files                                                                             | "${projectDir}/data/fastq/*_{R1,R2}.fastq.gz"                 |
| outdir       | Name of output folder to be created in the current directory                                    | output                                                        |
| strandedness | Whether RNA-seq library preparation retains strand information <forward / reverse / unstranded> | reverse                                                       |
| readType     | Whether sequencing generated single- or paired-end reads <single / paired>                      | paired                                                        |
| viral_fasta  | Path to viral FASTA file                                                                        | "${projectDir}/data/ref/palmdb_rdrp_seqs.fa"                  |
| viral_t2g    | Path to viral T2G file                                                                          | "${projectDir}/data/ref/palmdb_clustered_t2g.txt"             |
| cdna         | Path to host transcriptome reference.                                                         | "${projectDir}/data/ref/gencode.v46.transcripts.fa.gz"        |
| genome       | Path to host genomic reference                                                                | "${projectDir}/data/ref/GRCh38.primary_assembly.genome.fa.gz" |
| gtf          | Path to GTF annotation.                                                                          | "${projectDir}/data/ref/gencode.v46.annotation.gtf.gz"        |


The above parameters should be contained with a parameters file. An example is given in `params.json` for the file structure shown below.

```
.
├── data
│   ├── fastq
│   │   ├── SampleA_R1.fastq.gz
│   │   ├── SampleA_R2.fastq.gz
│   │   ├── SampleB_R1.fastq.gz
│   │   └── SampleB_R2.fastq.gz
│   └── ref
│       ├── gencode.v46.annotation.gtf.gz
│       ├── gencode.v46.transcripts.fa.gz
│       ├── GRCh38.primary_assembly.genome.fa.gz
│       ├── palmdb_clustered_t2g.txt
│       └── palmdb_rdrp_seqs.fa
├── kallistoViral.nf
├── launcher_slurm.sh
├── params.json
├── README.md
└── rocket.config
```


### Run workflow
The basic command for running the workflow is
```
micromamba activate nextflow
nextflow run kallistoViral.nf -config <config_file> -params-file <params_file>
```

However, it is advised to submit the nextflow job as a batch job to save resources on the login/head node. An example of this for use with the SLURM workload manager manager is given in `launcher_slurm.sh` which can be run using
```bash
sbatch launcher_slurm.sh
```

---
## Output
When finished successfully, the specified output folder should contain host (transcript + gene) and viral abundances for each sample. 
Information regarding the run + quantification is contained within JSON files.
```
output
├── host
│   ├── SampleA
│   │   ├── inspect.json
│   │   ├── kb_info.json
│   │   ├── quant_unfiltered
│   │   │   ├── abundance_1.tsv
│   │   │   └── abundance.gene_1.tsv
│   │   └── run_info.json
│   └── SampleB
│       ├── inspect.json
│       ├── kb_info.json
│       ├── quant_unfiltered
│       │   ├── abundance_1.tsv
│       │   └── abundance.gene_1.tsv
│       └── run_info.json
└── viral
    ├── SampleA
    │   ├── inspect.json
    │   ├── kb_info.json
    │   ├── quant_unfiltered
    │   │   ├── abundance_1.tsv
    │   │   ├── abundance_2.tsv
    │   │   ├── abundance.gene_1.tsv
    │   │   └── abundance.gene_2.tsv
    │   └── run_info.json
    └── SampleB
        ├── inspect.json
        ├── kb_info.json
        ├── quant_unfiltered
        │   ├── abundance_1.tsv
        │   ├── abundance_2.tsv
        │   ├── abundance.gene_1.tsv
        │   └── abundance.gene_2.tsv
        └── run_info.json
```

Paired-end sequencing produces two abundance files when quantifying viral reads as the use of amino acid sequences in pseudoalignment (`kb count -aa`) is not currently supported in paired-end reads.

Mapping the viral IDs to their respective taxonomies can be performed using [this file](https://github.com/pachterlab/LSCHWCP_2023/blob/main/PalmDB/ID_to_taxonomy_mapping.csv) for downstream analysis.


---
## To do
- Parse JSON log outputs
- Automatically map viral IDs to taxonomies
- Produce formatted count matrices
- Develop single cell workflow wrapper

---
## References
<sup>[1](https://www.biorxiv.org/content/10.1101/2023.12.11.571168v2)</sup> Luebbert, L., Sullivan, D.K., Carilli, M., Hjörleifsson, K.E., Winnett, A.V., Chari, T. & Pachter, L. (2024) ‘Efficient and accurate detection of viral sequences at single-cell resolution reveals putative novel viruses perturbing host gene expression’, bioRxiv: The Preprint Server for Biology, p. 2023.12.11.571168.

<sup>[2](https://www.nature.com/articles/s41587-021-00870-2)</sup> Melsted, P., Booeshaghi, A.S., Liu, L., Gao, F., Lu, L., Min, K.H. (Joseph), da Veiga Beltrame, E., Hjörleifsson, K.E., Gehring, J. & Pachter, L. (2021) ‘Modular, efficient and constant-memory single-cell RNA-seq preprocessing’, _Nature Biotechnology_, 39(7), pp. 813–818.

<sup>[3](https://www.nature.com/articles/nbt.3519)</sup> Bray, N.L., Pimentel, H., Melsted, P. & Pachter, L. (2016) ‘Near-optimal probabilistic RNA-seq quantification’, Nature Biotechnology, 34(5), pp. 525–527.
