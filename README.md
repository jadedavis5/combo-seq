## Combo-seq

Combo-seq is an analysis program which can process raw FASTA/FASTQ sequencing reads produced using the NEXTFLEX® Combo-Seq™ mRNA/miRNA Kit. It performs:

1. pair end read trimming
2. quality analysis
3. genome indexing and alignment 
4. miRDP2 miRNA identification
5. separates sRNA from mRNA 
6. creates a gene count matrix which can be used for further differential expression analysis. 


Usage of this pipeline on input reads prepared using the Combo-Seq kit presents a workflow for sRNA/mRNA combined analysis which replaces the tradtional need to separetly analyse mRNA and sRNA samples, reducing cost and time input. 

## Running the pipeline using `run.sh` (EXPERIMENTAL)
The `run.sh` script is preferred to run the pipeline, and has additional dependencies separate to the pipeline:
1. python3
2. python3-toml

The script loads "settings.toml", and uses two Python 3 scripts in the "bin" folder:
1. `bin/generate_configs.py`
2. `bin/toml_parser.py`

The `run.sh` script is functional without Python 3 or its TOML module, but cannot autogenerate the Nextflow configs, most likely preventing the pipeline from functioning. In case that it is unable to, the pipeline automatically loads the pre-generated configs under `conf/generated`, tailored for Pawsey's supercomputer, Setonix.


## Dependencies
- libarchive-tools (specifically `bsdtar`)


## Getting Started


## Pipeline architecture
An automatic configuration generation system, using TOML templates to configure the pipeline, is available over the built-in Nextflow configurations. Additionally, the pipeline makes use of DSL2 modules, stored under `./modules`.


## Contact
Jade Davis 20558259@student.curitn.edu.au

## License
### idk what liscence to use https://choosealicense.com/

