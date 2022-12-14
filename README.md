
## Combo-seq

Combo-seq is an analysis program which can process raw FASTA/FASTQ sequencing reads produced using the NEXTFLEX® Combo-Seq™ mRNA/miRNA Kit. It performs (1) pair end read trimming, (2) quality analysis, (3) genome indexing and alignment, (4) miRDP2 miRNA identification, (5) separates sRNA from mRNA 
and (6) creates a gene count matrix which can be used for further differential expression analysis. Usage of this pipeline on input reads prepared using the Combo-Seq kit presents a workflow for sRNA/mRNA combined analysis which replaces the tradtional need to separetly analyse mRNA and sRNA samples, reducing cost and time input. 


## Getting Started

``` 
$ git clone https://github.com/jadedavis5/combo-seq
$ cd combo-seq/
$ nano workflow.nf

# change parameters to desired input 
$ params.genome_file = "/path/to/genome.fa"
$ params.gtf_file = "/path/to/genome.gtf"
$ params.outdir = "/path/to/outdirectory"
$ params.reads = '/path/to/raw/reads/fa'
$ params.miRDP2 = "/path/to/miRDP2"

#run using:
$ nextflow run workflow.nf

#if the run is interrupted, resume from last point using:
$ nextflow run workflow.nf -resume
```

### Dependencies
The minimum required dependancies to run Combo-Seq are:
#### idk if i am supposed to list them out like this or use a package manager
* Bash >= 3.2 or Java >= 11
* [Nextflow](https://github.com/nextflow-io/nextflow) >= 22.04.3
* 



## Contact
- Email address
- Google Group/mailing list (if applicable)
- IRC or Slack (if applicable)

## License
### idk what liscence to use https://choosealicense.com/
