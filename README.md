
### Combo-seq

Combo-seq is an analysis program which can process raw FASTA/FASTQ sequencing reads produced using the NEXTFLEX® Combo-Seq™ mRNA/miRNA Kit. It performs (1) pair end read trimming, (2) quality analysis, (3) genome indexing and alignment, (4) miRDP2 miRNA identification, (5) separates sRNA from mRNA 
and (6) creates a gene count matrix which can be used for further differential expression analysis.  

its reproducible... saves your progress 
- State if it is out-of-the-box user-friendly, so it’s clear to the user.

- State its goals/what problem(s) it solves.
- Note and briefly describe any key concepts (technical, philosophical, or both) important to the user’s understanding.
- Link to any supplementary blog posts or project main pages.
- Note its development status.
- Include badges.
- If possible, include screenshots and demo videos.

### Core Technical Concepts/Inspiration

- Why does it exist?
- Frame your project for the potential user. 
- Compare/contrast your project with other, similar projects so the user knows how it is different from those projects.
- Highlight the technical concepts that your project demonstrates or supports. Keep it very brief.
- Keep it useful.

### Getting Started/Requirements/Prerequisites/Dependencies
# Dependencies 


``` 
$ git clone https://github.com/jadedavis5/combo-seq
$ cd combo-seq/
$ nano workflow.nf

# change parameters to desired input 
$ params.genome_file = "/path/to/genome.fa"
$ params.gtf_file = "/path/to/genome.gtf"
$ params.outdir = "/path/to/outdirectory"
$ params.reads = '/path/to/raw/reads/fa'
$ params.miRDP2 = "/path/to/miRDP2" # [Download here](https://sourceforge.net/projects/mirdp2/files/latest_version/)

#run using:
nextflow run workflow.nf

#if the run is interrupted, resume from last point using:
nextflow run workflow.nf -resume
```
[Download here](https://sourceforge.net/projects/mirdp2/files/latest_version/)
Include any essential instructions for:
- Getting it
- Installing It
- Configuring It
- Running it


### Contact
- Email address
- Google Group/mailing list (if applicable)
- IRC or Slack (if applicable)

### License
