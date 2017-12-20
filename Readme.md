Description
===========

Nimbus is a software suite for the analysis of amplicon based sequencing data. Nimbus includes tools for data preprocessing, alignment, variant calling, quality control and visualization. The source amplicons are tracked throughout alignment and variant calling allowing probable false positive variants present in a single amplicon to be distinguished from real variants present in all amplicons at the locus. Nimbus also determines the number of reads mapped to each amplicon.

A scientific publication describing Nimbus is currently in preparation.

Build instructions
------------------

The following libraries should be present and installed on your system.

* [BOOST libraries](http://www.boost.org/)
  * BOOST development libraries
* [HTSlib (version >= 1.1.; part of samtools)](http://www.htslib.org/) 
* [SAMTOOLS](http://www.htslib.org/)
* [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
* [Python 2.7](https://www.continuum.io/downloads)
* [pysam](http://pysam.readthedocs.io/en/latest/api.html)
* [Make](https://www.gnu.org/software/make/)

In most situations, the Make and the Boost libraries are easier to install via the 
package managers that are present on most linux systems (`yum` for Redhat
or `apt` for Ubuntu). For Python, we recommend using the Anaconda, but the system
wide Python may also work. 

To compile Nimbus, go to the directory with the freshly cloned git 
repository and export the locations of the HTS-lib include and 
library files. Next use `make` to compile.

```bash
cd ${nimbus_git_repository}
export HTS_INCLUDE=${path_to_samtools}/htslib-${version}/htslib/
export HTS_LIB=${path_to_samtools}/htslib-${version}/libhts.a
make
```

After a successfull build 3 tools should be present in a newly created `${nimbus_git_repository}/bin/` directory. These tools are `nimbus_align`, `nimbus_call` and `nimbus_trim`. The other tools of Nimbus are in the scripts directory.

Run instructions
----------------

### First time setup

To run Nimbus, the tools in the Nimbus repository should be successfully compiled (see below) and SAMtools and Annovar should be present. Additionally a genome reference and a amplicon design in BED format should be available. The genome reference should be indexed with the `samtools faidx` command.

In the default workflow, Annovar annotates the results with the refGene, genomicSuperDups, avsnp147, and gerp++gt2 databases. These should all be present at the Annovar database location. Other annotation databases can be added by downloading them via Annovar and adding them in the workflow at the `protocol` and `operation` variables.

### How to run the workflow

To run the Nimbus workflow, copy the FastQ files of the samples you wish to process in an directory. These FastQ files should be preferably named `${samplename}_R1.fastq` for read 1 and `${samplename}_R2.fastq`. These files can be compressed with `gzip`. The workflow can also be started with bam files that were previously aligned with `nimbus_align`.

After preparing the run folder, set the following export variables and run the workflow with `make`.

```bash
export PATH_NIMBUS=${path_to_nimbus}/
export PATH_SAMTOOLS=${path_to_samtools}/
export PATH_PYTHON=${path_to_python}/bin/
export PATH_ZCAT=/bin/
export PATH_GZIP=/bin/
export PATH_ANNOVAR=${path_to_annovar}/
export ANNOVARDB=${path_to_annovar_database}/
export AMPLICON_DESIGN=${path_to_amplicon_design}/amplicon_design.bed
export GENOME_REFERENCE=${path_to_genome_refernece}/genome_reference.fa

make -j 3 -f ${PATH_NIMBUS}/workflows/workflow.mak
```

After a successfull run, the following output files will have been created.

| File                        | Description |
|:----------------------------|:------------|
| ${samplename}.srt.bam       | A BAM file with all the alignments |
| ${samplename}.discarded.bam | A BAM file with the discarded alignments |
| ${samplename}.passed.bam    | A BAM file with the passed alignments |
| ${samplename}.flagstat.txt  | A text file with the alignment statistics |
| ${samplename}.srt.blck      | A blck file with the counts per amplicon for all the alignments |
| ${samplename}.discarded.blck| A blck file with the counts per amplicon for the discarded alignments  |
| ${samplename}.passed.blck   | A blck file with the counts per amplicon for the passed alignments |
| combined.var.gz             | A var file with the sequence content at positions with an alternate allele |
| combined.txt                | An Annovar input file with the variants for all samples combined |
| combined.hg19_multianno.txt | An Annovar output file with the variants for all samples combined |
| combined.hg19.anno.txt      | A tabular file with the annovar annotations and variants, but with the header corrected |
| combined.hg19.mut.txt       | The variants in the [mut](https://software.broadinstitute.org/software/igv/MutationData) format for [IGV](http://software.broadinstitute.org/software/igv/home) |

Custom file formats
-------------------

### The var format

The var format is produced by `nimbus_call` and is used to record alternate alleles from a BAM file. Each record in a var file consists out of 1 line to denote the position and several subsequent lines with the sequenced bases. Comments can be placed on the top of the file with two hash signs (`##`).

The position line is tab-delimited and holds the following fields: chromosome, position, reference base number of subsequent lines, total read depth, total quality.

```var
chr1	13378	A	3	5	127
	A	1	13	sample, forward, chr1:13310-13468(-)
	A	3	112	sample, forward, chr1:13324-13511(+)
	C	1	2	sample, reverse, chr1:13310-13468(-)
```

### The blck format

The blck format produced by `nimbus_count.py` is a simple tab-delimeted format with the number of reads per amplicon. In column 1, each amplicon in the design is present encoded as follows: `chromosome:start-end(strand)`. The strand can either be a `+` or `-`. The second column has the number of reads for that amplicon.

```blck
chr10:100003494-100003993(+)	8
chr10:100003776-100003936(+)	52
chr10:100003821-100003905(+)	49
chr10:100003826-100004019(+)	16
chr10:100008136-100008288(-)	0
chr10:100008564-100008746(-)	50
```

Guides
------

In the Nimbus manuscript, we describe an analysis of HaloPlex exome samples. For the HaloPlex exome, we expanded the design to include potential off-target amplicons. The methods for this expansion are shown [here](haloplex_design_expansion/Readme.md).

The analysis of a set of uveal melanoma samples from a study by [Koopmans *et al*](https://www.nature.com/articles/modpathol201443) is available [here](guides/custom_haloplex_design/Readme.md).

Copy-number variant calling with Nimbus block files and [ExomeDepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html) is described [here](guides/CNV_calling_with_exomedepth/Readme.md).

A guide on how to call structural variants with Nimbus is available [here](guides/structural_variant_calling/Readme.md).
