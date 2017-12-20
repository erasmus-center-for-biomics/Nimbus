Structural variant calling with Nimbus
==========================

Introduction
------------

For structural variant (SV) detection, amplicons are designed over breakpoints. With some adaptations, Nimbus can detect and quantify these SVs. In this guide, the adaptions to the protocol are discussed.

This is a purely theoretical guide as we have not yet had experiments in which specific SVs have been targeted.

Amplicon design considerations
------------------------------

To effectively call SVs with amplicon based technologies, amplicons must be designed over the break-points. In the case of a typical PCR, primers should be placed on next to both sites of the break point (Figure 1). Furthermore, primers can also be designed to identify the absence of a structural variants.

One should note that the sequence left and right from the breakpoint should differ between the SV amplicon and the normal cases. Otherwise these cannot be distinguished based on sequence.

With a clear amplicon design, the presence and absence of specific amplicons show the presence or absence of structural variants. The analysis of such a design can be facilitated by Nimbus.

```amplicon
break-point = //

primer a                  -->
chromosome z: +++++++++++++++++++//++++++++++++++++++++++
primer a*:                            <--

primer b*                   -->
chromosome q: ===================//======================
primer b                              <--

Yields:

if SV:
          +++++++//======
if no SV:
          +++++++//++++++
          =====//======
if "heterozygous" SV:
          +++++++//======
          +++++++//++++++
          =====//======
```

Figure 1 amplicon with a break point

Adapted reference files
-----------------------

To analyse the amplicons described in figure 1, breakpoints should be added to the reference sequence. The `samtools faidx` command can be used to obtain the sequence around a break point in FastA format.

```breakpoint_fasta
echo ">breakpoint01" > to_add.fasta
samtools faidx chr(z):primerstart(a)-breakpoint(z) | tr -d '\n' >> to_add.fasta
samtools faidx chr(q):breakpoint(q)-primerend(b) >> to_add.fasta
```

Next this sequence can be appended to (a copy of the reference) the reference sequence.

```add_fasta
cat to_add.fasta >> hg19.fasta
```

The breakpoint amplicon should also be added to the design BED file. For regular PCR enrichment both the amplicon on the + and - strand should be added.

```add_amplicon
echo -e "breakpoint01\t0\tlength(amplicon)\tbreakpoint01\t0\t+" >> design.bed
echo -e "breakpoint01\t0\tlength(amplicon)\tbreakpoint01\t0\t-" >> design.bed
```

Nimbus analysis
---------------

With the modified reference files, the analysis with Nimbus is reasonable straightforward. The reads should be trimmed, aligned and amplicons counted. A quick overview of these commands are as follows.

```bash
# make a directory for the output per tool
mkdir -p logs

# trim read 1
zcat {samplename}_R1_001.fastq.gz | bin/nimbus_trim \
    -i - \
    --maximum-mismatches 2 --minimum-matches 2 --minimum-bases-remaining 40 \
    -a AGATCGGAAGAG -a CTGTCTCTTATA \
    -o ${samplename}_R1.tr.fastq 2>> logs/trim.errors.log >> logs/trim.messages.errors.log

# trim read 2
zcat {samplename}_R2_001.fastq.gz | bin/nimbus_trim \
    -i - \
    --maximum-mismatches 2 --minimum-matches 2 --minimum-bases-remaining 40 \
    -a AGATCGGAAGAG -a CTGTCTCTTATA \
    -o ${samplename}_R2.tr.fastq 2>> logs/trim.errors.log >> logs/trim.messages.errors.log

# alignment
bin/nimbus_align align \
    --forward ${samplename}_R1.tr.fastq \
    --reverse ${samplename}_R2.tr.fastq \
    --design design.bed \
    --fasta hg19.fasta \
    --key-size 7 \
    --workers 5 \
    --maximum-amplicons 1000 \
    --sam ${samplename}.sam 2>> logs/nimbus.errors.log >> logs/nimbus.messages.log

# convert to BAM
samtools addreplacerg \
    -r "ID:${sample}" \
    -r "CN:center_label" \
    -r "LB:library_label" \
    -r "SM:${samplename}" \
    -r "PL:platform_label" \
    -r "DS:${fullpath}" \
    ${samplename}.sam | samtools sort \
         --threads 4 \
        -T ${samplename}_tmp_sort \
        -o ${samplename}.srt.bam \
        - ;

# count reads per amplicon
python scripts/nimbus_count.py \
    --input ${samplename}.srt.bam \
    --design design.bed \
    --output ${samplename}.blck >> logs/count.messages.log 2>> logs/count.errors.log
```

The result of this analysis is a block file which includes the breakpoint amplicon on the virtual chromosome `breakpoint01` and the amplicons to show the *normal* situation.

Downstream processing
---------------------

To compare these situations, one can use [R](https://www.r-project.org/) to analyse the results. Microsoft excel can also be used if only a few amplicons are sequenced. To load the blck file and parse the chromosome, start, end, and strand information, the following code can be used:

```R
library(tidyverse)
library(stringr)

fn_block <- "samplename.block"

block <- read_tsv(fn_block, col_names=c("key", "depth"))
matcher <- str_match(block$key, "^(.+):([0-9]+)-([0-9]+)\\(([+,-])\\)")

# add chromosome, start, end, and strand
block <- block %>%
    mutate(
        chromosome = parse_character(matcher[,2]),
        start = parse_integer(matcher[,3]) + 1,
        end = parse_integer(matcher[,4]),
        strand = parse_character(matcher[,5])
    )

# for non-stranded protocols summarise depth per amplicon
sblock <- block %>%
    group_by(chromosome, start, end) %>%
    summarise(depth = sum(depth))
```

Further comparisons can also be performed in R and are dependent on the protocols and methods used.

In conclusion
-------------

This guide is meant to provide some guidance on how to analyse structural variants with Nimbus. With some adaptations to reference files, Nimbus can provide a direct readout of amplicons targeting chromosomal breakpoints. These readouts can be processed further in statistical software like R.
