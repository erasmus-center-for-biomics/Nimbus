HaloPlex design expansion
=========================

Introduction
------------

In the HaloPlex DNA capture procedure, DNA fragments are captured by ligating partially double stranded probes to the ends of particular restriction fragments. The fragment specificity is achieved through the ligation reaction as this reaction has little tolerance for mismatches between the 2 fragments to be joined. By ligating a partially double stranded DNA molecule with 2 overhangs to specific restriction fragments are captured and circular DNA molecule is created (a "halo"). The halo is then amplified and paired-end reads are generated from it. More information on the HaloPlex procedure can be found on the [Agilent website](https://www.genomics.agilent.com/article.jsp?pageId=3081).

The HaloPlex DNA capture method is very specific, but not flawless. Off-target restriction fragments can be captured if their ends are similar to those in the design. To determine these fragments, we developed tools and procedures to determine possible off-target restriction fragments from the genome sequence and enzyme mixes. This procedure can also be employed to transfer a HaloPlex design from one genome build to another.

Prequisites
-----------

* Python version >= 2.6
* The HaloPlex amplicon design
* A FastA file with the genome sequence of the target organism
* A file with the restriction enzyme mix used in the HaloPlex reactions
* A GNU/Linux or UNIX server with cat and sort

Procedure
---------

To determine off-target DNA fragments, the genome is restricted *in silico* with a restriction enzyme mix and then the ends of amplicons are matched to the *in silico* fragments.

The enzyme mix needs to specified in a tab-delimited file with the name of the enzymes and the IUPAC code for the enzyme recognition site. The specific cut sites should be indicated with / for the leading strand and | for the lagging strand. The restriction enzymes used in an experiment should be obtained from Agilent.

Here follows the contents of an example enzyme file with HindIII, NlaIII, NlaIV and FokI:

```enzyme file
#name   site
HindIII A/AGCT|T
NlaIII  |CATG/
NlaIV   GGN/NCC
FokI    GGATGNNNNNNNNN/NNNN|
```

To generate a restriction map from a genome use the `map.py` script in the `bin` directory.

```bash
python bin/map.py --fasta genome.fasta --enzymes enzyme_mix_A.txt --out genome_restricted_with_enzyme_mix_A.bed
```

Next the ends of the fragments in the design are matched to those in the restriction map. To perform this matching, we need to specify how long the matching sequences must be and how many mismatches are allowed. Furthermore, we can specify the number of "safe bases" in which no mismatches are allowed. For the HaloPlex exome design expansion, the matching sequences were 11 bases long with 5 safe bases and only 1 mismatch allowed over both ends of the amplicon.

To determine the fragments with matching ends use the `homologs_2.py` script in the `bin` directory. This script can take quite some time to run.

```bash
python bin/homologs_2.py \
    --bed haloplex_amplicon_design.bed \
    --all genome_restricted_with_enzyme_mix_A.bed \
    --fasta genome.fasta \
    --ksize 11 \
    --mismatches 1 \
    --safebases 5 \
    --out expanded_haloplex_design_mix_A.bed
```

For a typical HaloPlex design multiple enzyme mixes are used. To combined these fragments in a single design use the following command.

```bash
cat expanded_haloplex_design_mix_A.bed expanded_haloplex_design_mix_B.bed | \
    sort -k1,1 -k2,2n > expanded_haloplex_design.bed
```

Trouble shooting
----------------

Q) I do not know the enzymes that were used in the HaloPlex design.

A) This information should be provided by Agilent. We are not allowed to disclose which enzymes are used in our designs.

Q) The design expansion procedure fails for a HaloPlex HS design.

A) This is correct. In the HaloPlex HS design the end of the amplicons should be expanded by 1 base on the leading strand and the start should be decreased on the lagging strand. The following Awk script performs this action.

```awk
//{
        if( $6 == "+" ) {
                print $1 "\t" $2 "\t" $3 + 1 "\t" $4 "\t" $5 "\t" $6 ;
        } else if( $6 == "-" ) {
                print $1 "\t" $2 - 1"\t" $3 "\t" $4 "\t" $5 "\t" $6 ;
        }
}
```

This modification of the design is also required for the Nimbus alignment.
