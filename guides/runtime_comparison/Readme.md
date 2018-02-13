Runtime comparison
=================

We compared the runtimes of the Nimbus alignment and variant calling to those of other commonly used tools. For the alignment, we compared to the [`bwa mem` aligner](https://arxiv.org/abs/1303.3997) which is an fast aligner suitable for short and long reads. The variant calling was compared that of [FreeBayes](https://arxiv.org/abs/1207.3907), the [GATK Haplotypecaller (version 3.7)](https://www.biorxiv.org/content/early/2017/11/14/201178.1), and [Vardict](https://www.ncbi.nlm.nih.gov/pubmed/27060149).

The run times of the tools were determined using the Linux `time` tool. The memory usage was monitored using the `top` process monitoring tool. Eventhough `top` does not provide overal memory statistics, it is still valid to use in this case, as all tools had a stable memory usage on our system. The tools were run on a DELL Precision tower 7910 running Ubuntu Linux 16.04.2 LTS with 2 Intel(R) Xeon(R) CPU E5-2640 processors (32 cores total), 128 GB of RAM memory and 10 TB of internal storage organised in a RAID5 array. During runtime testing no other demanding processes were performed on the machine.

The performance of `nimbus align` was compared to that of `bwa mem` with 2 datasets: a HaloPlex Exome dataset with ~40 million reads and 2.8 million amplicons expected and a HaloPlex custom panel with ~940 thousand reads over 4 thousand amplicons (Table 1).

Table 1. Test alignment cases

| Dataset               | Sample  | Reads    | Amplicons |
| :-------------------- | :-----: | :------: | :-------: |
| HaloPlex Exome        | NA15510 | 41417946 | 2803399   |
| HaloPlex custom panel | -       | 937468   | 4334      |

The alignment performed by Nimbus is substantially slower than `bwa mem` for the exome dataset (Table 2). For the HaloPlex custom panels, `bwa mem` is still faster than `nimbus align` (Table 2). However, in this case the difference is less than for the HaloPlex exome dataset. The Nimbus alignment slows down with the number of amplicons as each read-pair can be matched to a greater number of amplicons.

During the alignment, Nimbus requires approximately 50 % of the memory of `bwa mem`. The Nimbus alignment only loads the sequence underlying the amplicons allowing more `nimbus align` processes to be performed in parallel on a system.

In practice, the `nimbus align` alignment still takes place in a reasonable time (~1 day) after which the reads are annotated with their source amplicons. This functionality is vital for tracking amplicons during variant calling and amplicon quantification.

Table 2. Alignment runtimes

| Aligner      | Dataset               | clock time (m) | CPU time (m) | RAM used (GB) |
| :----------- | :-------------------: | :------------: | :----------: | :-----------: |
| Nimbus align | HaloExome             | 1387.9         | 5758.9       | ~ 3           |
| Nimbus align | HaloPlex custom panel | 3.2            | 6.5          | < 2           |
| bwa mem      | HaloExome             | 22.4           | 88.8         | ~ 6           |
| bwa mem      | HaloPlex custom panel | 1.0            | 3.8          | ~ 6           |

Variants were called in the HaloPlex exome dataset aligned with `nimbus align` (Table 3). The Nimbus variant calling was the fastest of the 4 tools tested here. Nimbus was ~ 10 minutes faster than the runner-up Vardict* and 1 hour faster than FreeBayes. The GATK Haplotypecaller is the slowest from the tools tested here and required almost 3 hours for variant calling.

Table 3. Variant calling runtimes

| Variant Caller       | preprocessing time (m) | calling time (m) | postprocessing time (m) | total time (m) |
| :------------------- | :--------------------: | :--------------: | :---------------------: | :------------: |
| Nimbus call          | 0.00                   | 62.66            | 15.75                   | 78.41          |
| FreeBayes            | 0.00                   | 138.61           | 0.00                    | 138.61         |
| GATK HaploTypeCaller | 110.53                 | 44.45            | 0.00                    | 154.98         |
| Vardict              | 0.00                   | 85.98            | 2.00                    | 87.99          |

\* Vardict seemed to be unsuitable for HaloPlex datasets as up to 50% of the variants called by other variant caller are not found in the Vardict output.