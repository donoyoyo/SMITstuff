# SMITstuff


1.)
Starting from the Cleveland dataset (the raw fastq files),
preliminary QC was done with fastqc on the raw fastq files.

link  here

http://intron.ucsc.edu/jen/jan2019/Cleveland/

This contains the raw fastq which has been deposited.
https://www.ncbi.nlm.nih.gov/sra/PRJNA975902

2.)
bbmap clumpify was used to remove duplicates and clump the reads.
That script is here (essentially default parameters)

http://intron.ucsc.edu/jen/april_23/clump

bbmap version 37.90

3.)
Cutadapt v1.11 was used to remove adapters and filter reads based on
presence of adapter, error rate, minimum length,
minimum overlap between read and adapter.
The script for that is here

http://intron.ucsc.edu/jen/april_23/cutadapt_smSMIT2.pl

with parameters

-O 8 -n 2 -m 23 -e 0.11

--discard-untrimmed


4.)
Alignment was done using  hisat2-2.1.0
with parameters --no-mixed --max-intronlen 10000 and aligned to sacCer3

Alignment statistics are here

http://intron.ucsc.edu/jen/april_23/alignment_stats

The script for that is here

http://intron.ucsc.edu/jen/april_23/raphisat2

The aligned reads are here

http://intron.ucsc.edu/jen/april_23/alignedreads/

5.)
The intron annotations were untouched from Tara and are contained here

http://intron.ucsc.edu/jen/april_23/SMIT_accessory_files.zip

Those annotation files anotate gene,intron and terminal exon length.

The R script to generate the splicing files is here

http://intron.ucsc.edu/jen/april_23/2SMIT_processing_local.R

as described in Tara's paper. This produced the files with the
splicing values as well as these plots

http://intron.ucsc.edu/jen/jan2019/processed_data100/smit_curves/panels/

As well as the SMIT_accessory_files, the gene list of 62 relevent genes that were
in the study is input.
that is here

http://intron.ucsc.edu/jen/april_23/genelist62

The input is the sacCer3 gene annotations , mapped read bed files , primed genes
and intron less genes as controls.
Only reads involved with the primed genes and intron less genes as controls are used.
Spliced reads and unspliced reads are counted by position relative to the 3'ss.
The raw splice value is (spliced count)/(unspliced count + spliced count).
The distribution and probability of insert length is determined from the data.
The raw splice value is then normalized using the position, insert length and lengths of the gene products (spliced and unspliced) as well as the insert length probabiltiy function.

The script is in this one SMIT_accessory_files/smitData.R

In particular, 

SMIT_accessory_files/smitData.R includes R data shaping, mapping, and positional R scripts 

SMIT_accessory_files/commonFunctions.R includes mostly statistical R functions

6.)
The splicing numbers were binned by 20 nt,  plotted and used to
calculate a wilcoxon paired  p-value to determine
if the difference was significant between conditions.
This was also used to create a score using the difference
between conditions for each gene for the region from start of signal to 200 nt after

wilcoxon two sample test is also called Mann-Whitney

That is the last step with the recent plots

http://intron.ucsc.edu/jen/april_23/finorm2200/

This being the script

http://intron.ucsc.edu/jen/april_23/finorm2200/normedaauc200.pl

