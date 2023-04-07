#!/usr/bin/perl -w 



#out1=1A_S50_R1.fastq  out2=1A_S50_R2.fastq 
#out1=1B_S51_R1.fastq  out2=1B_S51_R2.fastq 
#out1=1C_S52_R1.fastq  out2=1C_S52_R2.fastq 
#out1=1D_S53_R1.fastq  out2=1D_S53_R2.fastq 

@files = qw(1A_S50
1B_S51
1C_S52
1D_S53
);

for($i=0;$i<4;$i++){
$name = $files[$i];
print "R1 .. proceed with cutadapt \n";

#system("cutadapt -g CATTGATGGTGCCTACAG -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A CTGTAGGCACCATCAATG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -n 2 -m 23 -e 0.11  -o $name\_R1_filt1.fastq -p $name\_R2_filt1.fastq ../secondpass_clump_dupremoval/$name\_R1.fastq  ../secondpass_clump_dupremoval/$name\_R2.fastq");
#system("cutadapt --discard-untrimmed -g CATTGATGGTGCCTACAG -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A CTGTAGGCACCATCAATG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -n 2 -m 23 -e 0.11  -o $name\_R1_filt1.fastq -p $name\_R2_filt1.fastq ../secondpass_clump_dupremoval/$name\_R1.fastq  ../secondpass_clump_dupremoval/$name\_R2.fastq");
system("cutadapt --discard-untrimmed -g CATTGATGGTGCCTACAG -O 8 -n 2 -m 23 -e 0.11  -o $name\_R1_filt2.fastq -p $name\_R2_filt2.fastq ../secondpass_clump_dupremoval/$name\_R1.fastq  ../secondpass_clump_dupremoval/$name\_R2.fastq");


}




