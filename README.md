# RNA-seq
## Why we have to study RNA-seq?

Suppose we have two bunches of neural cells. A bunch of normal cells, and a bunch of mutated cells. **We want to know what genetic mechanism is causing the difference, which means we want to look at the differences in gene expression.**

As we know there is a group of chromosome in each neural cells. Each chromosome carries some genes. Among these genes, some of them are active, some are not. The active one has mRNA transcripts, while the inactive one does not have transcripts. **High throughput sequencing tells us which genes are active, and how much they are transcribed. We can use RNA-seq to measure gene expression in normal cells, and then use it to measure gene expression in mutated cells. After that, we can compare the two cell types and figure out what is different in the mutated cells.**

## Three steps for RNA-seq

There are three steps for RNA-seq: **First, prepare a sequencing library; Second, sequence; Third, data analysis.**

### Prepare a sequencing library

To prepare an RNA-seq library, we first isolate the RNA, then break the RNA into small fragments. The reason we do this is because RNA transcripts can be thousands of bases long, but the sequencing machine can only sequence short fragments. 

After that, we convert the RNA fragments into double stranded DNA. Double stranded DNA is more stable than RNA and can be easily amplified and modified. 

Then, we add sequencing adaptors. The adapters do two things: 

<ol>
  <li>Allow the sequencing machine to recognize the fragments;</li>
  <li>Allow you to sequence different samples at the same time, since different samples can use different adaptors.</li>
</ol>

Notice that adding sequencing adaptors does not work 100% of the time. Some fragments may not have adaptors.

After we add sequencing adaptors, we do PCR amplify. Only the fragments with sequencing adaptors are amplified.

The final step is do quality control. This step focuses on two aspects. Verify library concentration and library fragment lengths.

### Sequence

The fragment of DNA is vertical, because that is how it is inside the sequencing machine. Each sequencing machine has a lot of grids. We call each grid a "flow cell". Each flow cell contains eight lanes. Each lane contains multiple tiles. Each tile is imaged four times per cycle. The machine has fluorescent probes that are color coded according to the type of nucleotide they can bind to. The probes are attached to the first base in each sequence. Once the probes have attached, the machine takes a picture of the flow cell. Then the machine washes the color off of the probes. Then probes are bound to the next base in each fragment. The process repeats until the machine has determined each sequence of nucleotides.

Sometimes a probe will not shine as bright as it should and the machine is not super confident that it is calling the correct color. Quality scores, that are part of the output, reflect how confident the machine is that it correctly called a base. Another reason you might get a low quality score is when there are lots of probes the same color in the same region. This is called "low density", and the over abundance of a single color can make it hard to identify the individual sequences. "Low density" is especially a problem when the first few nucleotides that is when the machine determines where the DNA fragments are located on the grid.

The raw data after sequencing is composed of four lines of data. The first line is always start with "@" is a unique ID for the sequence that follows. The second line contains the bases called for the sequenced fragment. The third line is always a "+" character, and this line is always empty, do not know why. The fourth line contains the quality scores for each base in the sequenced fragment. A typical sequence run with 400,000,000 reads will generate a file containing 1.6 billion lines of data.

### Preprocessing

We need filter out garbage reads, align the high quality reads to a genome, and count the number of reads per gene.

### Garbage reads

Garbage reads are reads with low quality base calls and are clearly artifacts of the chemistry.

### Align the high quality reads to a genome

We first split the genome into small fragments. Then, we index of all the fragments and locations. After that we split the sequence reads into fragments. Match the read fragments to the genome fragments. The genome fragments that matched the read fragments will determine a location (chromosome and position) in the genome. The reason why we break the sequences up into small fragments is because it allows us to align reads even if they are not exact matches to the reference genome.

### Count the number of reads per gene

Once we know the chromosome and position for a read, we can see if it falls within the coordinates of a gene. After we count the reads per gene, we end up with a matrix of numbers. The first column of the matrix is the names of genes. The human genome has about 20,000 genes, so this matrix has about 20,000 rows. The remaining columns contain counts for each sample you sequenced.

The last thing we do before analysis is normalize the data. This is because each sample will have different number of reads assigned to it, due to the fact that one simple might have more low quality reads, or another sample might have a slightly higher concentration on the flow cell.

### Data Analysis

We use Principle Component Analysis (PCA) or something like it to plot the data.

<br />

# Variant Analysis
## Quality Assessment

### FastQC

Download Page: <http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc>

Unzip FastQC: `unzip fastqc_v0.11.9.zip`

Set environmental variables: `echo export PATH=$PATH:/data/notebook/Jerry/Tools/FastQC/fastqc >> ~/.bashrc`

`source ~/.bashrc`

Usage of FastQC: `fastqc -o ./ -t 6 V300035135_L03_531_1.clean.fq.gz ...`

Or you can process multiple .fq.gz files in batch with a .sh script:

```javascript
for id in *fastq
do
echo $id
/data/notebook/Jerry/Test/Data/Data20200323 $id
Done
```

Explanation: 
-o --outdir     

> Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it. If this option is not set then the output file for each sequence file is created in the same directory as the sequence file which was processed.

-t --threads    

> Specifies the number of files which can be processed simultaneously. Each thread will be allocated 250MB of memory so you shouldn't run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine.

After the execution, an .html file will be generated. This file is the fastqc report of the associate fq data. The explanation of this report can be found at: 
<http://www.bio-info-trainee.com/95.html>

<img src="https://github.com/JiayueASU/Variant_Analysis/blob/main/1FastQC.png" width=50% height=50%>

### Trimmomatic

Download Page: <http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip>

Unzip Trimmomatic: `unzip Trimmomatic-0.39.zip`

The tutorial of using Trimmomatic can be found at: 
<http://www.bio-info-trainee.com/1958.html>

<http://www.usadellab.org/cms/?page=trimmomatic>

Installation of Java: `yum install [java-1.8.0-openjdk.x86_64](java-1.8.0-openjdk.x86_64) `

For more details about using Trimmomatic: `java -jar Trimmomatic-0.39.jar -h`

Usage of Trimmomatic (Paired End): `java -jar trimmomatic-0.39.jar PE -phred33 /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_1.clean.fq.gz /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_2.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_1_paired.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_1_unpaired.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_2_paired.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_2_unpaired.clean.fq.gz /data/notebook/Jerry/Tools/Trimmomatic-0.39/adapters/ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

Output_paired: Usable data

Output_unpaired: Removed adapters, leading low quality, and trailing low quality

Usage of Trimmomatic (Single End): `java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

After the execution, we could obtain paired data for next step.

<img src="https://github.com/JiayueASU/Variant_Analysis/blob/main/2Trimmomatic.png" width=50% height=50%>

## Read Alignment

### BWA-mem

Download Page:
<https://sourceforge.net/projects/bio-bwa/files/>

Tutorial: <https://www.jianshu.com/p/1552cc6ac3be>

Unzip bwa-0.7.17.tar.bz2: `tar -jxvf bwa-0.7.17.tar.bz2`

Install bwa:`cd bwa-0.7.17` `make`

After `make`, execute bwa file: `./bwa`

Set environmental variables: `export PATH=$PATH:/data/notebook/Jerry/Tools/bwa-0.7.17/`

`echo export PATH=$PATH:/data/notebook/Jerry/Tools/bwa-0.7.17 >> ~/.bashrc` `source ~/.bashrc`

Create a new folder for mm10.fasta: `cd /data/notebook/Jerry/Test/Reference`

Create a new folder at Output: `mkdir bwa_index`

`mkdir mm10` `cd mm10`

Download mm10.fast: `wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz`

cd to: `/data/notebook/Jerry/Test/Reference/mm10/Mus_musculus/UCSC/mm10/Sequence/BWAIndex`

Copy genome.fa to the Output folder: `cp genome.fa /data/notebook/Jerry/Test/Output/bwa_index/`

Generate index sequence: `bwa index genome.fa`

<img src="https://github.com/JiayueASU/Variant_Analysis/blob/main/3BWAindex.png" width=50% height=50%>

After 613 iterations, five new files are generated: genome.fa.amb, genome.fa.ann, genome.fa.bwt, genome.fa.pac, and genome.fa.sa.

Use BWA-mem to obtain .sam file: `bwa mem genome.fa /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_1.clean.fq.gz /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_2.clean.fq.gz > aln-pe.sam`

<img src="https://github.com/JiayueASU/Variant_Analysis/blob/main/4BWAmem.png" width=50% height=50%>

After 11142.574 sec, aln-pe.sam file is generated.

<img src="https://github.com/JiayueASU/Variant_Analysis/blob/main/5AfterBWA.png" width=50% height=50%>

## Variant Identification

### Samtools

First download samtools: `wget -c https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2`

Then unzip samtools-1.9: `tar jxvf samtools-1.9.tar.bz2`

Set Configure: `./configure --prefix=/data/notebook/Jerry/Tools/samtools-1.9`

If it shows a bug like: "fatal error: curses.h: No such file or directory", please install libncurses5-dev (ubuntu) or curse-devel (centos).

Then `make` `make install`

If it shows a bug on "htslib-1.9", install this package.

Generate a new folder at Output: `mkdir samtools_bam`

Copy the generated .sam file to samtools_bam: `cp aln-pe.sam ~/samtools_bam`

Generate .bam file (15 min): `samtools view -bS aln-pe.sam > aln-pe.bam`

Sort the generated .bam file (25 min): `samtools sort -n aln-pe.bam -o aln-pe.sort.bam`

<img src="https://github.com/JiayueASU/Variant_Analysis/blob/main/8Samtools.png" width=50% height=50%>

### GATK

Tutorial: <http://www.bio-info-trainee.com/3144.html>

First download GATK: `cd /data/notebook/Jerry/Tools`

`wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip`

Then unzip gatk-4.0.2.1.zip: `unzip gatk-4.0.2.1`

Set environmental variables: `echo export PATH=$PATH:/data/notebook/Jerry/Tools/gatk-[4.0.2.1](4.0.2.1) >> ~/.bashrc ` `source ~/.bashrc`

Use GATK to mark duplicates (18 min): `java -jar gatk-package-4.0.2.1-local.jar MarkDuplicates \-I /data/notebook/Jerry/Test/Output/samtools_bam/aln-pe.sort.bam -O aln-pe.sort.markdup.bam -M aln-pe.sort.markdup.bam.metrics`

Sort the generated .bam file again: `samtools sort aln-pe.sort.markdup.bam -o aln-pe.sort1.markdup.bam`

Create index for sorted marked .bam file to generate aln-pe.sort1.markdup.bam.bai: `samtools index aln-pe.sort1.markdup.bam`

Genreate a new folder at Reference, named GATK. Then download .vcf files from: `wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz`

Switch to /data/notebook/Jerry/Test/Output/bwa_index, generate .fai file: `samtools faidx genome.fa`

Generate .dict file based on genome.fa: `samtools dict genome.fa > genome.dict`

Do recalibration at the folder of gate: `java -jar gatk-package-4.0.2.1-local.jar BaseRecalibrator -R /data/notebook/Jerry/Test/Output/bwa_index/chromosomes.fa -I /data/notebook/Jerry/Test/Output/samtools_bam/aln-pe-0406-sort1-markdup.bam -O /data/notebook/Jerry/Test/Output/samtools_bam/aln-pe-0406-sort1-markdup.recal.table --known-sites /data/notebook/Jerry/Test/Reference/mm10SNP/mgp.v3.snps.rsIDdbSNPv137.vcf --known-sites /data/notebook/Jerry/Test/Reference/mm10SNP/mgp.v3.indels.rsIDdbSNPv137.vcf`

Use IndexFeatureFile to generate vcd.idx file (This step solves the problem >= 1 but = 0): `java -jar gatk-package-4.0.2.1-local.jar IndexFeatureFile -F /data/notebook/Jerry/Test/Reference/GATK/00-All.vcf `

Check out this tutorial: <https://www.jianshu.com/p/fe0c876563b0>

Check out the detail on .bam file: `samtools view -H aln-pe-0406-sort1-markdup.bam`

Merge all the .fa file at Chromosomes to one .fa file: `cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chrM.fa chrX.fa chrY.fa > chromosomes.fa`
