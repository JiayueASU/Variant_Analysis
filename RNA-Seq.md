# RNA-seq
<br />

## Why we have to study RNA-seq?

Suppose we have two bunches of neural cells. A bunch of normal cells, and a bunch of mutated cells. **We want to know what genetic mechanism is causing the difference, which means we want to look at the differences in gene expression.**

As we know there is a group of chromosome in each neural cells. Each chromosome carries some genes. Among these genes, some of them are active, some are not. The active one has mRNA transcripts, while the inactive one does not have transcripts. **High throughput sequencing tells us which genes are active, and how much they are transcribed. We can use RNA-seq to measure gene expression in normal cells, and then use it to measure gene expression in mutated cells. After that, we can compare the two cell types and figure out what is different in the mutated cells.**

---

<br />

## Three steps for RNA-seq

There are three steps for RNA-seq: **First, prepare a sequencing library; Second, sequence; Third, data analysis.**

### Prepare a sequencing library

To prepare an RNA-seq library, we first isolate the RNA, then break the RNA into small fragments. The reason we do this is because RNA transcripts can be thousands of bases long, but the sequencing machine can only sequence short fragments. 

----

After that, we convert the RNA fragments into double stranded DNA. Double stranded DNA is more stable than RNA and can be easily amplified and modified. 

---

Then, we add sequencing adaptors. The adapters do two things: 

1) Allow the sequencing machine to recognize the fragments;

2) Allow you to sequence different samples at the same time, since different samples can use different adaptors.

Notice that adding sequencing adaptors does not work 100% of the time. Some fragments may not have adaptors.

---

After we add sequencing adaptors, we do PCR amplify. Only the fragments with sequencing adaptors are amplified.

---

The final step is do quality control. This step focuses on two aspects. Verify library concentration and library fragment lengths.

---

### Sequence

The fragment of DNA is vertical, because that is how it is inside the sequencing machine. Each sequencing machine has a lot of grids. We call each grid a "flow cell". Each flow cell contains eight lanes. Each lane contains multiple tiles. Each tile is imaged four times per cycle. The machine has fluorescent probes that are color coded according to the type of nucleotide they can bind to. The probes are attached to the first base in each sequence. Once the probes have attached, the machine takes a picture of the flow cell. Then the machine washes the color off of the probes. Then probes are bound to the next base in each fragment. The process repeats until the machine has determined each sequence of nucleotides.

---

Sometimes a probe will not shine as bright as it should and the machine is not super confident that it is calling the correct color. Quality scores, that are part of the output, reflect how confident the machine is that it correctly called a base. Another reason you might get a low quality score is when there are lots of probes the same color in the same region. This is called "low density", and the over abundance of a single color can make it hard to identify the individual sequences. "Low density" is especially a problem when the first few nucleotides that is when the machine determines where the DNA fragments are located on the grid.

---

The raw data after sequencing is composed of four lines of data. The first line is always start with "@" is a unique ID for the sequence that follows. The second line contains the bases called for the sequenced fragment. The third line is always a "+" character, and this line is always empty, do not know why. The fourth line contains the quality scores for each base in the sequenced fragment. A typical sequence run with 400,000,000 reads will generate a file containing 1.6 billion lines of data.

---

#### Preprocessing

We need filter out garbage reads, align the high quality reads to a genome, and count the number of reads per gene.

##### Garbage reads

Garbage reads are reads with low quality base calls and are clearly artifacts of the chemistry.

##### Align the high quality reads to a genome

We first split the genome into small fragments. Then, we index of all the fragments and locations. After that we split the sequence reads into fragments. Match the read fragments to the genome fragments. The genome fragments that matched the read fragments will determine a location (chromosome and position) in the genome. The reason why we break the sequences up into small fragments is because it allows us to align reads even if they are not exact matches to the reference genome.

##### Count the number of reads per gene

Once we know the chromosome and position for a read, we can see if it falls within the coordinates of a gene. After we count the reads per gene, we end up with a matrix of numbers. The first column of the matrix is the names of genes. The human genome has about 20,000 genes, so this matrix has about 20,000 rows. The remaining columns contain counts for each sample you sequenced.

---

The last thing we do before analysis is normalize the data. This is because each sample will have different number of reads assigned to it, due to the fact that one simple might have more low quality reads, or another sample might have a slightly higher concentration on the flow cell.

---

### Data Analysis

We use Principle Component Analysis (PCA) or something like it to plot the data.