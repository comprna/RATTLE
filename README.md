# RATTLE
Reference-free reconstruction and quantification of transcriptomes from long-read sequencing

* I. de la Rubia, A. Srivastava, W. Xue, J. A. Indi, S. Carbonell-Sala,  J. Lagarde,  M.M. Albà, E. Eyras. RATTLE: Reference-free reconstruction and quantification of transcriptomes from Nanopore long-read sequencing. bioRxiv doi: https://doi.org/10.1101/2020.02.08.939942


----------------------------
# Table of Contents
----------------------------

 * [Requirements](#requirements)
 * [Installation](#installation)
 * [Quick start](#quick-start)
 * [Running RATTLE](#running-rattle)
   * [Clustering](#clustering)
     * [Description of clustering parameters](#Description-of-clustering-parameters)
     * [Clustering summary](#clustering-summary)
     * [Cluster extraction](#cluster-extraction)
   * [Error correction](#error-correction)
     * [Description of error correction parameters](#Description-of-error-correction-parameters)  
   * [Polishing step](#polishing-step)
 * [Example datasets](#example-datasets)
   * [Human direct RNA sequencing](#human-direct-RNA-sequencing)
 * [Reference based benchmarking](#reference-based-benchmarking)
   * [ssCheck](#sscheck)
   * [Example dataset](#example-dataset)
   * [Cluster benchmark](#cluster-benchmark)
 * [Snakemake](#Snakemake)
# Requirements
GCC, G++ with C++14 suppport

# Installation
* Clone the repository
```
git clone --recurse-submodules https://github.com/comprna/RATTLE
```

* Build RATTLE
```
cd RATTLE
./build.sh
```
(this will generally take less than 1 minute)


# Quick start

We provide here some of the most common commands used when running RATTLE. 
**Note**: The commands and parameters are still under development and may be subject to changes in future versions

* Cluster cDNA Nanopore reads at gene level with 24 threads
```
$ ./rattle cluster -i reads.fq -t 24 -o . 
```

* Cluster cDNA Nanopore reads with labels at gene level with 24 threads
```
$ ./rattle cluster -i reads.fq -t 24 -o . 
```

* Cluster cDNA Nanopore reads at isoform level with 24 threads
```
$ ./rattle cluster -i reads.fq -t 24 --iso
```

* Cluster RNA Nanopore reads at isoform level with 24 threads
```
$ ./rattle cluster -i reads.fq -t 24 --iso --rna
```

* View clustering summary (csv with read_id,cluster_id)
```
$ ./rattle cluster_summary -i reads.fq -c clusters.out
```

* Extract 1 fastq file per cluster in clusters/ folder
```
$ mkdir clusters
$ ./rattle extract_clusters -i reads.fq -c transcripts.out -o clusters --fastq
```

* Correct reads with 24 threads using isoform clusters
```
$ ./rattle correct -i reads.fq -c clusters.out -t 24
```

* Polish RNA consensus sequences and build final transcriptome using 24 threads
```
$ ./rattle polish -i consensi.fq -t 24 --rna
```

# Running RATTLE

We provide here the details of each RATTLE command.

## Clustering

This is the first and most important step in RATTLE. This command will generate the first set of read clusters representing potential genes and transcripts.

```
$ ./rattle cluster -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file or compressed with .gz extension, will automatically check file extension (required)
    -l, --label
        labels for the files in order of entry (optional)
    -o, --output
        output folder (default: .)
    -t, --threads
        number of threads to use (default: 1)
    -k, --kmer-size
        k-mer size for gene clustering (default: 10, maximum: 16)
    -s, --score-threshold
        minimum score for two reads to be in the same gene cluster (default: 0.2)
    -v, --max-variance
        max allowed variance for two reads to be in the same gene cluster (default: 1000000)
    --iso
        perform clustering at the isoform level
    --iso-kmer-size
        k-mer size for isoform clustering (default: 11, maximum: 16)
    --iso-score-threshold
        minimum score for two reads to be in the same isoform cluster (default: 0.3)
    --iso-max-variance
        max allowed variance for two reads to be in the same isoform cluster (default: 25)
    -B, --bv-start-threshold
        starting threshold for the bitvector k-mer comparison (default: 0.4)
    -b, --bv-end-threshold
        ending threshold for the bitvector k-mer comparison (default: 0.2)
    -f, --bv-falloff
        falloff value for the bitvector threshold for each iteration (default: 0.05)
    -r, --min-reads-cluster
        minimum number of reads per cluster (default: 0)
    -p, --repr-percentile
        cluster representative percentile (default: 0.15)
    --rna
        use this mode if data is direct RNA (disables checking both strands)
    --verbose
        use this flag to print the progress of the run
    --raw
        set this flag to use all the reads without any length filtering (off by default)
    --lower-length
        filter out reads shorter than this value (default: 150)
    --upper-length
        filter out reads longer than this value (default: 100,000)
```

This clustering step will generate a file containing read clusters in binary format (clusters.out). To work with these clusters, the following commands are used.

### Description of clustering parameters 
**General parameters for the clustering step**
  
    --raw

If this flag is used, all reads from the input are used without any filtering for length, i.e. --lower-length and –upper-length parameters below are not used

    --lower-length (default: 150)

By default, we do not use reads that are shorter than 150nt. This limit can be increase to produce longer transcript models. Nanopore short-read sequencing is becoming possible, so this lower bound could be lowered to enable the reference-free reconstruction of small non-coding RNAs.

    --upper-length (default: 100,000)

Although very long transcripts are possible, we generally found reads longer than 100,000 nt not to be reliable, possibly resulting from experimental artifacts. As data improves, this parameter can be relaxed to identify ultra-long transcripts. 

**Parameters related to the bitvector comparison in the Clustering step**

    -B, --bv-start-threshold (default: 0.4)

This threshold is the minimal bitvector score to consider two reads to be potentially in the same gene cluster. The bitvector score defined as the fraction of unique k-mers that two reads have in common over the maximum of unique k-mers in the two reads. If the score is above this threshold, the two reads are compared using the LIS similarity score (see below –score-threshold). This threshold is the minimal score used. RATTLE performs multiple iterations of this test with all reads starting at the value of “-B” and decreasing by a step of “-f” until the threshold of “-b”. These multiple iterations makes it possible to test all reads under various conditions for clustering. To see the result of changing this parameter, please see Table S11 from RATTLE’s paper.

    -b, --bv-end-threshold (default: 0.2)

The ending threshold for the bitvector score in the iterations. A low value for -b makes possible to rescue reads have not been clustered in previous iterations. If -b is close to -B (or the same) only one or few iterations will be performed. This will make the clustering less sensitive, potentially resulting in many unclustered reads. To see the result of changing this parameter, please see Table S11 from RATTLE’s paper.

    -f, --bv-falloff (default: 0.05)

This is the step-change value between the first (-B) and final (-b) bitvector score thresholds and determines the number of iterations to perform clustering. A small value will provide more resolution in the definition of clusters but will result in more iterations, potentially leading to longer computational time. To see the result of changing this parameter, please see Table S11 from RATTLE’s paper.

    -r, --min-reads-cluster (default: 0)

Only clusters with more than this number of reads will be reported and used in the next step. The default means that also singletons (clusters composed of 1 single read) are also included.

    -p, --repr-percentile (default: 0.15)

In the iterative algorithm for clustering, reads are tested against representative of a cluster, rather than all the reads from that cluster. The value of -p is the position percentile position of the read in the ranking of reads sorted by length (from longest to shortest) in a cluster that is used as representative. The smaller the value, the closer to the top of the ranking. The longes read in a cluster may seem to be a better representative. However, during RATTLE optimization, we observed that this is not always the case, and using one few positions below (0.15 percentile, i.e. position 15th in a cluster of 100 reads) results in better performance.

**Parameters related to the LIS-based similarity in the Clustering step**

    -k, --kmer-size (default: 10, maximum: 16)

This is the size of k-mer used to compare two reads using the Longest Increasing Subsequence (LIS) algorithm (see RATTLE’s paper for details). A low value will enable a more sensitive comparison but will result in longer computing times. A higher value will make the comparison faster, but may miss cases due to sequencing errors. For reads with low error rate, this can be set to a value higher than the default. The maximum of 16 is used to ensure the efficiency of the algorithmic implementation.

    -s, --score-threshold (default: 0.2)

This parameter sets the minimum score threshold to include a read in a cluster. The similarity score is defined from the Longest Increasing Subsequence (LIS) algorithm as the number of bases covered by co-linearly matching k-mers (defined by –kmer-size) over the length of the shortest read in the pair. The score is calculated between the test read and the representative from the cluster. As any other similarity thresholds, a higher value will be more restrictive and result in more clusters, whereas a lower value would be more permissive and result in fewer clusters.

    -v, --max-variance (default: 1000000)

Two reads could present matching co-linear k-mers (defined by --kmer-size) that have gaps relative to each other. This parameter sets the maximum variance in the distribution of these gap lengths allowed to group two reads into the same gene cluster. It is set to be very large to enable clustering of reads into gene-like structures despite them having large gap differences. This first coarse-grained step can be refined in each gene cluster if the option “--iso” is used (see below). 

**Parameters related to the generation of transcript clusters**

    --iso

When this flag is used, RATTLE will process each gene cluster to determine whether the reads in the cluster can be separated into transcript clusters based on the properties of the co-linear matching k-mers. 

    --iso-kmer-size (default: 11, maximum: 16)

This is the size of k-mer used to compare two reads in the same gene cluster using the Longest Increasing Subsequence (LIS) algorithm to determine whether they should be split into separate transcript clusters. A low value will enable a more sensitive comparison but will result in longer computing times. A higher value will make the comparison faster, but may not be as accurate as separating reads that originate from different transcripts. The maximum of 16 is used to ensure the efficiency of the algorithmic implementation.

    --iso-score-threshold (default: 0.3)

Similar to –score-threshold but for the comparison of reads in the same gene cluster. This parameter sets the minimum score threshold to keep two reads in the same transcript cluster. The similarity score is defined from the Longest Increasing Subsequence (LIS) algorithm as the number of bases covered by co-linearly matching k-mers (defined by –kmer-size) over the length of the shortest read in the pair. The score is calculated between the test read and the representative from the cluster. As any other similarity thresholds, a higher value will be more restrictive and result in more clusters, whereas a lower value would be more permissive and result in fewer clusters.


    --iso-max-variance (default: 25)

Two reads belonging to the same gene cluster present matching co-linear k-mers (defined by --kmer-size) that may have gaps relative to each other. This parameter sets the maximum allowed variance in the distribution of these gap lengths. If the variance is above this threshold, the two reads are split into different transcript clusters. As described in RATTLE’s paper, a smaller value enables reconstruction of transcript isoforms that differ by smaller alternative exons. 



### Clustering summary
```
$ ./rattle cluster_summary -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -l, --label
        labels for the files in order of entry (optional)
    -c, --clusters
        clusters file. This is the output file from clustering step (required)
```

This command will output a CSV file to stdout containing two or three columns:
```
read_id, gene_cluster_id, [transcript_cluster_id]
```

### Cluster extraction
```
$ ./rattle extract_clusters -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -l, --label
        labels for the files in order of entry (optional)
    -c, --clusters
        clusters file. This is the output file from clustering step (required)
    -o, --output-folder
        output folder for fastx files (default: .)
    -m, --min-reads
        min reads per cluster to save it into a file
    --fastq
        whether input and output should be in fastq format (instead of fasta)
```

This command will generate 1 fastq/fasta file per cluster in the specified folder. This is useful to analyze particular clusters, perform alignments, etc.

## Error correction

This is the second main step in RATTLE. This command will perform error correction in each of the transcript clusters generated before. The command requires as input the previous clusters (`-c` option) and the original reads in fastq form (`-i` option), since the quality values are used in the correction step.

```
$ ./rattle correct -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -l, --label
        labels for the files in order of entry (optional)
    -c, --clusters
        clusters file. This is the output file from clustering step (required)
    -o, --output
        output folder (default: .)
    -g, --gap-occ
        gap-occ (default: 0.3)
    -m, --min-occ
        min-occ (default: 0.3)
    -s, --split
        split clusters into sub-clusters of size s for the multiple sequence alignment (default: 200)
    -r, --min-reads
        lower bound in the number of reads in a cluster to do the correction and produce a consensus (default: 5)
    -t, --threads
        number of threads to use (default: 1)
    --verbose
        use this flag to print the progress of the run
```

This step will generate 3 files:
* **corrected.fq**: contains reads that have been corrected by RATTLE
* **uncorrected.fq**: contains reads that have been left uncorrected by RATTLE. These will generally be those in clusters with fewer reads than specified with the option `-r`
* **consensi.fq**: contains one consensus sequence per transcript cluster

### Description of error correction parameters

    -g, --gap-occ (default: 0.3)

During the error correction step, RATTLE identifies small blocks (<10nt) followed by large gaps (larger than or equal to 20 positions) at both ends of each aligned read and removes them. A block is defined as a subsequence that has gaps occurring at a rate shorter to or equal to gap-occ (default 0.3, i.e. ~3nt in every 9-10 nt). RATTLE keeps removing blocks that satisfy these conditions at both ends of every aligned read until there no more such blocks are found

    -m, --min-occ (default: 0.3)

A consensus from each column in the MSA is then extracted in the following way: for each read and each base of the read, the base is changed to the consensus unless this base has an error  if the consensus base occurs with at least 60% frequency, but not if this base has an error probability (obtained from the FASTQ file) less than or equal to 1/3 times the average for the consensus base in that position of the alignment. Indels are treated similarly, but without the error constraint. This is only performed using aligned positions.


## Polishing step

In this last RATTLE step, transcript consensus sequences and gene definitions are reassessed to perform a cluster refinement using the same RATTLE metrics, and to obtain a final set of transcripts, consensus sequences (generated_from_transcript_clusters), abundances (total_reads) and labels.

```
$ ./rattle polish -h
    -h, --help
        shows this help message
    -i, --input
        input RATTLE consensi fasta/fastq file (required)
    -o, --output-folder
        output folder for fastx files (default: .)
    -l, --label
        labels for the files. if used, must be same as previous steps (optional)
    -t, --threads
        number of threads to use (default: 1)
    --rna
        use this mode if data is direct RNA (disables checking both strands)
    --verbose
        use this flag to print the progress of the run
```

Input must be the consensus sequences from the previous step.

# Example dataset

We provide below example datasets and rattle commands. We make available all input and output files in the folder **toyset/rna**. 

## Human direct RNA sequencing

We provide here as an example the analysis of a direct RNA sequencing dataset of 8306 reads. These are reads from an experiment in heart tissue used in the RATTLE manuscript and is available from the European Nucleotide Archive (ENA) under the umbrella study PRJEB40410 (https://www.ebi.ac.uk/ena/browser/view/PRJEB40410). The sample reads used in this example correspond to the subset that mapped to the human chromosome 20 with quality 60 (uniquely mapping reads). This input file is available at **./toyset/rna/input/sample.fastq**. We describe below the steps to build a reference-free transcriptome with RATTLE from these reads:

**cluster**
```
$ ./rattle cluster -i ./toyset/rna/input/sample.fastq -t 24 -o ./toyset/rna/output --rna
```
The output file generated from this step is clusters.out at ./toyset/rna/output/

**cluster_summary**
```
$ ./rattle cluster_summary -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out > ./toyset/rna/output/cluster_summary.tsv
```
The output generated from this step is a tsv file on the standard output which can be redirected to a file, for example ./toyset/rna/output/cluster_summary.tsv

**extract_clusters**
```
mkdir ./toyset/rna/output/clusters
$  ./rattle extract_clusters -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out -o ./toyset/rna/output/clusters --fastq
```
The output generated from this step is a clusters folder with 1 fastq file per cluster, for example ./toyset/rna/output/clusters

**correct**
```
$  ./rattle correct -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out -o ./toyset/rna/output/ -t 24
```
The output generated from this step are three files uncorrected.fq, corrected.fq and consensi.fq at ./toyset/rna/output/

**polish**
```
$  ./rattle polish -i ./toyset/rna/output/consensi.fq -o ./toyset/rna/output/ -t 24 --rna
```
The output generated from this step is the transcriptome.fq at ./toyset/rna/output/

**run time on example dataset** 

| command          | 1 thread      | 24 threads    |	
| :----------------|:-------------:| :-------------|
| cluster 	       | 0m21.302s     |0m8.457s       |
| cluster_summary  |  0m0.068s     |-              |
| extract_clusters |	0m0.079s	   |-              |  
| correct  		     | 1m16.180s	   | 0m15.939s     |
| polish  		     |  0m2.615s	   | 0m2.693s      |

# Reference based benchmarking

## ssCheck

We describe here the tool **ssCheck** for accuracy benchmarking of reads and transcripts mapped to a reference genome. ssCheck reads as input an alignment file in paf format and an annotation file in GTF format to compare with:

```
python ss_check.py [-h] [--beautiful] ref.gtf aln.paf

Calculate known/novel splice sites from PAF alignment and ref GTF file

positional arguments:
  ref.gtf      Reference GTF file (required)
  aln.paf      Alignment file in PAF format (required)

optional arguments:
  -h, --help   show this help message and exit
  --beautiful  Beautiful output (instead of csv)
```

ssCheck works by comparing annotation features in the mapped reads/sequences (exons, introns, exon-chains, intron-chains) with same features in the annotation and calculates the number of unique features as well as the total number of features predicted in the mapped reads. An intron-chain is defined as an ordered sequence of introns in an annotated transcript or mapped read. Similarly for exon-chains. The recall is calculated as the fraction of unique annotated features correctly found; precision is calculated as the fraction of unique predicted features that were in the annotation and read-precision is calculated as the fraction of the total number of predicted features in reads that corresponded to annotated features. Read-precision is affected by abundance levels but better reflects the accuracy per read. 

We developed ssCheck to be able to calculate the read-precision. Other methods merge identical intron-chains with different identifiers, which precludes this calculation.

## Example dataset

We provide below example datasets and ss_check.py commands. We make available all input and output files in the folder **toyset/sscheck**. The inputs are based on the same dataset described above for the RATTLE toyset. 
This input files are available at **./toyset/sscheck/input**. We describe below the steps to run ss_check.py:

**sam2paf**

We first use minimap2 paftools.js to convert the mapped reads in sam **format** to **paf** format. Please make sure minimap2 is installed and in yout PATH to use the below command.**
```
$ ./minimap2/misc/paftools.js sam2paf ./toyset/sscheck/input/sample.sam > ./toyset/sscheck/input/sample.paf
```
Both the files sample.sam and sample.paf are available at ./toyset/sscheck/input/

**ss_check with csv output**
```
$ python ./misc/ss_check.py  ./toyset/sscheck/input/sample_ref.gtf ./toyset/sscheck/input/sample.paf > ./toyset/sscheck/output/sample_output_sscheck.csv
```
The output generated from this step is a csv file for example ./toyset/sscheck/output/sample_output_sscheck.csv
The annotation file **sample_ref.gtf** is available at ./toyset/sscheck/input/ and corresponds to the GTF file selected for chr20 from Homo_sapiens.GRCh38.99.gtf downloaded from http://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/.

**ss_check with beautiful output**
```
$ python ./misc/ss_check.py --beautiful ./toyset/sscheck/input/sample_ref.gtf ./toyset/sscheck/input/sample.paf > ./toyset/sscheck/output/sample_output_sscheck.beautiful
```
The output generated in this case is easier to read. The output file for this command can be found in  ./toyset/sscheck/output/sample_output_sscheck.beautiful

For instance, for introns:

```
########################################
#             INTRON LEVEL             #
########################################
Introns in reference: 17769
Unique introns in reads: 1721
Reference introns found: 952/17769 (5.36%)
Total introns in reads: 15323
--> Known: 10121 (66.05%)
--> Novel: 5202 (33.95%)
```

## Cluster benchmark
We describe here the tool **cluster_benchmark** for cluster accuracy of predicted reads cluster label compared with reads mapped transciptome id. 

```
python cluster_benchmark.py [input] [output folder] [threads] [--rna]

positional arguments:
[input]  input fasta fie (required)
[output] output folder (required)

optional arguments:
threads  number of threads to run RATTLE
--rna    use this mode if data is direct RNA
```
Cluster_benchmark is designed to calculate the cluster accuracy scores (homogeneity, completeness) of the result from RATTLE clustering step. 

**example**  

We provide below example dataset and cluster_banchmark.py commands. We make available all input and output files under the folder ./toyset/cluster_benchmark. The inputs (sample.fa) are based on the same dataset described above for the RATTLE toyset with mapped transcriptome id label. And input (ref.fa) is a subset of the reference dataset downloaded from http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/. This input files are available at ./toyset/cluster_benchmark/input. We describe below the steps to run cluster benchmark.

**label the sample reads**  

We use minimap to label the sample reads with human transcriptome id. 

```
./minimap2 -ax map-ont --secondary=no ./toyset/cluster_benchmark/input/ref.fa ./toyset/rna/input/sample.fastq > ./toyset/cluster_benchmark/input/sample.sam
```
Then convert the sam file into a fasta file (sample.fa). The ref.fa and converted sample.fa file are both available at /toyset/cluster_benchmark/input folder. 

**cluster accuracy benchmark**  

We run the cluster_benchmark.py to benchmark the accuracy scores. The python file would automatically execute the RATTLE cluster and cluster summary step. Then use the summary tsv/csv file to calculate the accuracy score. 
```
python ./misc/cluster_benchmark.py ./toyset/cluster_benchmark/input/sample.fa ./toyset/cluster_benchmark/output 24 --rna
```

**result**  

This cluster accuracy benchmark result is generated from RATTLE v1.0.
| mode - reference    | homogeneity   | completeness   |
| :----------------|:-------------:| :-------------|
| RNA - transcriptome  | 84.56%    |92.22%     |
| RNA - genome  |  99.51%    |90.08%         |
| cDNA - transcriptome  | 84.52%   |92.21%     |
| cDNA - genome  |  99.47%    |90.07%         |

# Snakemake

Below we show how the user can run all steps in one go using Snakemake. We show this using the toy dataset as example. 
This can be adapted to any other dataset.
```
snakemake -s rattle_snakefile -p toyset/rna/snakemake_output/transcriptome.fq --cores 1
```