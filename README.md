# RATTLE
Reference-free reconstruction and quantification of transcriptomes from long-read sequencing

* I. de la Rubia, J. A. Indi, S. Carbonell-Sala,  J. Lagarde,  M.M. Albà, E. Eyras. Reference-free reconstruction and quantification of transcriptomes from Nanopore long-read sequencing. bioRxiv doi: https://doi.org/10.1101/2020.02.08.939942


----------------------------
# Table of Contents
----------------------------

 * [Requirements](#requirements)
 * [Installation](#installation)
 * [Quick start](#quick-start)
 * [Running RATTLE](#running-rattle)
   * [Clustering](#clustering)
     * [Clustering summary](#clustering-summary)
     * [Cluster extraction](#cluster-extraction)
   * [Error correction](#error-correction)
   * [Polishing step](#polishing-step)
 * [Example datasets](#example-datasets)
   * [Human direct RNA sequencing](#human-direct-RNA-sequencing)
 * [Appendix: reference-based benchmarking](#appendix-reference-based-benchmarking)
   * [ssCheck](#sscheck)
   * [Example dataset](#example-dataset)
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
$ ./rattle cluster -i reads.fq -t 24 -o . --fastq
```

* Cluster cDNA Nanopore reads at isoform level with 24 threads
```
$ ./rattle cluster -i reads.fq -t 24 --fastq --iso
```

* Cluster RNA Nanopore reads at isoform level with 24 threads
```
$ ./rattle cluster -i reads.fq -t 24 --fastq --iso --rna
```

* View clustering summary (csv with read_id,cluster_id)
```
$ ./rattle cluster_summary -i reads.fq -c clusters.out --fastq
```

* Extract 1 fastq file per cluster in clusters/ folder
```
$ mkdir clusters
$ ./rattle extract_clusters -i reads.fq -c transcripts.out -o clusters --fastq
```

* Correct reads with 24 threads using isoform clusters
```
$ ./rattle correct -i reads.fq -c clusters.out -t 24 --fastq
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
        input fasta/fastq file (required)
    --fastq
        whether input and output should be in fastq format (instead of fasta)
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
        k-mer size for isoform clustering (default: 11)
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
```

This clustering step will generate a file containing read clusters in binary format (clusters.out). To work with these clusters, the following commands are used.

### Clustering summary
```
$ ./rattle cluster_summary -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -c, --clusters
        clusters file (required)
    --fastq
        whether input and output should be in fastq format (instead of fasta)

```

This command will output a CSV file to stdout containing two columns:
```
read_id, cluster_id
```

### Cluster extraction
```
$ ./rattle extract_clusters -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -c, --clusters
        clusters file (required)
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
    -c, --clusters
        clusters file (required)
    -o, --output
        output folder (default: .)
    -g, --gap-occ
        gap-occ (default: 0.3)
    -m, --min-occ
        min-occ (default: 0.3)
    -s, --split
        split clusters into sub-clusters of size s for the multiple sequence alignment (default: 200)
    -r, --min-reads
        min reads to correct/output consensus for a cluster (default: 5)
    -t, --threads
        number of threads to use (default: 1)
    --fastq
        whether input should be in fastq format (instead of fasta)
```

This step will generate 3 files:
* **corrected.fq**: contains reads that have been corrected by RATTLE
* **uncorrected.fq**: contains reads that have been left uncorrected by RATTLE. These will generally be those in clusters with fewer reads than specified with the option `-r`
* **consensi.fq**: contains one consensus sequence per transcript cluster

## Polishing step

In this last RATTLE step, transcript consensus sequences and gene definitions are reassessed to perform a cluster refinement using the same RATTLE metrics, and to obtain a final set of transcripts, consensus sequences, and abundances.

```
$ ./rattle polish -h
    -h, --help
        shows this help message
    -i, --input
        input RATTLE consensi fasta/fastq file (required)
    -o, --output-folder
        output folder for fastx files (default: .)
    -t, --threads
        number of threads to use (default: 1)
    --rna
        use this mode if data is direct RNA (disables checking both strands)
```

Input must be the consensus sequences from the previous step.

# Example datasets

We provide below example datasets and rattle commands. We make available all input and output files in the folder **toyset/rna**. 

## Human direct RNA sequencing

We provide here as an example the analysis of a direct RNA sequencing dataset of 8306 reads. These are reads from an experiment in heart tissue used in the RATTLE manuscript and is available from the European Nucleotide Archive (ENA) under the umbrella study PRJEB40410 (https://www.ebi.ac.uk/ena/browser/view/PRJEB40410). The sample reads used in this example correspond to the subset that mapped to the human chromosome 20 with quality 60 (uniquely mapping reads). This input file is available at **./toyset/rna/input/sample.fastq**. We describe below the steps to build a reference-free transcriptome with RATTLE from these reads:

**cluster**
```
$ ./rattle cluster -i ./toyset/rna/input/sample.fastq -t 24 -o ./toyset/rna/output --fastq --rna
```
The output file generated from this step is clusters.out at ./toyset/rna/output/

**cluster_summary**
```
$ ./rattle cluster_summary -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out --fastq > ./toyset/rna/output/cluster_summary.tsv
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
$  ./rattle correct -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out -o ./toyset/rna/output/ -t 24 --fastq
```
The output generated from this step are three files uncorrected.fq, corrected.fq and consensi.fq at ./toyset/rna/output/

**polish**
```
$  ./rattle polish -i ./toyset/rna/output/consensi.fq -o ./toyset/rna/output/  -t 24 --rna
```
The output generated from this step is the transcriptome.fq at ./toyset/rna/output/

**run time on example dataset** 

| command          | 1 thread      | 24 threads    |	
| :----------------|:-------------:| :-------------|
| cluster 	       |	0m21.302s    |0m8.457s       |
| cluster_summary  | 0m0.068s      |-              |
| extract_clusters |	0m0.079s	   |-              |  
| correct  		     | 1m16.180s	   | 0m15.939s     |
| polish  		     | 0m2.615s	     | 0m2.693s      |

# Appendix: reference-based benchmarking

## ssCheck

We describe here the tool **ssCheck** for the accuracy benchmarking of reads and transcripts mapped to a reference genome. ssCheck reads as input an alignment file in paf format and an annotation file in GTF format to compare with:

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

