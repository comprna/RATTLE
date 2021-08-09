# RATTLE
Reference-free reconstruction and quantification of transcriptomes from long-read sequencing

# Requirements
GCC with C++14 suppport

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

# Running RATTLE
**Warning**: The commands and parameters are still under development and may be subject to changes in future versions

## Sample commands
Details on each command can be found below, but these are some of the most common commands used when running RATTLE.

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
$ ./rattle correct -i reads.fq -c clusters.out -t 24 
```

* Polish RNA consensus sequences and build final transcriptome using 24 threads
```
$ ./rattle polish -i consensi.fq -t 24 --rna
```

## Clustering
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
        k-mer size for gene clustering (default: 10)
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

## Clustering summary
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

## Clustering extraction
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
        split clusters into sub-clusters of size s for msa (default: 200)
    -r, --min-reads
        min reads to correct/output consensus for a cluster (default: 5)
    -t, --threads
        number of threads to use (default: 1)
```

This command will perform error correction based on a previous clustering. Input must be fastq (quality values are used in the correction step).
This step will generate 3 files:
* **corrected.fq** -> contains reads that have been corrected by RATTLE
* **uncorrected.fq** -> contains reads that have been left uncorrected by RATTLE (for example those that are in clusters with fewer reads than specified with `-r`
* **consensi.fq** -> contains one consensus sequence per cluster

## Polishing step
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

Input must be the consensus sequences from the previous step. A final clustering and correction is performed to output a final transcriptome with quantification

## Example datasets

We provide an example dRNA file of 8306 reads to test Rattle in the folder named **toyset**

## Steps

**cluster**
```
$ ./rattle cluster -i ./toyset/rna/input/sample.fastq -t 24 -o ./toyset/rna/output --fastq --rna
```
**cluster_summary**
```
$ ./rattle cluster_summary -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out --fastq > ./toyset/rna/output/cluster_summary.tsv
```
**extract_clusters**
```
mkdir ./toyset/rna/output/clusters
$  ./rattle extract_clusters -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out -o ./toyset/rna/output/clusters --fastq 
```
**correct**
```
$  ./rattle correct -i ./toyset/rna/input/sample.fastq -c ./toyset/rna/output/clusters.out -o ./toyset/rna/output/ -t 24
```
**polish**
```
$  ./rattle polish -i ./toyset/rna/input/sample.fastq  -t 24 --rna
```
