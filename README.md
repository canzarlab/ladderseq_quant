# kallisto-ls          

![Colorful reads in Ladder-seq](/LadderQuant400.png)



**Ladder-seq** is the concerted advancement of the RNA-seq protocol and its computational methods. It experimentally separates transcripts according to their length prior to sequencing to achieve a "coloring" of reads that points to its originating transcript.

In __kallisto-ls__ we have extended the Expectation-Maximization (EM) algorithm that was originally implemented in [kallisto](https://pachterlab.github.io/kallisto/) to use this coloring to more accurately quantify transcripts based on pseudoalignments of Ladder-seq reads. This software has been forked from the kallisto v0.44.0 [repository](https://github.com/pachterlab/kallisto).


## Installing kallisto-ls

__kallisto-ls__ can be installed through the CMake build manager:

```shell
git clone https://github.com/canzarlab/kallisto-ls.git
cd kallisto-ls
mkdir build
cd build
cmake ../
make
make install
```

## Running kallisto-ls

The main steps required to quantify the expression of transcripts from Ladder-seq data using __kallisto-ls__ are identical to conventional kallisto. A more detailed description of its usage can be found in the kallisto [manual](http://pachterlab.github.io/kallisto/manual.html).

__kallisto-ls__ assumes that Ladder-seq input reads for each band are stored in a common directory with the following file name convention: ```R1_band<i>.fastq``` and ```R2_band<i>.fastq``` for paired-end reads in band ```i```. __kallisto-ls__ currently assumes a separation of reads into 7 bands.

### Step 1: Building an index

First, an index needs to be built from a file of cDNA sequences in FASTA format.

```shell
kallisto-ls index -i <index_file> <cDNA_file.fasta>
```

### Step 2: Quantifying transcript abundances

__kallisto-ls__ estimates abundances of transcripts from Ladder-seq reads, which it outputs in file `abundances.tsv`. It expects the index computed in step 1 and (paired-end) reads from all 7 bands.

```shell
kallisto-ls quant -i <index_file> -o <output_directory> <R1_band1.fastq> <R2_band1.fastq> ... <R1_band7.fastq> <R2_band7.fastq>
```
The abundance estimation can also be run with a custom migration probabilities file which has to be provided in the command line using the flag -p.

Details on optional arguments and output files can be found in the kallisto [manual](https://pachterlab.github.io/kallisto/manual).


In addition to quantification of expression of transcripts, __kallisto-ls__ can also be used to calculate the migration probabilities from ladder-seq data which it outputs in file `migProb.tsv`.

```shell
kallisto-ls migProb -i <index_file> -o <output_directory> <R1_band1.fastq> <R2_band1.fastq> ... <R1_band7.fastq> <R2_band7.fastq>
```
