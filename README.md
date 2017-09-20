# Mini-metagenome Workshop September 25th 2017

Begin by logging into VM:

```
ssh -X ubuntu@137.205.69.49
```

Clone in the workshop repos:

```
mkdir ~/repos
cd repos
git clone https://github.com/chrisquince/WorkshopSept2017.git
```

Then we make ourselves a Projects directory:

```
mkdir ~/Projects
mkdir ~/Projects/AD
cd ~/Projects/AD
mkdir Reads
cd Reads
```

## Downloading the raw sequence reads

and download the anaerobic digester sequences:
```
cut -d"," -f7 ~/repos/WorkshopSept2017/data/metaFP1B.csv | sed '1d' > ForwardURL.txt
cut -d"," -f8 ~/repos/WorkshopSept2017/data/metaFP1B.csv | sed '1d' > ReverseURL.txt
```

```
while read line
do
    wget $line
done <  ForwardURL.txt
```

Can you work out how to download the reverse reads?

## Fastq file format

The reads are stored as pairs of fastq files.
```
head -n 10 S102_R1.fastq
```
Sometimes these will be zipped.

Lets have a look at the wikipedia page on [fastq format](https://en.wikipedia.org/wiki/FASTQ_format).

Count up number of reads in a fastq file:
```
cat S102_R1.fastq | echo $((`wc -l`/4))
```

What will the number be in S102_R2.fastq?

Discussion point: What are paired end reads?

Lets counts up reads in all files:
```
for file in *_R1.fastq; do cat $file | echo $((`wc -l`/4)); done > Counts.txt
```

And plot histogram, median read number is 4,917,354:
```
R
>library(ggplot2)
>Counts <- read.csv('Counts.txt',header=FALSE)
>summary(Counts$V1)
>pdf("Counts.pdf")
>qplot(Counts$V1, geom="histogram") 
>dev.off()
>q()
```

Now we run fastqc on one of the samples:
```
fastqc S102_R1.fastq
```

Look at the output files:
```
ls
firefox S102_R1_fastqc.html 
```

## Taxonomic profiling

For the taxonomic profiling we are going to subsample the fastq files to 1 million reads each 
for performance purposes.

```
cd ~/Projects/AD
mkdir ReadsSub
for file in Reads/*R1*fastq
do
    base=${file##*/}
    stub=${base%_R1.fastq}
    echo $stub
    seqtk sample -s100 $file 1000000 > ReadsSub/${stub}_Sub_R1.fastq
    seqtk sample -s100 Reads/${stub}_R2.fastq 1000000 > ReadsSub/${stub}_Sub_R2.fastq
done
```

We will use Kraken for profiling these reads but first lets convert them to interleaved fastq:

```
for file in ReadsSub/*R1*fastq
do
    
    stub=${file%_R1.fastq}
    echo $stub
    
    python ~/repos/WorkshopSept2017/scripts/Interleave.py $file ${stub}_R2.fastq ${stub}_R12.fastq
    
done

```

How does Kraken work?
![Kraken Figure1](Figures/KrakenFig.jpg)

Discussion point what is a kmer?

Now run kraken on the interleaved fastq:
```
mkdir Kraken
for file in ReadsSub/*R12*fastq
do
    base=${file##*/}
    stub=${base%_R12.fastq}
    echo $stub
    kraken --db ~/Databases/minikraken_20141208/ --threads 8 --preload --output Kraken/${stub}.kraken $file
done
```

Look at percentage of reads classified. Anaerobic digesters are under studied communities!

Discussion point what can we do about under representation in Database?

The output is just a text file:

```
head Kraken/S102_Sub.kraken
```

And we can generate a report:

```
kraken-report --db ~/Databases/minikraken_20141208/  Kraken/S102_Sub.kraken >  Kraken/S102_Sub.kraken.report
```

Some people prefer a different format:
```
kraken-mpa-report --db ~/Databases/minikraken_20141208/ Kraken/S102_Sub.kraken > Kraken/S102_Sub.kraken.mpa.report
```

We can get a report of the predicted genera:
```
cat  Kraken/S102_Sub.kraken.report | awk '$4=="G"'
```

Now lets get reports on all samples:
```
for file in Kraken/*.kraken
do
    stub=${file%.kraken}
    echo $stub
    kraken-report --db ~/Databases/minikraken_20141208/ $file >  ${stub}.kraken.report
done
```

Having done this we want to get one table of annotations at the genera level for community comparisons:

```
for file in Kraken/*.kraken.report
do
    stub=${file%.kraken.report}
    cat  $file | awk '$4=="G"' > $stub.genera
done
```

And then run associated script:
```
./CollateK.pl Kraken > GeneraKraken.csv
```

Put in NMDS plot

## Functional metagenome profiling


## Assembly based analysis

##Assembly based metagenomics analysis

We are now going to perform a basic assembly based metagenomics analysis of these same samples. 
This will involve a collection of different software programs:

1. megahit: A highly efficient metagenomics assembler currently our default for most studies

2. bwa: Necessary for mapping reads onto contigs

3. [samtools] (http://www.htslib.org/download/): Utilities for processing mapped files

4. CONCOCT: Our own contig binning algorithm

5. [prodigal] (https://github.com/hyattpd/prodigal/releases/): Used for calling genes on contigs

[gnu parallel] (http://www.gnu.org/software/parallel/): Used for parallelising rps-blast

[standalone blast] (http://www.ncbi.nlm.nih.gov/books/NBK52640/): Needs rps-blast

[COG RPS database] (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/): Cog databases

[GFF python parser] (https://github.com/chapmanb/bcbb/tree/master/gff)

#Co-assembly

We begin by performing a co-assembly of these samples using a program called megahit:

```
ls Reads/*R1.fastq | tr "\n" "," | sed 's/,$//' > R1.csv
ls Reads/*R2.fastq | tr "\n" "," | sed 's/,$//' > R2.csv
```

```
nohup megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 8 -o Assembly > megahit.out&
```

cat MetaTutorial/*R12.fasta > MetaTutorial/All_R12.fasta
megahit -r MetaTutorial/All_R12.fasta --presets meta -o Coassembly -t 8
We can have a look at how good the assembly was:

contig-stats.pl < Coassembly/final.contigs.fa


## Software installation

Going to make an installation directory:
```
mkdir Installation
cd Installation
```

1. Install R on the VM [R-base](https://www.r-bloggers.com/how-to-install-r-on-linux-ubuntu-16-04-xenial-xerus/):

```
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

and now some R packages:
```
sudo R
install.packages("ggplot2")
```

2. Install FastQC
```
sudo apt-get install fastqc
```
Requires a bug fix:
```
cd ~/Installation
sudo mkdir /etc/fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
sudo cp -r FastQC /etc/fastqc
```

3. Install viewer evince:
```
sudo apt install evince
```

4. Install html viewer firefox:
```
sudo apt install firefox
```

5. Install Metaphlan2:
```
sudo apt install mercurial
cd ~/Installation
hg clone https://bitbucket.org/biobakery/metaphlan2
```

Requires python numpy:
```
sudo apt-get install python-pip
sudo pip install numpy scipy
``
and also [BowTie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3)

6. Install Kraken

Get the mini-kraken database:
```
wget http://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
```

```
git clone https://github.com/DerrickWood/kraken.git
```

7. Seqtk

```
cd ~/Installation
git clone https://github.com/lh3/seqtk.git
cd seqtk; make
cp seqtk ~/bin/
```

8. Biopython
```
sudo apt-get update
sudo apt-get install python-biopython
```

9. Centrifuge
```
git clone https://github.com/infphilo/centrifuge.git
```

10. Megahit
```
git clone https://github.com/voutcn/megahit.git
cd megahit
make
./megahit -1 pe_1.fq.gz -2 pe_2.fq.gz -o megahit_out
cp megahit* ~/bin
```