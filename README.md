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

We will use Kraken for profiling these reads but first lets convert them to interleaves fastq:

```

```

How does Kraken work?
![Kraken Figure1](Figures/KrakenFig.png)


Discussion point what is a kmer?

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