# Mini-metagenome Workshop September 25th 2017

## Getting started

Begin by logging into VM:

```
ssh -X ubuntu@137.205.69.49
```

Add this line to .profile using vi:
```
PATH=$HOME/repos/WorkshopSept2017/scripts:$PATH
```

Discussion point environment variables and configuration files.

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
    seqtk sample -s100 $file 1000000 > ReadsSub/${stub}_Sub_R1.fastq&
    seqtk sample -s100 Reads/${stub}_R2.fastq 1000000 > ReadsSub/${stub}_Sub_R2.fastq&
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

## Functional gene profiling

To perform functional gene profiling we will use Diamond to map against the KEGG database. 
First we will set an environmental variable to point to our copy of the Kegg:

export KEGG_DB=~/Databases/keggs_database/KeggUpdate/
mkdir KeggD
for file in MetaTutorial/{C,H}*R12.fasta
do 
   
   stub=${file%_R12.fasta}
   stub=${stub#MetaTutorial\/}
   echo $stub
   if [ ! -f KeggD/${stub}.m8 ]; then
    echo "KeggD/${stub}.m8"
    diamond blastx -d $KEGG_DB/genes/fasta/kegg_genes_dmd -q $file -p 8 -a KeggD/${stub}.dmd
    diamond view -a KeggD/${stub}.dmd -o KeggD/${stub}.m8
   fi
done
This is a slow process even using the very efficient Diamond aligner. We recommend therefore stopping the above process and copying across the prerun samples:

rm -r KeggD
cp -r ~/Archive/KeggD .
Having mapped reads to the KEGG genes we can collate these into ortholog coverages:

for file in KeggD/*.m8
do
    stub=${file%.m8}

    echo $stub
    
    python ~/bin/CalcKOCov.py $file $KEGG_DB/ko_genes_length.csv $KEGG_DB/genes/ko/ko_genes.list > ${stub}_ko_cov.csv

done
We collate these into a sample table:

mkdir FuncResults
Collate.pl KeggD _ko_cov.csv KeggD/*_ko_cov.csv > FuncResults/ko_cov.csv
and also KEGG modules:

for file in KeggD/*ko_cov.csv
do
    stub=${file%_ko_cov.csv}

    echo $stub
    python ~/bin/MapKO.py $KEGG_DB/genes/ko/ko_module.list $file > ${stub}_mod_cov.csv 
done
Collate those across samples:

Collate.pl KeggD _mod_cov.csv KeggD/*_mod_cov.csv > FuncResults/mod_cov.csv
It turns out that neither the Kegg orthologs or modules are significant between the two groups:

SigTest.R -c FuncResults/ko_cov.csv -m MetaTutorial/Meta.csv
SigTest.R -c FuncResults/mod_cov.csv -m MetaTutorial/Meta.csv



## Assembly based metagenomics analysis

We are now going to perform a basic assembly based metagenomics analysis of these same samples. 
This will involve a collection of different software programs:

1. megahit: A highly efficient metagenomics assembler currently our default for most studies

2. bwa: Necessary for mapping reads onto contigs

3. [samtools] (http://www.htslib.org/download/): Utilities for processing mapped files

4. CONCOCT: Our own contig binning algorithm

5. [prodigal] (https://github.com/hyattpd/prodigal/releases/): Used for calling genes on contigs

6. [gnu parallel] (http://www.gnu.org/software/parallel/): Used for parallelising rps-blast

7. [standalone blast] (http://www.ncbi.nlm.nih.gov/books/NBK52640/): Needs rps-blast

8. [COG RPS database] (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/): Cog databases

9. [GFF python parser] (https://github.com/chapmanb/bcbb/tree/master/gff)

### Co-assembly

We begin by performing a co-assembly of these samples using a program called megahit:

```
ls ReadsSub/*R1.fastq | tr "\n" "," | sed 's/,$//' > R1.csv
ls ReadsSub/*R2.fastq | tr "\n" "," | sed 's/,$//' > R2.csv
```

```
nohup megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 8 -o Assembly > megahit.out&
```

```
contig-stats.pl < Assembly/final.contigs.fa
```

Should see results like:
```
sequence #: 469120	total length: 412545660	max length: 444864	N50: 1124	N90: 375
```

Discussion point what is N50?

If the assembly takes too long download the results instead:
```
mkdir Assembly
cd Assembly
wget https://septworkshop.s3.climb.ac.uk/final.contigs.fa
cd ..
```

### Read mapping

Then cut up contigs and place in new dir:

```bash
cd ~/DesmanExample/Example
mkdir contigs
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m Assembly/final.contigs.fa > contigs/final_contigs_c10K.fa
```

Having cut-up the contigs the next step is to map all the reads from each sample back onto them. First index the contigs with bwa:

```bash
cd contigs
bwa index final_contigs_c10K.fa
cd ..
```

Then perform the actual mapping you may want to put this in a shell script:

```bash
mkdir Map

for file in *R1.fastq
do 
   
   stub=${file%_R1.fastq}

   echo $stub

   file2=${stub}_R2.fastq

   bwa mem -t 32 contigs/final_contigs_c10K.fa $file $file2 > Map/${stub}.sam
done
```

Discussion point how do mappers differ from aligners? Can we list examples of each?

How does (this)[https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform]  help DNA sequence analysis!

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

6. Install Kraken. Get the mini-kraken database:
    ```
    wget http://ccb.jhu.edu/software/kraken/dl/minikraken.tgz   
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
    cp megahit* ~/bin
    ```
11. [bwa](https://github.com/lh3/bwa): Necessary for mapping reads onto contigs
    ```
    cd ~/repos
    git clone https://github.com/lh3/bwa.git
    cd bwa; make
    cp bwa ~/bin
    ```

12. [bam-readcount](https://github.com/genome/bam-readcount): Used to get per sample base frequencies at each position

    ```
    cd ~/repos
    sudo apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev patch
    git clone https://github.com/genome/bam-readcount.git
    mkdir bam-readcount-build
    cd bam-readcount-build/
    cmake ../bam-readcount
    make
    cp bin/bam-readcount ~/bin/
    ```

13. [samtools](http://www.htslib.org/download/): Utilities for processing mapped files. The version available through apt will *NOT* work instead...

    ```
    cd ~/repos
    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
    tar xvfj samtools-1.3.1.tar.bz2 
    cd samtools-1.3.1/ 
    sudo apt-get install libcurl4-openssl-dev libssl-dev
    ./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.3.1
    make all plugins-htslib
    cp samtools ~/bin/  
    ```

14. [bedtools](http://bedtools.readthedocs.io/en/latest/): Utilities for working with read mappings

    ```
    sudo apt-get install bedtools
    ```

15. [prodigal](https://github.com/hyattpd/prodigal/releases/): Used for calling genes on contigs

    ```
    wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux 
    cp prodigal.linux ~/bin/prodigal
    chmod +rwx ~/bin/prodigal
    ```

16. [gnu parallel](http://www.gnu.org/software/parallel/): Used for parallelising rps-blast

    ```
    sudo apt-get install parallel
    ```

17. [standalone blast](http://www.ncbi.nlm.nih.gov/books/NBK52640/): Need a legacy blast 2.5.0 which we provide as a download:

    ```
    wget https://desmandatabases.s3.climb.ac.uk/ncbi-blast-2.5.0+-x64-linux.tar.gz
    
    tar -xvzf ncbi-blast-2.5.0+-x64-linux.tar.gz
    
    cp ncbi-blast-2.5.0+/bin/* ~/bin
    ```
    
18. [diamond](https://github.com/bbuchfink/diamond): BLAST compatible accelerated aligner

    ```
    cd ~/repos
    mkdir diamond
    cd diamond
    wget http://github.com/bbuchfink/diamond/releases/download/v0.8.31/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz
    cp diamond ~/bin/
    ```
    
19. We then install both the [CONCOCT](https://github.com/BinPro/CONCOCT) and 
[DESMAN]((https://github.com/chrisquince/DESMAN)) repositories. 
These are both Python 2.7 and require the following modules:
    ```
        sudo apt-get -y install python-pip
        sudo pip install cython numpy scipy biopython pandas pip scikit-learn pysam bcbio-gff
    ```

    They also need the GSL
    ```
        sudo apt-get install libgsl2-dev
        sudo apt-get install libgsl-dev
    ```

    Then install the repos and set their location in your .bashrc:
    ```
    cd ~/repos

    git clone https://github.com/BinPro/CONCOCT.git

    cd CONCOCT
    
    git fetch
    
    git checkout SpeedUp_Mp
    
    sudo python ./setup.py install
    
    ```
    Then DESMAN
    ```
        cd ~/repos

        git clone https://github.com/chrisquince/DESMAN.git

        cd DESMAN

        sudo python ./setup.py install
    ```

    Then add this lines to .bashrc:

    ```
        export CONCOCT=~/repos/CONCOCT
        export DESMAN=~/repos/DESMAN
    ```
    