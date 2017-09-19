# Mini-metagenome Workshop September 25th 2017

Begin by logging into VM:

```
ssh ubuntu@137.205.69.49
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

## Software installation
