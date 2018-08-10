[//]: ========================================
# Software for Summary Statistics Imputation
[//]: ========================================


This command-line software enables summary statistics imputation (SSimp) for GWAS summary statistics. 

The only input needed from the user are the **GWAS summary statistics** and a **reference panel** (e.g. 1000 genomes, needed for LD computation).


## Current version: 0.4
[//]: -------------------------------

## Installation
[//]: -------------------------------


### Folder with compiled version (binary)

We recommend to download the full folder. This will give you access to all toy examples. 

(1) Download

`wget https://github.com/sinarueeger/ssimp_software/archive/master.zip`

`unzip master.zip`

(2) Access the folder

`cd ssimp_software-master`

(3) Rename the binary

Depending on what OS you are on:

`cp compiled/ssimp-linux-0.4 ssimp`

`cp compiled/ssimp-osx-0.4 ssimp`

### Compiled version (binary)

[Linux](compiled/ssimp-linux-0.4) - static version

[MacOS](compiled/ssimp-osx-0.4) - dynamic version. You need to install GSL 1.16 (GNU Scientific Library) from here: [http://ftp.gnu.org/gnu/gsl/](http://ftp.gnu.org/gnu/gsl/). 


### Compile from source 

Check your gcc version with `gcc --version`. It needs to be 5.0.0 or higher.

Mac user? You need to install GSL: `brew install gsl`.

(1) Download the zip file and unpack the zip file (or clone the github folder)

`wget https://github.com/sinarueeger/ssimp_software/archive/master.zip`

`unzip master.zip`

(2) Access the folder

`cd ssimp_software-master`

(3) Run the makefile (source compilation)

`make`

This will create a file `bin/ssimp`.

`cp bin/ssimp ssimp`

(4) Run your first summary statistics imputation on a test file (uses a small toy reference panel)

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt`

(5) To run your own summary statistics imputation, you need to have a reference panel ready (in the folder `reference_panels` in your home directory). Check out [Download 1000 genomes reference panel](#download-1000-genomes-reference-panel) below to download. 




## :warning: Size of the software 
[//]: -------------------------------


While the github folder itself is rather small (136 MB), the folder `reference_panels` takes up about 17 GB storage on the **hard disk** when 1KG is installed.

When running, the software uses about 5 GB of memory. Make sure you have enough RAM when running the software in parallel!


## Runtime
[//]: -------------------------------

It takes roughly ~200 CPU hours to impute to full genome (~20M variants) with 500 individuals from 1000 Genomes reference panel. Runtime depends on the number of variants imputed and the size of the reference panel. 


## Example
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` 

will impute the Z-statistics, using a selected reference panel (see section above) and generate a file `output.txt`. The txt file assigned to `--gwas` contains at least the following columns: SNP-id, reference allele, risk allele, Z-statistic, and at least one row. The imputed summary statistics are stored in `output.txt`.


## Documentation
[//]: -------------------------------
Run `ssimp` with no arguments to see the [usage message](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt). 

Check out [examples](https://github.com/sinarueeger/ssimp_software/blob/master/docu/examples.md).

We also provide a [detailed manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/manual.md) that contains information not present in the usage message.


## Download 1000 genomes reference panel
[//]: -------------------------------

**Important!** *Reference panels provided in folder `ref` are toy reference panels for testing and examples and not made to use for proper usage.*

### Automatic download
[//]: -------------------------------
By running `ssimp` with a special argument - `1KG/SUPER-POP` - assigned to the reference panel option, it will automatically download 1000 genomes reference panel and use the specified super population for imputation. 

Make sure that:
1) you have `wget` installed and
2) no folder `reference_panels` exists.

For example, if we want 1000 genomes to be downloaded and EUR population used for imputation, we type: 

`ssimp gwas/small.random.txt 1KG/EUR output.txt` 

It will then download 1000 genomes reference panel **and** a file that enables fast imputation of single SNPs. 

When using `ssimp` you can then pass this directory name to 'ssimp', and specify that only
as subset of individuals (here AFR) should be used:

`ssimp gwas/small.random.txt ~/reference_panels/1000genomes/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz output.txt --sample.names ~/reference_panels/1000genomes/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=AFR`

More info on handling reference panel data can be found starting from line 49 in [usage
message](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt).


## Run tests
[//]: -------------------------------
`stu @all.tests` or check [this file](tests/howto.txt). 

## Compile new version
[//]: -------------------------------

1. Follow instructions in [this file](compiled/howto.txt).
2. Update VERSION file.

## Contributors
[//]: -------------------------------
[Aaron McDaid](https://github.com/aaronmcdaid) (implementation, method development)

[Sina R&uuml;eger](https://github.com/sinarueeger) (coordination, method development)

[Zoltan Kutalik](https://github.com/zkutalik) (supervision, method development)

We have used code (with permission, and under the GPL) from [libStatGen](https://genome.sph.umich.edu/wiki/C%2B%2B_Library:_libStatGen) and [stu](https://github.com/kunegis/stu).

## Citation
[//]: -------------------------------

If you use this software in your scientific research, please consider citing

Rüeger, S., McDaid, A., Kutalik, Z. (2017). Evaluation and application of summary statistic imputation to discover new height-associated loci. PLOS Genetics. https://doi.org/10.1371/journal.pgen.1007371

Rüeger, S., McDaid, A., Kutalik, Z. (2017). Improved imputation of summary statistics for realistic settings. bioRxiv. https://doi.org/10.1101/203927

Copyright 2018.
(Aaron McDaid, Sina Rüeger, Zoltán Kutalik, https://github.com/zkutalik/ssimp_software)

## Bug report
[//]: -------------------------------
Report bugs here: (https://github.com/sinarueeger/ssimp_software/issues)[https://github.com/sinarueeger/ssimp_software/issues]


<!--## Wiki / FAQ
[//]: -------------------------------
We are currently working on a Wiki.-->


## Contact
[//]: -------------------------------
<zoltan.kutalik@gmail.com>

