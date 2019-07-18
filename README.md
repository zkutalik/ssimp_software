[//]: ========================================
# Software for Summary Statistics Imputation
[//]: ========================================


This command-line software enables summary statistics imputation (SSimp) for GWAS summary statistics. 

The only input needed from the user are the **GWAS summary statistics** and a **reference panel** (e.g. 1000 genomes, needed for LD computation).


## Current version: 0.5.6
[//]: -------------------------------

If this is the **first time** that you download SSIMP run the following line:
```
ssimp --download.build.db
```

This will download a database with all positions on different builds.

Alternatively, run the following, if you are planning to use 1000 genomes as the reference panel.

```
ssimp --download.1KG
```
This will also download the build database along with 1000 genomes reference panel.


If you had previous versions of SSIMP installed, you **NEED to delete and redownload** the database:
```
rm ~/reference_panels/database.of.builds.1kg.uk10k.hrc.2018.01.18.bin
ssimp --download.build.db
```
The download option is equivalant to:
```
# wget https://drive.switch.ch/index.php/s/uOyjAtdvYjxxwZd/download -O ~/reference_panels/database.of.builds.1kg.uk10k.hrc.2018.01.18.bin
```

## Troubleshooting

You may get assertion errors when running this, even after following the *Installation* Instructions. Please also check the *Bug Reports* section below. Here are a couple of things you can check:
1) the size of the build database should be `1'789'839'360` bytes: `~/reference_panels/database.of.builds.1kg.uk10k.hrc.2018.01.18.bin`. If not, see *Bug Reports* below on how to download an updated version
2) The simplest way to call this is with: `ssimp gwas/small.random.txt 1KG/EUR output.txt`. Otherwise, you can select the any (super) population, e.g. AFR, with: `ssimp gwas/small.random.txt ~/reference_panels/1000genomes/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz output.txt --sample.names ~/reference_panels/1000genomes/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=AFR`. If you wish to impute only a particular chromosome, you should use `--impute.range 17`.


## Installation
[//]: -------------------------------


### Folder with compiled version (binary)

We recommend to download the full folder. This will give you access to all toy examples. 

(1) Download

`wget https://github.com/zkutalik/ssimp_software/archive/master.zip`

`unzip master.zip`

(2) Access the folder

`cd ssimp_software-master`

(3) Rename the binary

Depending on what OS you are on :

`cp compiled/ssimp-linux-0.5.6 ssimp`

`cp compiled/ssimp-osx-0.5.6 ssimp`

### Compiled version (binary)

[Linux](compiled/ssimp-linux-0.5.6) - static version

[MacOS](compiled/ssimp-osx-0.5.6) - dynamic version. You need to install GSL 1.16 (GNU Scientific Library) from here: [http://ftp.gnu.org/gnu/gsl/](http://ftp.gnu.org/gnu/gsl/).


### Compile from source 

Check your gcc version with `gcc --version`. It needs to be 5.0.0 or higher.

Mac user? You need to install GSL: `brew install gsl`.

(1) Download the zip file and unpack the zip file (or clone the github folder)

`wget https://github.com/zkutalik/ssimp_software/archive/master.zip`

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

See more examples [here](doc/examples.md).


## Documentation
[//]: -------------------------------
Run `ssimp` with no arguments to see the [usage message](doc/usage.txt). 

Check out [examples](doc/examples.md).

We also provide a [detailed manual](doc/manual.md) that has overlapping information with the usage message, but is a bit more detailed in more aspects.


## Download 1000 genomes reference panel
[//]: -------------------------------

Run 

```
ssimp --download.1KG
```

**Important!** *Reference panels provided in folder `ref` are toy reference panels for testing and examples and not made to use for proper usage.*

More info on handling reference panel data can be found in [usage message](doc/usage.txt).



## Run tests
[//]: -------------------------------
`stu @all.tests` or check [this file](tests/howto.txt). 

## Compile new version
[//]: -------------------------------

1. Follow instructions in [this file](compiled/howto.txt).
2. Make a pull request and merge it.
3. Add tag with `git pull master && git tag v0.5.6 master && git push --tags`

## Contributors
[//]: -------------------------------
[Aaron McDaid](https://github.com/aaronmcdaid) (implementation, method development)

[Sina R&uuml;eger](https://github.com/sinarueeger) (coordination, method development)

[Zoltan Kutalik](https://github.com/zkutalik) (supervision, method development)

Ninon Mounier (Bugfix [#104](https://github.com/zkutalik/ssimp_software/pull/104/files))

We have used code (with permission, and under the GPL) from [libStatGen](https://genome.sph.umich.edu/wiki/C%2B%2B_Library:_libStatGen) and [stu](https://github.com/kunegis/stu).

## Citation
[//]: -------------------------------

If you use this software in your scientific research, please consider citing

Rüeger, S., McDaid, A., Kutalik, Z. (2018). Evaluation and application of summary statistic imputation to discover new height-associated loci. PLOS Genetics. https://doi.org/10.1371/journal.pgen.1007371

Rüeger, S., McDaid, A., Kutalik, Z. (2018). Improved imputation of summary statistics for realistic settings. bioRxiv. https://doi.org/10.1101/203927 

Copyright 2018.
(Aaron McDaid, Sina Rüeger, Zoltán Kutalik, https://github.com/zkutalik/ssimp_software)

## Bug report
[//]: -------------------------------
Report bugs here: [https://github.com/zkutalik/ssimp_software/issues](https://github.com/zkutalik/ssimp_software/issues)

If you are having trouble with SSIMP versions >= 0.4, make sure that you delete your database in folder `~/reference_panels/database*` and redownload it using the following command

```
wget https://drive.switch.ch/index.php/s/uOyjAtdvYjxxwZd/download -O ~/reference_panels/database.of.builds.1kg.uk10k.hrc.2018.01.18.bin
```

<!--## Wiki / FAQ
[//]: -------------------------------
We are currently working on a Wiki.-->


## Contact
[//]: -------------------------------
<zoltan.kutalik@gmail.com>

