[//]: ========================================
# Software for Summary Statistics Imputation
[//]: ========================================

This command-line software enables summary statistics imputation (SSimp) for GWAS summary statistics. 

The only input needed from the user are the **GWAS summary statistics** and a **reference panel** (e.g. 1000 genomes, needed for LD computation).

## Installation
[//]: -------------------------------

### Compiled version

Download for Linux 
* [ssimp v0 - Linux](compiled/ssimp_linux_v0.zip) (compiled on a CentOS 6.5)

### Compile from source 
(1) Download the zip file and unpack the zip file (or clone the github folder)

`wget https://github.com/sinarueeger/ssimp_software/archive/master.zip`
`unzip ssimp_software-master.zip`

(2) access the folder

`cd ssimp_software-master`

(3) run the makefile (source compilation)

`make`

(4) Run your first summary statistics imputation on a test file (uses a small toy reference panel)

`bin/ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt`

(5) To run your own summary statistics imputation, you need to have a reference panel ready. Check out the next section to see how to download and install a reference panel. 

## Example
[//]: -------------------------------

`bin/ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` 

will impute the Z-statistics, using a selected reference panel (see section above) and generate a file `output.txt`. The txt file assigned to `--gwas` contains at least the following columns: SNP-id, reference allele, risk allele, Z-statistic, and at least one row. The imputed summary statistics are stored in `output.txt`. 

## Documentation
[//]: -------------------------------
Run `ssimp` with no arguments to see the [usage message](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt). 

Check out [examples](https://github.com/sinarueeger/ssimp_software/blob/master/docu/examples.md).

We also provide a [detailed manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/manual.md) that contains information not present in the usage message.


## Download 1000 genomes reference panel
[//]: -------------------------------

**Important!** *Reference panels provided in folder `ref` are toy reference panels for testing and examples and not made to use for proper usage!*

### Automatical download
By running `ssimp` with a special argument - `1KG/SUPER-POP` - assigned to the reference panel option, it will automatically download 1000 genomes reference panel and use the specified super population for imputation. 

For example, if we want 1000 genomes to be downloaded and EUR population used for imputation, we type: 

`bin/ssimp gwas/small.random.csv 1KG/EUR output.txt` 

Make sure that:
1) you have `wget` installed and
2) no folder `ref_panels` exists.

### Manual download

Download the files in 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
to a directory on your computer. 

    cd ~                # to your home directory
    mkdir -p ref_panels/1000genomes
    cd       ref_panels/1000genomes
    wget -nd -r 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'

When using `ssimp` you can then pass this directory name to 'ssimp', and specify that only
as subset of individuals (here AFR) should be used:

`ssimp gwas.txt output.txt ~/ref_panels/1000genomes --sample.names super_pop=AFR`

More info on handling reference panel data can be found starting from line 49 in [usage message](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt).

