
# LDSC (LD SCore) `v2.2.0`

`ldsc` is a command line tool for estimating heritability and genetic correlation from GWAS summary statistics. `ldsc` also computes LD Scores.

## Getting Started

There are several ways to install ldsc (v2.0.2). Docker and Pip are preferred but the user can also install from the github repository. Each of these ways are explained below:

### Docker/Singularity (preferred method)

ldsc is available as an image on docker and has been built for amd64 and arm64. This method of installation has all required dependencies install already. This image can be pulled using the following command:

```
docker pull jtb114/ldsc-test:latest
```

The working directory of this image is "app". There is a subdirectory called "ldsc" that has all the scripts used by ldsc such as "ldsc.py", "munge_sumstats.py", and "make_annot.py" (You will have to mount the appropriate data directories to the images). Users can then run this image in interactive mode using the following command:

```
docker run -it jtb114/ldsc-test
```

Sometimes in HPC environments, other containerization software is preferred to docker. One option is Singularity. Singularity is able to pull images from Dockerhub so the following command will pull the ldsc image and build a singularity image:

```
singularity pull docker://jtb114/ldsc-test
```

This command will create a file called ldsc-test-latest.sif in the directory in which the command is run. Users can then run the singularity image using the following command:

```
singularity shell ldsc-test-latest.sif
```

This will open the container in an interactive mode. Users can find the scripts for ldsc in the root directory but running a command such as:

```
ls /app/ldsc/
```

### Pip:
ldsc is also on PYPI and can be installed by Pip using the following command (It is preferable that the program be installed into a virtual environment created using venv or conda):

```
pip install ldsc
```

*warning:* The "make_annot.py" script requires that bedtools be installed. There have been issues when running ldsc on an HPC environment when bedtools has not been detect even when loaded/installed. If this error occurs, it is required to use the docker image to run ldsc.

### Github:
In order to download `ldsc`, you should clone this repository via the commands
```  
git clone https://github.com/belowlab/ldsc.git
cd ldsc
```

In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager or you can use pip. *warning*: you will have to install bedtools into the system python because this is not a python dependency. 

***If using anaconda:***
After installing Anaconda, run the following commands to create an environment with LDSC's dependencies:

```
conda env create --file environment.yml
source activate ldsc
```

***If using pip in a virtual environment (ldsc has been tested with python 3.9)***
First create the virtualenv using python. I have specified the python version in my command but you do not have to 

```
python3.9 -m venv venv
```

Nest use the requirements.txt file to install the correct requirements after activating the environment.

```
source venv/bin/activate

pip install -r requirements.txt
```

***after install dependencies using pip/anaconda***
Once the above has completed, you can run:

```
./ldsc.py -h
./munge_sumstats.py -h
```
to print a list of all command-line options. If these commands fail with an error, then something as gone wrong during the installation process. 

Short tutorials describing the four basic functions of `ldsc` (estimating LD Scores, h2 and partitioned h2, genetic correlation, the LD Score regression intercept) can be found in the wiki. If you would like to run the tests, please see the wiki.

*warning:* The "make_annot.py" script requires that bedtools be installed. There have been issues when running ldsc on an HPC environment when bedtools has not been detect even when loaded/installed. If this error occurs, it is required to use the docker image to run ldsc.

## Updating LDSC

If using docker/singularity, you can update ldsc to a new version by by pulling the new image down. If installed from github. You can update to the newest version of `ldsc` using `git`. First, navigate to your `ldsc/` directory (e.g., `cd ldsc`), then run
```
git pull
```
If `ldsc` is up to date, you will see 
```
Already up-to-date.
```
otherwise, you will see `git` output similar to 
```
remote: Counting objects: 3, done.
remote: Compressing objects: 100% (3/3), done.
remote: Total 3 (delta 0), reused 0 (delta 0), pack-reused 0
Unpacking objects: 100% (3/3), done.
From https://github.com/bulik/ldsc
   95f4db3..a6a6b18  master     -> origin/master
Updating 95f4db3..a6a6b18
Fast-forward
 README.md | 15 +++++++++++++++
 1 file changed, 15 insertions(+)
 ```
which tells you which files were changed. If you have modified the `ldsc` source code, `git pull` may fail with an error such as `error: Your local changes to the following files would be overwritten by merge:`. 

In case the Python dependencies have changed, you can update the LDSC environment with

```
conda env update --file environment.yml
```

## Where Can I Get LD Scores?

You can download [European](https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2) and [East Asian LD Scores](https://data.broadinstitute.org/alkesgroup/LDSCORE/eas_ldscores.tar.bz2) from 1000 Genomes [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/). These LD Scores are suitable for basic LD Score analyses (the LD Score regression intercept, heritability, genetic correlation, cross-sex genetic correlation). You can download partitioned LD Scores for partitioned heritability estimation [here](http://data.broadinstitute.org/alkesgroup/LDSCORE/).


## Support

Before contacting us, please try the following:

1. The [wiki](https://github.com/bulik/ldsc/wiki) has tutorials on [estimating LD Score](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial), [heritability, genetic correlation and the LD Score regression intercept](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) and [partitioned heritability](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability).
2. Common issues are described in the [FAQ](https://github.com/bulik/ldsc/wiki/FAQ)
2. The methods are described in the papers (citations below)

If that doesn't work, you can get in touch with us via the [google group](https://groups.google.com/forum/?hl=en#!forum/ldsc_users).

Issues with LD Hub?  Email ld-hub@bristol.ac.uk


## Citation

If you use the software or the LD Score regression intercept, please cite

[Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)

For genetic correlation, please also cite

[Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3406) Preprint available on bioRxiv doi: http://dx.doi.org/10.1101/014498

For partitioned heritability, please also cite

[Finucane, HK, et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3404) Preprint available on bioRxiv doi: http://dx.doi.org/10.1101/014241

For stratified heritability using continuous annotation, please also cite

[Gazal, S, et al. Linkage disequilibrium–dependent architecture of human complex traits shows action of negative selection. Nature Genetics, 2017.](https://www.nature.com/articles/ng.3954) 

If you find the fact that LD Score regression approximates HE regression to be conceptually useful, please cite

Bulik-Sullivan, Brendan. Relationship between LD Score and Haseman-Elston, bioRxiv doi: http://dx.doi.org/10.1101/018283

For LD Hub, please cite

[Zheng, et al. LD Hub: a centralized database and web interface to perform LD score regression that maximizes the potential of summary level GWAS data for SNP heritability and genetic correlation analysis. Bioinformatics (2016)](https://doi.org/10.1093/bioinformatics/btw613)


## License

This project is licensed under GNU GPL v3.


## Authors

Brendan Bulik-Sullivan (Broad Institute of MIT and Harvard)

Hilary Finucane (MIT Department of Mathematics)
