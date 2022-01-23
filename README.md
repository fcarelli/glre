# Analysis of germline-specific elements

Full code to reproduce the analysis of C. elegans germline regulatory elements. It involves two separate steps:
- annotation of regulatory elements in C. elegans and C. briggsae using an adaptation of the relmapping pipeline by Jänes et al, 2018 ("relmapping" step)
- analysis of germline-specific regulatory elements ("glre" step)

To reproduce the analysis, clone this repository in your home directory
```
git clone https://github.com/fcarelli/glre.git
```
then create two conda environments using the relmapping.yml and the glre.yml files. Use the following commands to create them:
```
conda env create --file relmapping.yml
conda env create --file glre.yml
```

## Regulatory elements annotation (relmapping step)
To run the relmapping step, you will need to upload the ATAC-seq and the longCap-seq data in the glre/relmapping/samples directory. Then activate the relmapping environment and launch the analysis on C. elegans and C. briggsae (you can to add ```-n``` to test a snakemake command)
```
source activate relmapping_fnc
cd ~/glre/relmapping
snakemake --use-conda --cores 50 -s Snakefile.cb annot_cb --cluster sbatch
snakemake --use-conda --cores 50 -s Snakefile.ce annot_ce --cluster sbatch
```
If you plan to use different sets of data, modify the [config.yaml](relmapping/workflows/config.yaml) file to update the filenames

## germline-specific elements analysis (glre step)
After annotating the C. elegans and C. briggsae regulatory elements, you can reproduce the analysis of germline-specific regulatory elements described in Carelli et al \[in prep\]. 

You will first have to upload in the /glre/data folder the ChIP-seq and RNA-seq datasets deposited in GEO under the accession GSE192540. Then, you can launch the glre step
```
source activate glre
cd ~/glre
snakemake --cores 10
```
