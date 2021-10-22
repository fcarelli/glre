# glre

Full code to reproduce the analysis of C. elegans germline regulatory elements. It involves two separate steps:
- annotation of regulatory elements in C. elegans and C. briggsae using an adaptation of the relmapping pipeline by Jänes et al, 2018 ("relmapping" step)
- analysis of germline-specific regulatory elements ("glre" step)

To reproduce the analysis, clone this repository on your system
```
git clone 
```
then create two conda environments using the relmapping.yml and the glre.yml files. Use the following commands to create them:
- conda env create --file relmapping.yml
- conda env create --file glre.yml

To run the relmapping step, you will need to upload the ATAC-seq and the longCap-seq data in the glre/relmapping/samples directory. Then activate the relmapping environment and launch the analysis on C. elegans and C. briggsae:
- conda activate relmapping_fnc
- cd glre/relmapping
- snakemake --use-conda --cores 50 -s Snakefile.cb annot_cb --cluster sbatch
- snakemake --use-conda --cores 50 -s Snakefile.ce annot_ce --cluster sbatch

If you plan to use different sets of data, modify the glre/relmapping/workflows/config.yaml file to update the filenames
