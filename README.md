# Overview
This project scrapes Zurich rental property data from Homegate.ch, visualizes the data on a map, and performes an analysis to understand the characteristics and deciding factors of rental prices in Zurich. 

# Dependencies
The workflow manager snakemake handles the installation of the required dependencies into a local virtual environment. The only external dependencies you need are:

1. Install [anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) based on your operating system
2. Install [snakemake](https://snakemake.github.io/), ideally in its own separate conda virtual environment:
   ```bash
   conda create -c conda-forge -c bioconda -n snakemake snakemake
   ```

# Compiling
The steps to build the project are described in its snakemake file. The project can be compiled from scratch by running the snakemake command in its root directory:
```bash
cd /path/to/project_pp4rs
conda activate snakemake
snakemake --cores all --use-conda --conda-frontend conda
```

# Output
You can find the output with the main summary statistics and graphs of this analysis in ```pp4r_final_assignments/website.html```.






