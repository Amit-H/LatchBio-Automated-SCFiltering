# Latch Bio - Automated Single Cell Filtering
LatchBio Implementation of Dropkick for automated single cell filtering

This workflow automates single cell filtering using [Dropkick](https://doi.org/10.1101/gr.271908.120) - [source](https://github.com/KenLauLab/dropkick).

You can find the deployed version of this workflow on the [Latch Console](https://console.latch.bio/explore/66656/info). This was built using the tutorial from the Dropkick repo and the docs from the Latch SDK.

## Basic Usage:

### Input

A directory containing 10x Single Cell data, which should contain:
1. barcodes.tsv
2. genes.tsv
3. matrix.mtx

Please also include a name for the run - e.g Human_PBMCs, Mouse_TH17

### Output

A HTML file in the root directory named RUN_NAME__SCDK_Report.html, which is a converted Jupyter notebook containing:
1. Quality Control plots
2. Score plots
3. Gene coefficients list and plots
4. Prinicipal Component Analyses Plots
5. UMAP comparisons to other methods

## Usage Notes:

- The cell cycle genes used are as defined in [Tirosh et al, Science, 2016](https://www.science.org/doi/10.1126/science.aad0501), these are not currently changable.
- The software packages that power this workflow (Scanpy & Dropkick) were initially designed to be used as interactive tools. To overcome this, this workflow programatically constructs notebooks using nbformat and nbconvert.

## TO DO:
1. Integration of the outputs of this workflow with the inputs of downstream analysis.
2. Cleaner output HTML file.

## Contributions
Please feel free to contribute to this repo and improve it! 
