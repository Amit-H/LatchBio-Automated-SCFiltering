"""
Preprocess and Cluster 10x Single Cell Data
"""

from pathlib import Path
import subprocess
from latch import large_task, workflow
from latch.types import LatchFile, LatchDir
import nbformat as nbf


@large_task
def run_scanpy_analysis(run_name: str, data_dir: LatchDir) -> LatchFile:
    '''
    This function runs the scanpy clustering steps as outlined in scanpy_wrapper.py

    Input params:
    data_dir: the path to the 10x folder

    Output params:
    A HTML notebook of the qc plot, score plot, coefficients plot, PCA plot and UMAP plot
    '''

    data_folder = Path(data_dir).resolve()

    nb = nbf.v4.new_notebook()

    cells = [
        nbf.v4.new_markdown_cell(f"""# Automated {run_name} Filtering Report"""),
        nbf.v4.new_markdown_cell("""## Document Setup"""),
        nbf.v4.new_code_cell("""import scanpy as sc; sc.set_figure_params(color_map="viridis", frameon=False); import dropkick as dk""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell(f"""adata = sc.read_10x_mtx('{data_folder}', var_names='gene_symbols', cache=True)""", metadata = {"collapsed" : True}),
        nbf.v4.new_markdown_cell("""## Quality Control"""),
        nbf.v4.new_code_cell("""adata = dk.recipe_dropkick(adata, n_hvgs=None, X_final="raw_counts")""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""qc_plt = dk.qc_summary(adata)"""),
        nbf.v4.new_markdown_cell("""## Model Generation"""),
        nbf.v4.new_code_cell("""adata_model = dk.dropkick(adata, n_jobs=5)"""),
        nbf.v4.new_markdown_cell("""## Filtering"""),
        nbf.v4.new_code_cell("""adata_filtered = adata""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""adata_filtered = dk.recipe_dropkick(adata_filtered, X_final="arcsinh_norm", filter=True, n_hvgs=2000, verbose=True)""", metadata = {"collapsed" : True}),
        nbf.v4.new_markdown_cell("""## Cell Cycle Gene Loading"""),
        nbf.v4.new_code_cell("""ccg = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"]"""),
        nbf.v4.new_code_cell("""cc_genes = ccg""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""s_genes = cc_genes[:43]""",metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""g2m_genes = cc_genes[43:]""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""sc.tl.score_genes_cell_cycle(adata_filtered, s_genes=s_genes, g2m_genes=g2m_genes)""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""sc.tl.pca(adata_filtered, n_comps=50, use_highly_variable=True)""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""sc.pp.neighbors(adata_filtered, n_neighbors=30, random_state=1, n_pcs=10)""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""sc.tl.leiden(adata_filtered)""", metadata = {"collapsed" : True}),
        nbf.v4.new_code_cell("""sc.tl.umap(adata_filtered, random_state=1)""", metadata = {"collapsed" : True}),
        nbf.v4.new_markdown_cell("""## Outputs"""),
        nbf.v4.new_markdown_cell("""### Ambient Percentage VS Total Genes Detected Per Cell"""),
        nbf.v4.new_code_cell("""score_plt = dk.score_plot(adata)"""),
        nbf.v4.new_markdown_cell("""### Gene Coefficients"""),
        nbf.v4.new_code_cell("""dk.coef_inventory(adata)"""),
        nbf.v4.new_markdown_cell("""### Gene Coefficients & Cross Validation Scores"""),
        nbf.v4.new_code_cell("""coef_plt = dk.coef_plot(adata)"""),
        nbf.v4.new_markdown_cell("""### Prinicpal Component Analysis"""),
        nbf.v4.new_code_cell("""pca_plt = sc.pl.pca_overview(
            adata_filtered, 
            color=[
                "arcsinh_total_counts",
                "pct_counts_mito",
                "phase"])"""),
        nbf.v4.new_markdown_cell("""### UMAP"""),
        nbf.v4.new_code_cell("""umap_plt = sc.pl.umap(
            adata_filtered,
            color=[
                "arcsinh_total_counts",
                "pct_counts_mito",
                "phase",
                "leiden",
                "dropkick_label",
                "dropkick_score",
                ],
            legend_fontsize="large",
            ncols=4)""")
    ]
    nb.cells.extend(cells)
    
    with open('Report.ipynb', 'w') as f:
        nbf.write(nb, f)

    _html_convert_cmd  = [
        "jupyter",
        "nbconvert",
        "--to",
        "HTML",
        "--execute",
        "/root/Report.ipynb",
    ]

    subprocess.run(_html_convert_cmd)
    
    return LatchFile("/root/Report.html", f"latch:///{run_name}_SCDK_Report.html")


@workflow
def scanpyXdropkick(run_name: str, data_dir: LatchDir) -> LatchFile:
    """Description...

    Scanpy & Dropkick - Automated Single Cell Filtering  
    ----

    This workflow automates single cell filtering using [Dropkick](https://doi.org/10.1101/gr.271908.120) - [source](https://github.com/KenLauLab/dropkick).

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
    - The cell cycle genes used are as defined in [Tirosh et al, Science, 2016](https://doi.org/10.1126/science.aad0501), these are not currently changable.
    - The software packages that power this workflow (Scanpy & Dropkick) were initially designed to be used as interactive tools. To overcome this, this workflow programatically constructs notebooks using [nbformat](https://nbformat.readthedocs.io/en/latest/) and [nbconvert](https://nbconvert.readthedocs.io/en/latest/).

    ## Future work:
    1. Integration of the outputs of this workflow with the inputs of downstream analysis.
    2. Cleaner output HTML file.

    __metadata__:
        display_name: Preprocess and Cluster 10x Single Cell Data
        author:
            name: Amit Halkhoree
            email: amithalkhoree@gmail.com
            github: https://github.com/Amit-H
        repository: https://github.com/Amit-H/LatchBio-Automated-SCFiltering
        license:
            id: MIT

    Args:
        data_dir:
            Directory of data where .mtx file is stored.

            __metadata__:
                display_name: Data Directory
        
        run_name:
            Name of analyis - i.e cell types: PBMC, TH17 
            
            __metadata__:
                display_name: Run Name
                
    """
    from latch.resources.launch_plan import LaunchPlan

    LaunchPlan(
    scanpyXdropkick,
    "hg19-PBMCs",
    {"run_name": "human_PBMCs", "c": LatchFile("s3://latch-public/test-data/hg19")}
    )

    return run_scanpy_analysis(data_dir=data_dir, run_name = run_name)
