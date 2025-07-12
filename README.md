# Tumor-immune-atlases
Reproducible pipeline for benchmarking tumour-immune reference atlases and evaluating supervised vs. unsupervised annotation strategies.

## Quick start
```bash
git clone https://github.com/your-lab/ImmuneAtlasBenchmark.git
cd ImmuneAtlasBenchmark
conda env create -f environment.yml
conda activate immune-atlas
Rscript data/download_data.R
Rscript src/01_QC_and_preprocessing.R
Rscript src/02_Atlas_similarity.R
Rscript src/03_Supervised_annotation.R
Rscript src/04_Benchmark_metrics.R
Rscript src/05_Immunotherapy_analysis.R
## Citation
