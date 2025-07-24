# Tumor-immune-atlases
Reproducible pipeline for benchmarking tumour-immune reference atlases and evaluating supervised vs. unsupervised annotation strategies.

## Quick start 
```bash
git clone https://github.com/your-lab/ImmuneAtlasBenchmark.git
cd ImmuneAtlasBenchmark
conda env create -f environment.yml
conda activate immune-atlas
Rscript data/download_data.R
Rscript src/1_Atlas_similarity.R
Rscript src/2_Assessment_of_manually_labeled_datasets.R
Rscript src/3_Benchmark_metrics.R
Rscript src/4_immune_response_analysis.R
```
## Citation

