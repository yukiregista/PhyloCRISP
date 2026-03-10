# Supplementary Data

This directory contains curated data archives used in the manuscript.

```text
delivery_20260314_v2/
├── simulation.tar.gz
├── mammals_full.tar.gz
├── mammals_subsampled.tar.gz
└── hiv_full.tar.gz
```

Each archive includes its own `README.md` with dataset-specific details.

## Archives in this folder

- `simulation.tar.gz`: simulation data used in both bootstrap and Bayesian experiments, across four scenarios, with a true tree and sequence alignments for each replicate.
- `mammals_full.tar.gz`: full mammals data, including the full alignment and the main full-dataset reference/bootstrap trees.
- `mammals_subsampled.tar.gz`: 100 subsampled mammals replicates, each with a subsampled alignment and the corresponding NCBI reference subtree.
- `hiv_full.tar.gz`: full HIV alignment together with its maximum-likelihood reference tree.

## External data links

### Bayesian benchmark data

The Bayesian benchmark datasets (`Coal320` and `Yule400`) are not included here.
The original posterior tree sets (before subsampling) and corresponding true trees are available from the following dataset records:

- Berling L, Klawitter J, Bouckaert R, Xie W, Gavryushkin A, Drummond A (2024) Posterior tree sets from Serial coalescent with 320 taxa. The University of Auckland. Dataset. doi:10.17608/k6.auckland.25338313
- Berling L, Klawitter J, Bouckaert R, Xie W, Gavryushkin A, Drummond A (2024) Posterior tree sets from Yule with 400 taxa. The University of Auckland. Dataset. doi:10.17608/k6.auckland.25337818

Reference paper:

- Berling L, Klawitter J, Bouckaert R, Xie D, Gavryushkin A, Drummond AJ (2025) Accurate Bayesian phylogenetic point estimation using a tree distribution parameterized by clade probabilities. PLOS Computational Biology 21(2): e1012789. doi:10.1371/journal.pcbi.1012789

### Original HIV and mammals data

The original HIV and mammals data used in this study are from Lemoine et al. (2018).
The complete release, including full bootstrap tree sets, is available at:

- https://github.com/evolbioinfo/booster-workflows/releases/tag/v0.1.0

Reference:

- Lemoine F., Domelevo Entfellner J.B., Wilkinson E., Correia D., Dávila Felipe M., De Oliveira T., Gascuel O. 2018. Renewing Felsenstein’s phylogenetic bootstrap in the era of big data. Nature. 556:452–456.
