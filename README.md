# Genetic Correlation Visualizations for Cardiac Imaging Traits

Code to generate Figure 2 and Figure 3 from the systematic review: **"Genetic Architecture of Imaging Derived Cardiac Phenotypes and Their Association with Cardiomyopathies"**

## Figures Overview

### Figure 2: Trait-Trait Genetic Correlation Heatmap
Shows the pairwise genetic overlap between 18 cardiac imaging trait categories based on shared genes identified across 84 studies.

**What it shows:**
- Each cell represents the Jaccard similarity between two trait categories
- Jaccard similarity = (shared genes) / (total unique genes across both traits)
- Values range from 0 (no shared genes) to 1 (identical gene sets)
- Hierarchical clustering groups traits with similar genetic architecture

### Figure 3: Gene-Trait Association Heatmap
Visualizes pleiotropic genes (genes associated with multiple cardiac traits) and their specific trait associations.

**What it shows:**
- Rows = pleiotropic genes (genes associated with ≥2 traits)
- Columns = cardiac imaging trait categories
- Gray tiles = significant gene-trait association
- Genes are clustered by similar patterns of trait associations
- Genes are annotated by functional category (sarcomeric, regulatory, transcription, etc.)


## Requirements

### R packages
```r
install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", "dendextend"))
```

## Usage

```bash
# Clone the repository
git clone https://github.com/yourusername/cardiac-genetics-figures.git
cd cardiac-genetics-figures

# Generate Figure 2
Rscript scripts/generate_figure2.R

# Generate Figure 3
Rscript scripts/generate_figure3.R
```

## Data Source

Data extracted from 84 studies identified in our systematic review of genetic associations with imaging-derived cardiac phenotypes. See the full manuscript for details:

- **Protocol:** [OSF Registration](https://doi.org/10.17605/OSF.IO/VPKS2)
- **Manuscript:** [Link when available]

## Trait Categories

The 18 cardiac imaging trait categories analyzed:

1. LV structure (mass, volumes, geometry)
2. LV function (ejection fraction, strain, contractility)
3. LV geometry & shape
4. LV morphology (trabeculation, wall thickness)
5. LV segmental geometry
6. LV outflow tract
7. LV tissue characterization
8. Derived LV phenotypes
9. LA structure (volume, dimensions)
10. LA function (emptying fractions, strain)
11. LA-LV coupling
12. RA structure
13. RA function
14. RV structure
15. RV function
16. RV-LV coupling
17. Aortic structure
18. Cardiomyopathy outcome (HCM, DCM, HF)

## Citation

If you use this code or data, please cite:

```bibtex
@article{dhabalia2025cardiac,
  title={Genetic Architecture of Imaging Derived Cardiac Phenotypes and Their 
         Association with Cardiomyopathies: A Systematic Review},
  author={Dhabalia Ashok, A and Benesch Vidal, ML and Lackner, MK and 
          Blankenberg, S and Magnussen, C and Ziegler, A and Becher, PM},
  year={2025},
  doi={10.17605/OSF.IO/VPKS2}
}
```

## Contact

**Amra Dhabalia Ashok:** amra.dhabaliaashok@cardio-care.ch or amra.dhabalia@gmail.com


## License

CC-By Attribution-NonCommercial-NoDerivatives 4.0 International
