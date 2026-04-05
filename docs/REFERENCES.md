# References

Complete citation list for the TCGA-BRCA PAM50 Multi-Subtype Transcriptomics Pipeline.

---

## Data Sources

**UCSC Xena / TOIL RNA-seq Compendium**
Goldman, M. J., Craft, B., Hastie, M., Repečka, K., McDade, F., Kamath, A., ... & Haussler, D.
(2020). Visualizing and interpreting cancer genomics data via the Xena platform.
*Nature Biotechnology*, 38(6), 675–678. https://doi.org/10.1038/s41587-020-0546-8

**TCGA Breast Invasive Carcinoma**
Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours.
*Nature*, 490(7418), 61–70. https://doi.org/10.1038/nature11412

**GTEx**
GTEx Consortium. (2013). The Genotype-Tissue Expression (GTEx) project.
*Nature Genetics*, 45(6), 580–585. https://doi.org/10.1038/ng.2653

---

## Breast Cancer Biology and PAM50

**Molecular portraits of breast cancer**
Perou, C. M., Sørlie, T., Eisen, M. B., van de Rijn, M., Jeffrey, S. S., Rees, C. A., ...
& Botstein, D. (2000). Molecular portraits of human breast tumours.
*Nature*, 406(6797), 747–752. https://doi.org/10.1038/35021093

**Intrinsic subtypes and prognosis**
Sørlie, T., Perou, C. M., Tibshirani, R., Aas, T., Geisler, S., Johnsen, H., ... &
Børresen-Dale, A.-L. (2001). Gene expression patterns of breast carcinomas distinguish
tumor subclasses with clinical implications.
*Proceedings of the National Academy of Sciences*, 98(19), 10869–10874.
https://doi.org/10.1073/pnas.191367098

**PAM50 intrinsic subtype classifier**
Parker, J. S., Mullins, M., Cheang, M. C. U., Leung, S., Voduc, D., Vickery, T., ...
& Perou, C. M. (2009). Supervised risk predictor of breast cancer based on intrinsic subtypes.
*Journal of Clinical Oncology*, 27(8), 1160–1167. https://doi.org/10.1200/JCO.2008.18.1370

**Field cancerization in breast**
Deng, G., Lu, Y., Zlotnikov, G., Thor, A. D., & Smith, H. S. (1996). Loss of heterozygosity in
normal tissue adjacent to breast carcinomas. *Science*, 274(5295), 2057–2059.
https://doi.org/10.1126/science.274.5295.2057

---

## Statistical Methods

**limma — linear models and moderated t-statistics**
Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
*Nucleic Acids Research*, 43(7), e47. https://doi.org/10.1093/nar/gkv007

**Stability selection — ElasticNet bootstrap**
Meinshausen, N., & Bühlmann, P. (2010). Stability selection.
*Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 72(4), 417–473.
https://doi.org/10.1111/j.1467-9856.2010.01740.x

**glmnet — ElasticNet and LASSO**
Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized
linear models via coordinate descent.
*Journal of Statistical Software*, 33(1), 1–22. https://doi.org/10.18637/jss.v033.i01

**fgsea — fast gene set enrichment analysis**
Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artem Loboda, M. N. K., & Sergushichev, A.
(2021). Fast gene set enrichment analysis.
*bioRxiv*. https://doi.org/10.1101/060012

**MSigDB Hallmark gene sets**
Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P.
(2015). The Molecular Signatures Database Hallmark gene set collection.
*Cell Systems*, 1(6), 417–425. https://doi.org/10.1016/j.cels.2015.12.004

**KEGG pathway database**
Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes.
*Nucleic Acids Research*, 28(1), 27–30. https://doi.org/10.1093/nar/28.1.27

**msigdbr R package**
Dolgalev, I. (2022). msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format.
R package. https://igordot.github.io/msigdbr/

**WGCNA — weighted gene co-expression network analysis**
Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network
analysis. *BMC Bioinformatics*, 9, 559. https://doi.org/10.1186/1471-2105-9-559

**UMAP — Uniform Manifold Approximation and Projection**
McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform manifold approximation and
projection for dimension reduction. *arXiv preprint* arXiv:1802.03426.
https://arxiv.org/abs/1802.03426

---

## R Packages and Bioconductor

**Bioconductor**
Gentleman, R. C., Carey, V. J., Bates, D. M., Bolstad, B., Dettling, M., Dudoit, S., ...
& Zhang, J. (2004). Bioconductor: open software development for computational biology and
bioinformatics. *Genome Biology*, 5(10), R80. https://doi.org/10.1186/gb-2004-5-10-r80

**org.Hs.eg.db — human gene annotation**
Carlson, M. (2023). org.Hs.eg.db: Genome wide annotation for Human.
R package version 3.18.0. Bioconductor. https://bioconductor.org/packages/org.Hs.eg.db

**AnnotationDbi**
Pagès, H., Carlson, M., Falcon, S., & Li, N. (2023). AnnotationDbi: Manipulation of SQLite-based
annotations in Bioconductor. R package. https://bioconductor.org/packages/AnnotationDbi

**survival R package**
Therneau, T. M. (2023). A Package for Survival Analysis in R.
R package version 3.5-7. https://CRAN.R-project.org/package=survival

**survminer**
Kassambara, A., Kosinski, M., & Biecek, P. (2021). survminer: Drawing Survival Curves using
ggplot2. R package. https://CRAN.R-project.org/package=survminer

**uwot — UMAP in R**
Melville, J. (2023). uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for
Dimensionality Reduction. R package. https://CRAN.R-project.org/package=uwot

**UCSCXenaTools**
Wang, S., & Liu, X. (2019). The UCSCXenaTools R package: a toolkit for accessing genomics data
from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq.
*Journal of Open Source Software*, 4(40), 1627. https://doi.org/10.21105/joss.01627

**ComplexHeatmap**
Gu, Z. (2022). Complex heatmap visualization.
*iMeta*, 1(3), e43. https://doi.org/10.1002/imt2.43

**ggplot2**
Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York.
https://ggplot2.tidyverse.org

**renv — reproducible R environments**
Ushey, K., & Wickham, H. (2023). renv: Project Environments.
R package. https://CRAN.R-project.org/package=renv
