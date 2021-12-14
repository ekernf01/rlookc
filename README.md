## An application-oriented R package for constructing and validating knockoffs

The framework of *Model-X knockoffs* and the corresponding R package [`knockoff`](https://cran.r-project.org/web/packages/knockoff/index.html) together provide FDR control for subset selection in sparse regression models with multivariate Gaussian covariates (1). Model-X knockoffs have also been applied for structure learning via graphical lasso (2) by leaving out each variable in turn, and they have been mathematically extended to enable testing of composite hypotheses (3,4,5) and mixture models for covariates (6). 

To facilitate progress in applications, especially in causal modeling of transcription, this package provides efficient free software for leave-one-out knockoff generation with both simple and composite null hypotheses, plus an implementation of Gaussian mixture model knockoffs. We also include certain features for checking model assumptions and FDR calibration, notably the K-nearest neighbors exchangeability test from (7).

#### Getting started

To install, run this. (You need the [`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html) package.)

`devtools::install_github("ekernf01/rlookc")`

The most informative vignettes for new users are:

- `vignette_calibration.R` (how to check your models)
- `vignette_knockoff_construction.R` (how to generate knockoffs)

#### References

1. Candes, E., Fan, Y., Janson, L., & Lv, J. (2016). Panning for gold: Model-X knockoffs for high-dimensional controlled variable selection. arXiv preprint arXiv:1610.02351.
2. Zheng, Z., Zhou, J., Guo, X., & Li, D. (2018). Recovering the graphical structures via knockoffs. Procedia Computer Science, 129, 201-207.
3. Dai, R., & Barber, R. (2016, June). The knockoff filter for FDR control in group-sparse and multitask regression. In International Conference on Machine Learning (pp. 1851-1859). PMLR.
4. Multi-resolution localization of causal variants across the genome. M. Sesia, E. Katsevich, S. Bates, E. Candès, C. Sabatti
Nature Communications, 11, 1093 (2020). https://www.nature.com/articles/s41467-020-14791-2 
5. Katsevich, E., & Sabatti, C. (2019). Multilayer knockoff filter: Controlled variable selection at multiple resolutions. The annals of applied statistics, 13(1), 1.
6. Gimenez, J. R., Ghorbani, A., & Zou, J. (2019, April). Knockoffs for the mass: new feature importance statistics with false discovery guarantees. In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 2125-2133). PMLR.
7. Romano, Y., Sesia, M., & Candès, E. (2020). Deep knockoffs. Journal of the American Statistical Association, 115(532), 1861-1872.