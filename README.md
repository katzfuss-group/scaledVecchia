# scaledVecchia
R implementation of the scaled Vecchia approximation of Gaussian processes (Katzfuss et al., 2020)

## Basic use
To use the scaled Vecchia method on your own data, source the file *vecchia_scaled.R*, including installing the R packages *GpGp* (**>= 0.2.2**) and *GPvecchia* listed at the top of that file. Some example code:

```{r}
source('https://raw.githubusercontent.com/katzfuss-group/scaledVecchia/master/vecchia_scaled.R')

## generate some toy data
inputs=matrix(runif(200),ncol=2)
y=sin(rowSums(inputs*5))
inputs.test=matrix(runif(100),ncol=2)

## estimate parameters and make predictions using scaled Vecchia
fit=fit_scaled(y,inputs)
preds=predictions_scaled(fit,inputs.test)

## plot results
plot(rowSums(inputs.test),preds)
```
For details on the functionality of `fit_scaled()` and `predictions_scaled()`, see the *vecchia_scaled.R* file.

## Reproducing results
The other R files here reproduce the figures and results from Katzfuss et al. (2020):
- Figure 1: scaled_NN_illustration.R
- Figure 2: matern_comparison.R
- Figure 3: borehole_comparison.R
- Table 1: testfun_comparison.R
- Figure 4: satelliteDrag_comparison.R

## Reference
Katzfuss, M., Guinness, J., & Lawrence, E. (2020). Scaled Vecchia approximation for fast computer-model emulation. [*arXiv:2005.00386*](https://arxiv.org/abs/2005.00386).
