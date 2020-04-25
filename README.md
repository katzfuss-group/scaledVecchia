# scaledVecchia
Scaled Vecchia approximation of Gaussian processes (Katzfuss et al., 2020)

## Basic use
To use the scaled Vecchia method on your own data, download and source the file vecchia_scaled.R, including installing the R packages GpGp and GPvecchia listed at the top of that file. Some example code:

```{r}
source('vecchia_scaled.R')

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
For details on the functionality of fit_scaled() and predictions_scaled(), see the vecchia_scaled.R file.

## Reproducing results
The other R files here reproduce the figures and results from Katzfuss et al. (2020):

## Reference
<!---
[Katzfuss, M.... (2020). Title. *arXiv:20.02*.](https://arxiv.org/abs/...)
--->
