# SpPDCC R package
An R package for estimating sparse and positive definite basis covariance matrices from compositional data using the method described in [Direct covariance matrix estimation with compositional data](https://arxiv.org/abs/2212.09833). 

This package will be updated sporadically. Please contact amolstad@umn.edu with any questions or comments. 

### Installation
SpPDCC can be loaded directly into R through the `devtools` package:
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("ajmolstad/SpPDCC")
```
### Citation instructions
Please cite the most recent version of the article mentioned above. As of October 2023, this was the following (in bibtex): 
```
@article{molstad2024direct,
  title={Direct covariance matrix estimation with compositional data},
  author={Molstad, Aaron J and Ekvall, Karl Oskar and Suder, Piotr M},
  journal={Electronic Journal of Statistics},
  volume={18},
  number={1},
  pages={1702--1748},
  year={2024},
  publisher={The Institute of Mathematical Statistics and the Bernoulli Society}
}
```
### Vignette
Please visit [this example page](https://ajmolstad.github.io/docs/SpPDCCExample.html) for details on implementation and usage. 


### Reproducing simulation study results
Code to reproduce simulation results from the article can be found at [this repository](https://github.com/ajmolstad/CompositionalCovariance).
