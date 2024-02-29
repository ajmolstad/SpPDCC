# SpPDCC
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
@misc{molstad2022direct,
      title={Direct covariance matrix estimation with compositional data}, 
      author={Aaron J. Molstad and Karl Oskar Ekvall and Piotr M. Suder},
      year={2022},
      eprint={2212.09833},
      archivePrefix={arXiv},
      primaryClass={stat.ME}
}
```
### Vignette
Please visit [this example page](https://ajmolstad.github.io/docs/SpPDCCExample.html) for details on implementation and usage. 


### Reproducing simulation study results
Code to reproduce simulation results from the article can be found at [this repository](https://github.com/ajmolstad/CompositionalCovariance).
