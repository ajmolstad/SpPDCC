# SpPDCC
An R package for estimating sparse and positive definite basis covariance matrices from compositional data using method described in [Direct covariance matrix estimation with compositional data](https://arxiv.org/abs/2212.09833). Please contact amolstad@umn.edu with any questions or comments. 

Note that the package will be sporadically updated. We hope that subsequent updates will improve computing times. 

### Installation
SpPDCC can be loaded directly into R through the the `devtools` package:
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
### Usage directions
Please visit [this example page](https://ajmolstad.github.io/docs/SpPDCCExample.html) for details on implementation and usage. 
