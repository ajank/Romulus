### Romulus: Robust multi-state identification of transcription factor binding sites from DNase-seq data

Romulus is a computational method to accurately identify individual
transcription factor binding sites from genome sequence information and
cell-type--specific experimental data, such as DNase-seq. It combines the
strengths of its predecessors, CENTIPEDE and Wellington, while keeping
the number of free parameters in the model robustly low. The method is
unique in allowing for multiple binding states for a single transcription
factor, differing in their cut profile and overall number of DNase I cuts.

#### Documentation

A PDF version of the Romulus package manual is available from
http://www.mimuw.edu.pl/~ajank/Romulus/Romulus-manual.pdf.

#### Installation

The most recent stable version of Romulus is available from
http://www.mimuw.edu.pl/~ajank/Romulus/Romulus_1.0.2.tar.gz.
Please install it using the following R command:

```R
install.packages("Romulus_1.0.2.tar.gz", repos = NULL, type = "source")
```

Alternatively, the development version can be installed directly from GitHub
with the following two R commands:

```R
library(devtools)
install_github("ajank/Romulus")
```

If you do not have the package `devtools´ installed, please install it first:

```R
install.packages("devtools")
```

#### Citation

At the moment, please cite Romulus as:

Jankowski, A., Tiuryn, J. and Prabhakar, S. (2016) Romulus: Robust
multi-state identification of transcription factor binding sites
from DNase-seq data. Bioinformatics. doi: 10.1093/bioinformatics/btw209
