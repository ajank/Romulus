MOCCA is a novel computational method to accurately identify TF
footprints from genome sequence information and cell-type--specific
experimental data, such as DNase-seq data. Our approach combines the
strengths of CENTIPEDE and Wellington, while keeping the number of free
parameters in the model reasonably low. For a given TF, we first
identify candidate binding sites that have reasonable sequence affinity,
using a position weight matrix. Then, like CENTIPEDE, we employ an
Expectation-Maximization-based approach to simultaneously learn the
DNase I cut profiles and classify the binding sites as bound or unbound.

Our method is unique in allowing for multiple bound states for a single
TF, differing in their cut profile and overall number of DNase I cuts.
To make the model robust, we employ a systematic approach to group the
DNase I cuts, according to their location and strand. Inspired by
Wellington, we take the forward strand DNase I cuts only upstream and
within the cut site, while the reverse strand DNase I cuts -- within the
cut site and downstream. We model the total number of cuts as a negative
binomial component, while the cut distribution (regularized by binning
outside the cut site) is modeled as a multinomial component. Overall,
MOCCA predictions agree well with experimental ChIP-seq measurements of
TF binding at candidate motif sites.
