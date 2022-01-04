# ek-ptv
Particle tracking algorithm based on package "trackpy" [1]

Objective: To obtain individual particle velocity from a sequence of images. The velocity measures are then be used to calculate an electrokinetic (EK) mobility via statistical analysis (R code available upon request).

Quick description

1) Read image sequence file (TIFF format) using PIMS [2]
2) Detect location of particles in every frame of the sequence
3) Link relevant locations through multiple frames and create trace information for every particle detected
4) Filter out noises
5) Convert trace information into physical quantities (units in microns and seconds)
6) Save converted trace data for statistical analysis

(Details are described in [1])

References

[1] Allan, Daniel B., Caswell, Thomas, Keim, Nathan C., van der Wel, Casper M., & Verweij, Ruben W. (2021). soft-matter/trackpy: Trackpy v0.5.0 (v0.5.0). Zenodo. https://doi.org/10.5281/zenodo.4682814

[2] Allan, Daniel B., Caswell, Thomas, Keim, Nathan C., van der Wel, Casper M., & Dimiduk, Thomas (2020). PIMS: Python Image Sequence. http://soft-matter.github.io/pims/v0.5/index.html#
