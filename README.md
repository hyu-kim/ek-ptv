# ek-ptv
### Single-cell tracking algorithm based on a package "trackpy".

- Objective: To obtain individual particle velocity from a sequence of images. The velocity measures are then be used to calculate an electrokinetic (EK) mobility via statistical analysis (R code available upon request).
- Now in public at bioRxiv..
> Wang, Q.*; Kim, H.*; Halvorsen, T. M.; Buie, C. R., Leveraging microfluidic dielectrophoresis to distinguish compositional variations of lipopolysaccharide in E. coli, bioRxiv 2022, 2022.08.19.504570

### Description

1) Read image sequence file (TIFF format) using PIMS [2]
2) Detect location of particles in every frame of the sequence
3) Link relevant locations through multiple frames and create trace information for every particle detected
4) Filter out noises
5) Convert trace information into physical quantities (units in microns and seconds)
6) Save converted trace data for statistical analysis via R

### References

1) Allan, Daniel B., Caswell, Thomas, Keim, Nathan C., van der Wel, Casper M., & Verweij, Ruben W. (2021). soft-matter/trackpy: Trackpy v0.5.0 (v0.5.0). Zenodo. https://doi.org/10.5281/zenodo.4682814
2) Allan, Daniel B., Caswell, Thomas, Keim, Nathan C., van der Wel, Casper M., & Dimiduk, Thomas (2020). PIMS: Python Image Sequence. http://soft-matter.github.io/pims/v0.5/index.html#
