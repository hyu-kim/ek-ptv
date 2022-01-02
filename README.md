# ek-ptv
Particle tracking algorithm based on package "trackpy" [1]

Objective: To obtain individual particle velocity from a sequence of images. The velocity measures are then be used to calculate an electrokinetic (EK) mobility via statistical analysis (R code available upon request).

Quick description (details are described in [1])
1) Read image sequence file (TIFF format) using PIMS [2]
2) Detect location of particles in every frame of the sequence
3) Link relevant locations through multiple frames and create trace information for every particle detected
4) Filter out noises
5) Convert trace information into physical quantities (units in microns and seconds)
6) Save converted trace data for statistical analysis
