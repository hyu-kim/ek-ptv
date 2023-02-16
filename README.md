# ek-ptv
### Single-cell tracking algorithm based on a package "trackpy" <sup>1</sup>.

Now published in Frontiers in Bioengineering and Biotechnology.
> Wang Q, Kim H, Halvorsen TM, Chen S, Hayes CS and Buie CR (2023) Leveraging microfluidic dielectrophoresis to distinguish compositional variations of lipopolysaccharide in E. coli. Front. Bioeng. Biotechnol. 11:991784. doi: 10.3389/fbioe.2023.991784

Abstract
> Lipopolysaccharide (LPS) is the unique feature that composes the outer leaflet of the Gram-negative bacterial cell envelope. Variations in LPS structures affect a number of physiological processes, including outer membrane permeability, antimicrobial resistance, recognition by the host immune system, biofilm formation, and interbacterial competition. Rapid characterization of LPS properties is crucial for studying the relationship between these LPS structural changes and bacterial physiology. However, current assessments of LPS structures require LPS extraction and purification followed by cumbersome proteomic analysis. This paper demonstrates one of the first high-throughput and noninvasive strategies to directly distinguish Escherichia coli with different LPS structures. Using a combination of three-dimensional insulator-based dielectrophoresis (3DiDEP) and cell tracking in a linear electrokinetics assay, we elucidate the effects of structural changes in E. coli LPS oligosaccharides on electrokinetic mobility and polarizability. We show that our platform is sufficiently sensitive to detect the LPS structural variations at molecular levels. To correlate the electrokinetic properties of LPS with the outer membrane permeability, we further examined the effects of LPS structural variations on bacterial susceptibility to colistin, an antibiotic known to disrupt the outer membrane by targeting LPS. Our results suggest that microfluidic electrokinetic platforms employing 3DiDEP can be a useful tool for isolating and selecting bacteria based on their LPS glycoforms. Future iterations of these platforms could be leveraged for rapid profiling of pathogens based on their surface LPS structural identity.

### Description
1) Read image sequence file (TIFF format) using PIMS <sup>2</sup>
2) Detect location of particles in every frame of the sequence
3) Link relevant locations through multiple frames and create trace information for every particle detected
4) Filter out noises
5) Convert trace information into physical quantities (units in microns and seconds)
6) Save converted trace data for statistical analysis via R

### References
1) Allan, Daniel B., Caswell, Thomas, Keim, Nathan C., van der Wel, Casper M., & Verweij, Ruben W. (2021). soft-matter/trackpy: Trackpy v0.5.0 (v0.5.0). Zenodo. https://doi.org/10.5281/zenodo.4682814
2) Allan, Daniel B., Caswell, Thomas, Keim, Nathan C., van der Wel, Casper M., & Dimiduk, Thomas (2020). PIMS: Python Image Sequence. http://soft-matter.github.io/pims/v0.5/index.html#
