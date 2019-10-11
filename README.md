# TreeQSM

**Version 2.3.1**
**Reconstruction of quantitative structure models of trees from point cloud data**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.844626.svg)](https://doi.org/10.5281/zenodo.844626) (DOI for version 2.3.0)

TreeQSM is a modelling method that reconstructs quantitative structure models (QSMs) of trees from point clouds. A QSM consists of a hierarchical collection of cylinders which estimate topological, geometrical and volumetric details of the woody structure of the tree. The input point cloud, which is usually produced by a terrestrial laser scanner, must contain only one tree, but the point cloud may contain also some points from the ground and understory. Much more details of the method and QSMs can be found from the manual that is part of the code distribution.

Web: http://math.tut.fi/inversegroup/
Some published papers about the method and applications: 	
Raumonen et al. 2013, Remote Sensing https://www.mdpi.com/2072-4292/5/2/491
Calders et al. 2015, Methods in Ecology and Evolution https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12301
Raumonen et al. 2015, ISPRS Annals https://www.isprs-ann-photogramm-remote-sens-spatial-inf-sci.net/II-3-W4/189/2015/
Åkerblom et al. 2015, Remote Sensing https://www.mdpi.com/2072-4292/7/4/4581
Åkerblom et al. 2017, Remote Sensing of Environment https://www.sciencedirect.com/science/article/abs/pii/S0034425716304746
de Tanago Menaca et al. 2017, Methods in Ecology and Evolution https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12904
Åkerblom et al. 2018, Interface Focus http://dx.doi.org/10.1098/rsfs.2017.0045 
Disney et al. 2018, Interface Focus http://dx.doi.org/10.1098/rsfs.2017.0048 


The TreeQSM is written in Matlab.
The main function is _treeqsm.m_, which takes in a point cloud and a structure array specifying the needed parameters. Refer to the manual or the help documentation of a particular function for further details.