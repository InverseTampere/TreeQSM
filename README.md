# TreeQSM

**Version 2.30**
**Reconstruction of quantitative structure models of trees from point cloud data**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.844626.svg)](https://doi.org/10.5281/zenodo.844626)

TreeQSM is a modelling method that reconstructs quantitative structure models (QSMs) of trees from point clouds. A QSM consists of a hierarchical collection of cylinders which estimate topological, geometrical and volumetric details of the woody structure of the tree. The input point cloud, which is usually produced by a terrestrial laser scanner, must contain only one tree, but the point cloud may contain also some points from the ground and understory. Much more details of the method and QSMs can be found from the manual that is part of the code distribution.

Web: http://math.tut.fi/inversegroup/
Some published papers about the method: 	
Raumonen et al. 2013, Remote Sensing
Calders et al. 2015, Methods in Ecology and Evolution
Raumonen et al. 2015, ISPRS Annals
Åkerblom et al. 2015, Remote Sensing

The TreeQSM is written in Matlab.
The main function is _treeqsm.m_, which takes in a point cloud and a structure array specifying the needed parameters. Refer to the manual or the help documentation of a particular function for further details.