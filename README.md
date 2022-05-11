# TreeQSM

**Version 2.4.1**
**Reconstruction of quantitative structure models for trees from point cloud data**

[![DOI](https://zenodo.org/badge/100592530.svg)](https://zenodo.org/badge/latestdoi/100592530)

![QSM image](https://github.com/InverseTampere/TreeQSM/blob/master/Manual/fig_point_cloud_qsm.png)


### Description

TreeQSM is a modelling method that reconstructs quantitative structure models (QSMs) for trees from point clouds. A QSM consists of a hierarchical collection of cylinders estimating topological, geometrical and volumetric details of the woody structure of the tree. The input point cloud, which is usually produced by a terrestrial laser scanner, must contain only one tree, which is intended to be modelled, but the point cloud may contain also some points from the ground and understory. Moreover, the point cloud should not contain significant amount of noise or points from leaves as these are interpreted as points from woody parts of the tree and can therefore lead to erroneous results. Much more details of the method and QSMs can be found from the manual that is part of the code distribution.

The TreeQSM is written in Matlab.
The main function is _treeqsm.m_, which takes in a point cloud and a structure array specifying the needed parameters. Refer to the manual or the help documentation of a particular function for further details.

### References

Web: https://research.tuni.fi/inverse/
Some published papers about the method and applications:  
Raumonen et al. 2013, Remote Sensing https://www.mdpi.com/2072-4292/5/2/491  
Calders et al. 2015, Methods in Ecology and Evolution https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12301  
Raumonen et al. 2015, ISPRS Annals https://www.isprs-ann-photogramm-remote-sens-spatial-inf-sci.net/II-3-W4/189/2015/  
Åkerblom et al. 2015, Remote Sensing https://www.mdpi.com/2072-4292/7/4/4581  
Åkerblom et al. 2017, Remote Sensing of Environment https://www.sciencedirect.com/science/article/abs/pii/S0034425716304746  
de Tanago Menaca et al. 2017, Methods in Ecology and Evolution https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12904  
Åkerblom et al. 2018, Interface Focus http://dx.doi.org/10.1098/rsfs.2017.0045   
Disney et al. 2018, Interface Focus http://dx.doi.org/10.1098/rsfs.2017.0048   


### Quick guide

Here is a quick guide for testing the code and starting its use. However, it is highly recommended that after the testing the user reads the manual for more information how to best use the code.  

1) Start MATLAB and set the main path to the root folder, where _treeqsm.m_ is located.\
2) Use _Set Path_ --> _Add with Subfolders_ --> _Open_ --> _Save_ --> _Close_ to add the subfolders, where all the codes of the software are, to the paths of MATLAB.\
3) Import a point cloud of a tree into the workspace. Let us name it P.\
4) Define suitable inputs:\
  &nbsp; &nbsp; >> inputs = define_input(P,1,1,1);\
5) Reconstruct QSMs:\
  &nbsp; &nbsp; >> QSM = treeqsm(P,inputs);  
