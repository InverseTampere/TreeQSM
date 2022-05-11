
% Creates input parameter structure array needed to run "treeqsm" function
% and "filtering" function.
% NOTE: use this code to define all the parameters but the PatchDiam and
% BallRad parameters can be conveniently defined by "define_input"
% function.
% 
% Last update 11 May 2022

clear inputs

%% QSM reconstruction parameters
%%% THE THREE INPUT PARAMETERS TO BE OPTIMIZED.
% These CAN BE VARIED AND SHOULD BE OPTIMIZED 
% One possibility to define these is to use "define_input" code
% (These can have multiple values given as vectors, e.g. [0.01 0.02]).
% Patch size of the first uniform-size cover:
inputs.PatchDiam1 = [0.08 0.12]; 
% Minimum patch size of the cover sets in the second cover:
inputs.PatchDiam2Min = [0.02 0.03]; 
% Maximum cover set size in the stem's base in the second cover:
inputs.PatchDiam2Max = [0.07 0.1]; 

%%% ADDITIONAL PATCH GENERATION PARAMETERS.
% The following parameters CAN BE VARIED BUT CAN BE USUALLY KEPT AS SHOWN
% (i.e. little bigger than PatchDiam parameters).
% One possibility to define these is to use "define_input" code
% Ball radius in the first uniform-size cover generation:
inputs.BallRad1 = inputs.PatchDiam1+0.015; 
% Maximum ball radius in the second cover generation:
inputs.BallRad2 = inputs.PatchDiam2Max+0.01; 

% The following parameters CAN BE USUALLY KEPT FIXED as shown.
% Minimum number of points in BallRad1-balls, generally good value is 3:
inputs.nmin1 = 3; 
% Minimum number of points in BallRad2-balls, generally good value is 1:
inputs.nmin2 = 1; 
% Does the point cloud contain points only from the tree (if 1, then yes):
inputs.OnlyTree = 1; 
% Produce a triangulation of the stem's bottom part up to the first main
% branch (if 1, then yes):
inputs.Tria = 0; 
% Compute the point-model distances (if 1, then yes):
inputs.Dist = 1; 

%%% RADIUS CORRECTION OPTIONS FOR MODIFYING TOO LARGE AND TOO SMALL CYLINDERS.
% These parameters CAN BE USUALLY KEPT FIXED as shown.
% Traditional TreeQSM choices:
% Minimum cylinder radius, used particularly in the taper corrections:
inputs.MinCylRad = 0.0025; 
% Radius correction based on radius of the parent. If 1, radii in a branch 
% are always smaller than the radius of the parent in the parent branch:
inputs.ParentCor = 1; 
% Parabola taper correction of radii inside branches. If 1, use the
% correction:
inputs.TaperCor = 1; 

% Growth volume correction approach introduced by Jan Hackenberg, 
% allometry: Radius = a*GrowthVol^b+c
inputs.GrowthVolCor = 0; % If 1, use growth volume (GV) correction 
% fac-parameter of the GV-approach, defines upper and lower bound. When 
% using GV-approach, consider setting TaperCorr = 0, ParentCorr = 0, 
% MinCylinderRadius = 0.
inputs.GrowthVolFac = 1.5; % Defines the allowed radius: 
%  1/fac*predicted_radius <= radius <= fac*predicted_radius
%  However, the radii of the branch tip cylinders are not increased.

%% Filtering parameters
% NOTE: These are all optional, but needed to run the "filtering" function.
% Statistical k-nearest neighbor distance outlier filtering, applied if 
% filter.k > 0. The value filter.k is the number of nearest neighbors. 
inputs.filter.k = 10;
% Statistical point density outlier filtering, applied if filter.radius > 0. 
% The value filter.radius is the radius of the ball neighborhood. This is
% usually meant as alternative to the above knn-filtering.
inputs.filter.radius = 0.00;
% The value filter.nsigma is the multiplier of the standard deviation of
% the kth-nearest neighbor distance/point density and points whose 
% kth-nearest neighbor distance/point density is larger/lower than the 
% average +/- filter.nsigma * std are removed:
inputs.filter.nsigma = 1.5;
% Small component filtering is applied if filter.ncomp > 0. This filter is
% based on cover whose patches are defined by filter.PatchDiam1 and 
% filter.BallRad1. The points which are included in components that have
% less than filter.ncomp patches are removed:
inputs.filter.PatchDiam1 = 0.05;
inputs.filter.BallRad1 = 0.075;
inputs.filter.ncomp = 2;
% Cubical downsampling is applied if filter.EdgeLength > 0. 
% The value filter.EdgeLength is the length of the cube edges:
inputs.filter.EdgeLength = 0.004;
% Plot the filtering results automatically after the filtering if
% filter.plot > 0
inputs.filter.plot = 1;

%% Other inputs
% These parameters don't affect the QSM-reconstruction but define what is
% saved, plotted, and displayed and how the models are named/indexed
% Name string for saving output files and naming models:
inputs.name = 'tree'; 
% Tree index. If modelling multiple trees, then they can be indexed uniquely:
inputs.tree = 1;
% Model index, can separate models if multiple models with the same inputs:
inputs.model = 1; 
% Save the output struct QSM as a matlab-file into \result folder. 
% If name = 'pine', tree = 2, model = 5, the name of the saved file is 
% 'QSM_pine_t2_m5.mat':
inputs.savemat = 1; 
% Save the models in .txt-files (check "save_model_text.m"):
inputs.savetxt = 1; 
% What are plotted during reconstruction process: 
% 2 = plots the QSM, the segmentated point cloud and distributions, 
% 1 = plots the QSM and the segmentated point cloud
% 0 = plots nothing
inputs.plot = 2; 
% What are displayed during the reconstruction: 2 = display all; 
% 1 = display name, parameters and distances; 0 = display only the name:
inputs.disp = 2; 
