% This file is part of TREEQSM.
% 
% TREEQSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TREEQSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.


% creates input parameter structure array needed to run TreeQSM function

clear inputs

%% QSM reconstruction parameters
% The following parameters can be varied and should be optimised (they 
% can have multiple values given as vectors, e.g. [4 6]):
inputs.PatchDiam1 = [0.1 0.15]; % Patch size of the first uniform-size cover
inputs.PatchDiam2Min = [0.02 0.03]; % Minimum patch size of the cover sets in the second cover
inputs.PatchDiam2Max = [0.06 0.08]; % Maximum cover set size in the stem's base in the second cover
inputs.lcyl = [3 5]; % Relative (length/radius) length of the cylinders
inputs.FilRad = 3; % Relative radius for outlier point filtering

% The following parameters can be varied and but usually can be kept as
% shown (i.e. little bigger than PatchDiam parameters):
inputs.BallRad1 = inputs.PatchDiam1+0.02; % Ball radius in the first uniform-size cover generation
inputs.BallRad2 = inputs.PatchDiam2Max+0.01; % Maximum ball radius in the second cover generation

% The following parameters can be usually kept fixed as shown:
inputs.nmin1 = 3; % Minimum number of points in BallRad1-balls, generally good value is 3
inputs.nmin2 = 1; % Minimum number of points in BallRad2-balls, generally good value is 1
inputs.OnlyTree = 1; % If 1, point cloud contains points only from the tree
inputs.Tria = 1; % If 1, produces a triangulation 
inputs.Dist = 1; % If 1, computes the point-model distances

% Different cylinder radius correction options for modifying too large and
% too small cylinders:
% Traditional TreeQSM choices:
inputs.MinCylRad = 0.0025; % Minimum cylinder radius, used particularly in the taper corrections
inputs.ParentCor = 1; % Radii in a child branch are always smaller than the radii of the 
                      % parent cylinder in the parent branch
inputs.TaperCor = 1; % Use partially linear (stem) and parabola (branches) taper corrections
% Growth volume correction approach introduced by Jan Hackenberg, 
% allometry: GrowthVol = a*Radius^b+c
inputs.GrowthVolCor = 0; % Use growth volume (GV) correction 
inputs.GrowthVolFac = 2.5; % fac-parameter of the GV-approach, defines upper and lower bound
                           % When using GV-approach, consider setting:
                           % TaperCorr = 0, ParentCorr = 0, MinCylinderRadius = 0.

%% Other inputs
% These parameters don't affect the QSM-reconstruction but define what is
% saved, plotted, and displayed and how the models are named/indexed
inputs.name = 'pine'; % Name string for saving output files and naming models
inputs.tree = 1; % Tree index. If modelling multiple trees, then they can be indexed uniquely
inputs.model = 1; % Model index, can separate models if multiple models with the same inputs
inputs.savemat = 1; % If 1, saves the output struct QSM as a matlab-file into \result folder 
                    % If name = 'pine', tree = 2, model = 5, 
                    % the name of the saved file is 'QSM_pine_t2_m5.mat'
inputs.savetxt = 1; % If 1, saves the models in .txt-files
inputs.plot = 1; % If 1, plots the model, the segmentation of the point cloud and distributions
inputs.disp = 2; % Defines what is displayed during the reconstruction:
                 % 2 = display all; 1 = display name, parameters and distances;
                 % 0 = display only the name
