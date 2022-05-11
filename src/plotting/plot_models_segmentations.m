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

function plot_models_segmentations(P,cover,segment,cylinder,trunk,triangulation)

% ---------------------------------------------------------------------
% PLOT_MODELS_SEGMENTATION.M   Plots the segmented point clouds and
%                               cylinder/triangulation models
%
% Version 1.1.0
% Latest update     13 July 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------

% Inputs:
% P             Point cloud
% cover         cover-structure array
% segment       segment-structure array
% cylinder      cylinder-structure array
% trunk         point cloud of the trunk
% triangulation triangulation-structure array

% Changes from version 1.0.0 to 1.1.0, 13 July 2020:
% 1) plots now figure 1 and 2 with two subplots; in the first the colors 
%    are based on branching order and in the second they are based on
%    branch

%% figure 1: branch-segmented point cloud 
% colors denote the branching order and branches
figure(1)
subplot(1,2,1)
plot_branch_segmentation(P,cover,segment,'order')
subplot(1,2,2)
plot_branch_segmentation(P,cover,segment,'branch')

%% figure 2: cylinder model 
% colors denote the branching order and branches
Sta = cylinder.start;
P = P-Sta(1,:);
if nargin > 5
    trunk = trunk-Sta(1,:);
    Vert = double(triangulation.vert);
    Vert = Vert-Sta(1,:);
end
Sta = Sta-Sta(1,:);
cylinder.start = Sta;
figure(2)
subplot(1,2,1)
plot_cylinder_model(cylinder,'order',2,10)
subplot(1,2,2)
plot_cylinder_model(cylinder,'branch',2,10)

%% figure 3, segmented point cloud and cylinder model
plot_branch_segmentation(P,cover,segment,'order',3,1)
hold on
plot_cylinder_model(cylinder,'order',3,10,0.7)
hold off

if nargin > 4 
    %% figure 4, triangulation model (bottom) and cylinder model (top) 
    % of the stem
    Facets = double(triangulation.facet);
    CylInd = triangulation.cylind;
    fvd = triangulation.fvd;
    if max(size(Vert)) > 5
        Bran = cylinder.branch;
        nc = size(Bran,1);
        ind = (1:1:nc)';
        C = ind(Bran == 1);
        n = size(trunk,1);
        I = logical(round(0.55*rand(n,1)));
        figure(4)
        point_cloud_plotting(trunk(I,:),4,3)
        patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,...
            'FaceColor','flat')
        alpha(1)
        hold on
        plot_cylinder_model(cylinder,'order',4,20,1,(CylInd:C(end)))
        axis equal
        hold off
    else
        disp('No triangulation model generated!')
    end
end