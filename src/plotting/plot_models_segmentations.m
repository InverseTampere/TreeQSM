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


%% figure 1,  segmented point cloud, colors denote the branching order
plot_tree_structure(P,cover,segment,1,1,0,1)

%% figure 2,  cylinder model, colors denote the branching order
Sta = cylinder.start;
P = mat_vec_subtraction(P,Sta(1,:));
if nargin > 5
    trunk = mat_vec_subtraction(trunk,Sta(1,:));
    Vert = double(triangulation.vert);
    Vert = mat_vec_subtraction(Vert,Sta(1,:));
end
Sta = mat_vec_subtraction(Sta,Sta(1,:));
cylinder.start = Sta;
plot_cylinder_model(cylinder,2,10,1)

%% figure 3, segmented point cloud and cylinder model in the same figure
plot_tree_structure(P,cover,segment,3,3,0,1)
hold on
plot_cylinder_model(cylinder,3,10,0.7)
hold off

if nargin > 4 
    %% figure 4, triangulation model (bottom) and cylinder model (top) of the stem
    Facets = double(triangulation.facet);
    CylInd = triangulation.cylind;
    fvd = triangulation.fvd;
    if max(size(Vert)) > 5
        Bran = cylinder.branch;
        nc = size(Bran,1);
        ind = (1:1:nc)';
        C = ind(Bran == 1);
        %fvd = ones(size(Facets,1),1);
        n = size(trunk,1);
        I = logical(round(0.55*rand(n,1)));
        figure(4)
        point_cloud_plotting(trunk(I,:),4,3)
        patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
        alpha(1)
        hold on
        plot_cylinder_model(cylinder,4,20,1,(CylInd:C(end)))
        axis equal
        hold off
    else
        disp('No triangulation model generated!')
    end
end