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

function [treedata,triangulation] = tree_data(cylinder,branch,inputs,trunk,Disp)

% ---------------------------------------------------------------------
% TREE_DATA.M       Calculates some tree attributes from cylinder QSM. 
%
% Version 2.0.2
% Latest update     26 Nov 2019
%
% Copyright (C) 2013-2019 Pasi Raumonen
% ---------------------------------------------------------------------

% Inputs:
% cylinder:
%   radius (Rad)    Radii of the cylinders
%   length (Len)    Lengths of the cylinders
%   start (Sta)     Starting points of the cylinders
%   axis (Axe)      Axes of the cylinders
% branch:
%   order (BOrd)    Branch order data
%   volume (BVol)   Branch volume data
%   length (BLen)   Branch length data
% trunk     Point cloud of the trunk
% Disp      If one ot true, prints the tree data
%
% Output:
% treedata     Tree data/attributes in a struct
% ---------------------------------------------------------------------

% Changes from version 2.0.1 to 2.0.2, 26 Nov 2019:  
% 1) Bug fix: Added a statement "C < nc" for a while command that makes sure 
%    that the index "C" does not exceed the number of stem cylinders, when 
%    determining the index of cylinders up to first branch.
% 2) Bug fix: Changed "for i = 1:BO" to "for i = 1:max(1,BO)" where 
%    computing branch order data.
% 3) Added the plotting of the triangulation model

% Changes from version 2.0.0 to 2.0.1, 9 Oct 2019:  
% 1) Bug fix: Changed the units (from 100m to 1m) for computing the branch 
%    length distribution: branch length per branch order.



% Define variables from CylData and BranchData
Rad = cylinder.radius;
Len = cylinder.length;
Sta = cylinder.start;
Axe = cylinder.axis;
CBra = cylinder.branch;
BOrd = branch.order;
BVol = branch.volume;
BLen = branch.length;

clear treedata
%% Tree attributes from cylinders
% Trunk cylinders
Trunk = CBra == 1;

% Volumes, lengths, branches and area
treedata.TotalVolume = 1000*pi*Rad.^2'*Len;
treedata.TrunkVolume = 1000*pi*Rad(Trunk).^2'*Len(Trunk);
treedata.BranchVolume = 1000*pi*Rad(~Trunk).^2'*Len(~Trunk);
bottom = min(Sta(:,3));
[top,i] = max(Sta(:,3));
if Axe(i,3) > 0
    top = top+Len(i)*Axe(i,3);
end
treedata.TreeHeight = top-bottom;
treedata.TrunkLength = sum(Len(Trunk));
treedata.BranchLength = sum(Len(~Trunk));
NB = length(BOrd)-1; % number of branches
treedata.NumberBranches = NB;
BO = max(BOrd); % maximum branch order
treedata.MaxBranchOrder = BO;
treedata.TotalArea = 2*pi*sum(Rad.*Len);

%% Diameter at breast height (dbh)
% Dbh from the QSM
i = 1;
n = nnz(Trunk);
nc = length(Rad);
ind = (1:1:nc)';
T = ind(Trunk);
while i < n && sum(Len(T(1:i))) < 1.3
    i = i+1;
end
DBHqsm = 2*Rad(T(i));
treedata.DBHqsm = DBHqsm;

% Determine DBH from cylinder fitted particularly to the correct place
% Select the trunk point set
A = mat_vec_subtraction(trunk,Sta(1,:));
h = A*Axe(1,:)';
I = h < 1.5;
J = h > 1.1;
I = I&J;
if nnz(I) > 100
    T = trunk(I,:);
    % Fit cylinder
    [R,~,~,Axis,~,conv,rel] = least_squares_cylinder(T,mean(T),Axe(i,:),DBHqsm);
    if 2*R > 0.8*DBHqsm && 2*R < 1.2*DBHqsm && abs(Axe(i,:)*Axis') > 0.9 && conv && rel
        treedata.DBHcyl = 2*R;
    else
        treedata.DBHcyl = DBHqsm;
    end
else
    treedata.DBHcyl = DBHqsm;
end

%% Trunk volume and DBH from triangulation
% Determine suitable cylinders up to first branch
C = 1;
while C < nc && cylinder.branch(C) == 1% && Rad(C) > 0.75*Rad(1)
    C = C+1;
end

n = nnz(Trunk);
i = 2;
while i < n && Sta(i,3) < Sta(C,3)
    i = i+1;
end
CylInd = i;
TrunkLenTri = Sta(CylInd,3)-Sta(1,3);

if inputs.Tria
    EmptyTriangulation = false;
    % Calculate the volumes
    if size(trunk,1) > 1000 && TrunkLenTri >= 1
        
        % Set the parameters for triangulation
        PointDensity = zeros(CylInd-1,1);
        for i = 1:CylInd-1
            I = trunk(:,3) >= Sta(i,3) & trunk(:,3) < Sta(i+1,3);
            PointDensity(i) = pi*Rad(i)*Len(i)/nnz(I);
        end
        d = max(PointDensity);
        if DBHqsm > 10
            MinTriaHeight = 0.05;
        else
            MinTriaHeight = 0.02;
        end
        TriaHeight = max(MinTriaHeight,3*sqrt(d));
        TriaWidth = TriaHeight;
        
        % Select the trunk point set used for triangulation
        I = trunk(:,3) <= Sta(CylInd,3);
        Stem = trunk(I,:);
        
        % Do the triangulation
        clear triangulation
        triangulation = zeros(1,0);
        j = 0;
        while isempty(triangulation) && j < 10
            triangulation = curve_based_triangulation(Stem,TriaHeight,TriaWidth);
            j = j+1;
            if isempty(triangulation) && j < 10
                disp('  Try triangulation again')
            end
        end
        
        if ~isempty(triangulation)
            triangulation.cylind = CylInd;
            % Dbh from triangulation
            Vert = triangulation.vert;
            h = Vert(:,3)-triangulation.bottom;
            [~,I] = min(abs(h-1.3));
            H = h(I);
            I = abs(h-H) < triangulation.triah/2;
            V = Vert(I,:);
            V = V([2:end 1],:)-V(1:end,:);
            d = sqrt(sum(V.*V,2));
            treedata.DBHtri = sum(d)/pi;
            % volumes from the triangulation
            treedata.TriaTrunkVolume = triangulation.volume;
            TrunkVolMix = treedata.TrunkVolume-...
                1000*pi*sum(Rad(1:CylInd-1).^2.*Len(1:CylInd-1))+triangulation.volume;
            treedata.MixTrunkVolume = TrunkVolMix;
            treedata.MixTotalVolume = TrunkVolMix+treedata.BranchVolume;
            treedata.TriaTrunkLength = TrunkLenTri;
            
            if Disp
                figure(5)
                Vert = triangulation.vert;
                Tria = triangulation.facet;
                fvd = triangulation.fvd;
                plot3(Vert(1,1),Vert(1,2),Vert(1,3))
                point_cloud_plotting(trunk,5,6)
                patch('Vertices',Vert,'Faces',Tria,'FaceVertexCData',fvd,'FaceColor','flat')
                axis equal
                alpha(0.8)
                pause(0.01)
            end
        else
            treedata.DBHtri = DBHqsm;
            treedata.TriaTrunkVolume = treedata.TrunkVolume;
            treedata.MixTrunkVolume = treedata.TrunkVolume;
            treedata.MixTotalVolume = treedata.TotalVolume;
            treedata.TriaTrunkLength = 0;
            EmptyTriangulation = true;
        end
    else
        treedata.DBHtri = DBHqsm;
        treedata.TriaTrunkVolume = treedata.TrunkVolume;
        treedata.MixTrunkVolume = treedata.TrunkVolume;
        treedata.MixTotalVolume = treedata.TotalVolume;
        treedata.TriaTrunkLength = 0;
        EmptyTriangulation = true;
    end
    
    if EmptyTriangulation
        disp('  No triangulation model produced')
        clear triangulation
        triangulation.vert = zeros(0,3);
        triangulation.facet = zeros(0,3);
        triangulation.fvd = zeros(0,1);
        triangulation.volume = 0;
        triangulation.bottom = 0;
        triangulation.top = 0;
        triangulation.triah = 0;
        triangulation.triaw = 0;
        triangulation.cylind = 0;
    end
else
    triangulation = 0;
end

%% Tree Location
treedata.location = Sta(1,:);

%% Stem taper
R = Rad(Trunk);
L = Len(Trunk);
n = length(R);
Taper = zeros(n+1,2);
Taper(1,2) = 2*R(1);
Taper(2:end,1) = cumsum(L);
Taper(2:end,2) = [2*R(2:end); 2*R(n)];
treedata.StemTaper = Taper';

%% Wood part size distributions
% Volume and length of wood parts as functions of the cylinder diameter
% (in 1cm diameter classes)
MaxDiam = ceil(max(200*Rad));
PartSizeDistri = zeros(MaxDiam,2);
for i = 1:MaxDiam
    if i > 1
        J = Rad <= i*0.005;
        K = I&J;
        I = ~J;
    else
        K = Rad <= i*0.005;
        I = ~K;
    end
    PartSizeDistri(i,1) = 1000*pi*sum(Rad(K).^2.*Len(K)); % volumes in litres
    PartSizeDistri(i,2) = sum(Len(K)); % lengths in meters
end
treedata.VolumeCylDiam = PartSizeDistri(:,1)';
treedata.LengthCylDiam = PartSizeDistri(:,2)';

%% Branch order data
BranchOrdDistri = zeros(BO,3);
for i = 1:max(1,BO)
    I = BOrd == i;
    BranchOrdDistri(i,1) = sum(BVol(I)); % volumes
    BranchOrdDistri(i,2) = sum(BLen(I)); % lengths
    BranchOrdDistri(i,3) = nnz(I); % number of ith-order branches
end
treedata.VolumeBranchOrder = BranchOrdDistri(:,1)';
treedata.LengthBranchOrder = BranchOrdDistri(:,2)';
treedata.NumberBranchOrder = BranchOrdDistri(:,3)';

if Disp
figure(6)
subplot(2,3,1)
plot(Taper(:,1),Taper(:,2),'-b')
title('Stem taper')
xlabel('Distance from base (m)')
ylabel('Diameter (m)')
grid on

subplot(2,3,2)
bar(1:1:MaxDiam,PartSizeDistri(:,1))
title('Tree segment volume per diameter class')
xlabel('Cylinder diameter class (cm)')
ylabel('Volume (L)')
axis tight
grid on

subplot(2,3,3)
bar(1:1:MaxDiam,PartSizeDistri(:,2))
title('Tree segment length per diameter class')
xlabel('Cylinder diameter class (cm)')
ylabel('Length (m)')
axis tight
grid on

subplot(2,3,4)
bar(1:BO,BranchOrdDistri(:,1))
title('Branch volume per branch order')
xlabel('Branch order')
ylabel('Volume (L)')
axis tight
grid on

subplot(2,3,5)
bar(1:BO,BranchOrdDistri(:,2))
title('Branch length per branch order')
xlabel('Branch order')
ylabel('Length (m)')
axis tight
grid on

subplot(2,3,6)
bar(1:BO,BranchOrdDistri(:,3))
title('Number of branch per branch order')
xlabel('Branch order')
ylabel('Number of branches')
axis tight
grid on
end

% change into single-format
Names = fieldnames(treedata);
n = size(Names,1);
for i = 1:n
    treedata.(Names{i}) = single(treedata.(Names{i}));
end
Units = zeros(n,3);
for i = 1:n
    if strcmp(Names{i},'DBHcyl') && ~inputs.Tria
        m = i;
    elseif strcmp(Names{i},'DBHcyl') && inputs.Tria
        m = i+5;
    end
    if strcmp(Names{i}(1:3),'DBH')
        Units(i,:) = 'm  ';
    elseif strcmp(Names{i}(end-2:end),'ume')
        Units(i,:) = 'L  ';
    elseif strcmp(Names{i}(end-2:end),'ght')
        Units(i,:) = 'm  ';
    elseif strcmp(Names{i}(end-2:end),'gth')
        Units(i,:) = 'm  ';
    elseif strcmp(Names{i}(1:3),'vol')
        Units(i,:) = 'L  ';
    elseif strcmp(Names{i}(1:3),'len')
        Units(i,:) = 'm  ';
    elseif strcmp(Names{i}(end-2:end),'rea')
        Units(i,:) = 'm^2';
    elseif strcmp(Names{i}(1:3),'loc')
        Units(i,:) = 'm  ';
    end
end

if Disp
    %% Display treedata
    disp('------------')
    disp('  Tree attributes:')
    for i = 1:m
        v = change_precision(treedata.(Names{i}));
        disp(['  ',Names{i},' = ',num2str(v),' ',Units(i,:)])
    end
    disp('  -----')
end
