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

function [branch,cylinder] = branches(segment,cylinder)

% ---------------------------------------------------------------------
% BRANCHES.M        Determines the branching structure and computes branch
%                       attributes
%
% Version 2.00
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Determines the branches (cylinders in a segment define a branch), their order
% and topological parent-child-relation. Branch number one is the trunk and
% its order is zero. Notice that branch number does not tell its age in the 
% sense that branch number two would be the oldest branch and the number 
% three the second oldest. 
%
% Inputs:
% segment   Segments, structure array
% cylinder  Cylinders, structure array
%
% Outputs:
% branch    Branch structure array, contains fields:
%             Branch order, parent, volume, length, angle, height, azimuth and diameter
% cylinder  Updated cylinder structure array, contains new fields:
%             Branch of the cylinder, branch order, position inside the branch

Segs = segment.segments;
SChi = segment.ChildSegment;
Rad = cylinder.radius;
Len = cylinder.length;
Sta = cylinder.start;
Axe = cylinder.axis;
CPar = cylinder.parent;
Added = cylinder.added;
CChi = cylinder.ChildCyls;
CiS = cylinder.CylsInSegment;
cylinder = rmfield(cylinder,{'ChildCyls','CylsInSegment'});

%% Branches
nc = size(Rad,1);  % number of cylinder
ns = size(Segs,1);  % number of segments
BOrd = zeros(ns,1); % Branch order of the branch
CiB = cell(ns,1); % Cylinders in each branch
BoC = zeros(nc,1); % Branch of cylinder
BOrdC = zeros(nc,1); % Branch order of the cylinder
CPiB = zeros(nc,1); % Cylinder's position inside the branch
C = CiS{1};
CiB{1} = C;
BoC(C) = 1;
CPiB(C) = (1:length(C))';
BVol = zeros(ns,1); % branch volume
BLen = zeros(ns,1); % branch length
BAng = zeros(ns,1); % branch angle
BHei = zeros(ns,1); % branch height
BAzi = zeros(ns,1); % azimuth of branch
BDia = zeros(ns,1); % branch diameter at the base
BDia(1) = 2*Rad(1);
BVol(1) = pi*sum(Len(C).*Rad(C).^2);
BLen(1) = sum(Len(C));

S = SChi{1};    % segments under inspection
b = 1;          % branches determined so far
BO = 0;         % branch order under inspection
while ~isempty(S)
    BO = BO+1;
    n = length(S);
    for j = 1:n
        C = CiS{S(j)};
        if ~isempty(C)
            b = b+1;
            CiB{b} = C;
            BOrd(b) = BO;
            BoC(C) = b;
            BOrdC(C) = BO;
            CPiB(C) = (1:length(C))';
            BVol(b) = pi*sum(Len(C).*Rad(C).^2);
            BLen(b) = sum(Len(C));
            BHei(b) = Sta(C(1),3)-Sta(1,3);
            BAzi(b) = 180/pi*atan2(Axe(C(1),2),Axe(C(1),1));
            BDia(b) = 2*Rad(C(1));
            I = Added(C(1));  % if the first cylinder is added to fill gap, 
                              % use the second cylinder
            if I
                FC = C(2);          % first cyl in the branch
                PC = CPar(C(I));    % parent cylinder of the branch
            else
                FC = C(1);
                PC = CPar(FC);
            end
            if PC > 0
                BAng(b) = 180/pi*acos(Axe(FC,:)*Axe(PC,:)');
            end
        end
    end
    S = vertcat(SChi{S});
end
CiB = CiB(1:b);
clear branch

%% Branching structure (topology, parent-child-relation)
branch.order = uint8(BOrd(1:b));
BPar = zeros(b,1);
for i = 1:b
    C = CiB{i};
    ChildCyls = unique(vertcat(CChi{C}));
    CB = unique(BoC(ChildCyls));  % Child branches
    BPar(CB) = i; 
end
if b <= 2^16
    branch.parent = uint16(BPar);
else
    branch.parent = uint32(BPar);
end

%% Finish the definition of branch
branch.volume = single(1000*BVol(1:b));  % volumes in liters
branch.length = single(BLen(1:b));  % lengths in meters
branch.angle = single(BAng(1:b)); % angles in degrees
branch.height = single(BHei(1:b)); % heights in meters
branch.azimuth = single(BAzi(1:b)); % aximuths in angles
branch.diameter = single(BDia(1:b));  % diameters in meters

%% Add cylinder outputs
if b <= 2^16
    cylinder.branch = uint16(BoC);
else
    cylinder.branch = uint32(BoC);
end
cylinder.BranchOrder = uint8(BOrdC);
if max(CPiB) <= 2^8
    cylinder.PositionInBranch = uint8(CPiB);
else
    cylinder.PositionInBranch = uint16(CPiB);
end
