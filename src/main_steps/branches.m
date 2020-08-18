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
% Version 2.1.0
% Latest update     25 Jan 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
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
% ---------------------------------------------------------------------

% Changes from version 2.0.0 to 2.1.0, 25 Jan 2020:
% 1) Changed the coding to simplify and shorten the code
% 2) Added branch area and zenith direction as new fields in the 
%    branch-structure array
% 3) Removed the line were 'ChildCyls' and'CylsInSegment' fields are
%    removed from the cylinder-structure array 

Rad = cylinder.radius;
Len = cylinder.length;
Axe = cylinder.axis;

%% Branches
nc = size(Rad,1);  % number of cylinder
ns = size(segment.segments,1);  % number of segments
CiB = cell(ns,1); % Cylinders in each branch
C = cylinder.CylsInSegment{1};
CiB{1} = C;
CData = zeros(nc,3); % cylinder's branch, b-order, position in the branch
CData(C,2) = 1;
CData(C,3) = (1:length(C))';
BData = zeros(ns,9); % branch ord, dia, vol, are, len, ang, hei, azi, zen
BData(1,2) = 2*Rad(1); % diameter
BData(1,3) = 1000*pi*sum(Len(C).*Rad(C).^2); % volume
BData(1,4) = 2*pi*sum(Len(C).*Rad(C)); % area
BData(1,5) = sum(Len(C)); % length
BData(1,8) = 180/pi*atan2(Axe(1,2),Axe(1,1)); % azimuth
BData(1,9) = 180/pi*acos(Axe(1,3)); % zenith
S = segment.ChildSegment{1};    % segments under inspection
b = 1;          % branches determined so far
BO = 0;         % branch order under inspection
while ~isempty(S)
    BO = BO+1;
    n = length(S);
    for j = 1:n
        C = cylinder.CylsInSegment{S(j)};
        if ~isempty(C)
            b = b+1;
            CiB{b} = C;
            CData(C,1) = BO;
            CData(C,2) = b;
            CData(C,3) = (1:length(C))';
            
            BData(b,1) = BO; % branch order
            BData(b,2) = 2*Rad(C(1)); % branch diameter
            BData(b,3) = 1000*pi*sum(Len(C).*Rad(C).^2); % branch volume
            BData(b,4) = 2*pi*sum(Len(C).*Rad(C)); % branch area
            BData(b,5) = sum(Len(C)); % branch length
            
            % if the first cylinder is added to fill a gap, then
            % use the second cylinder:
            added = cylinder.added(C(1));
            if added
                FC = C(2);  % first cyl in the branch
                PC = cylinder.parent(C(added)); % parent cylinder of the branch
            else
                FC = C(1);
                PC = cylinder.parent(FC);
            end
            if PC > 0
                BData(b,6) = 180/pi*acos(Axe(FC,:)*Axe(PC,:)'); % branch angle
            end
            
            BData(b,7) = cylinder.start(C(1),3)-cylinder.start(1,3); % branch height
            BData(b,8) = 180/pi*atan2(Axe(C(1),2),Axe(C(1),1)); % branch azimuth
            BData(b,9) = 180/pi*acos(Axe(C(1),3)); % branch zenith
        end
    end
    S = vertcat(segment.ChildSegment{S});
end
CiB = CiB(1:b);
BData = single(BData(1:b,:));
clear branch

%% Branching structure (topology, parent-child-relation)
branch.order = uint8(BData(:,1));
BPar = zeros(b,1);
for i = 1:b
    C = CiB{i};
    ChildCyls = unique(vertcat(cylinder.ChildCyls{C}));
    CB = unique(CData(ChildCyls,2));  % Child branches
    BPar(CB) = i; 
end
if b <= 2^16
    branch.parent = uint16(BPar);
else
    branch.parent = uint32(BPar);
end

%% Finish the definition of branch
branch.diameter = BData(:,2);   % diameters in meters
branch.volume = BData(:,3);     % volumes in liters
branch.area = BData(:,4);       % areas in square meters
branch.length = BData(:,5);     % lengths in meters
branch.angle = BData(:,6);      % angles in degrees
branch.height = BData(:,7);     % heights in meters
branch.azimuth = BData(:,8);    % azimuth directions in angles
branch.zenith = BData(:,9);     % zenith directions in angles

%% Add cylinder outputs
if b <= 2^16
    cylinder.branch = uint16(CData(:,2));
else
    cylinder.branch = uint32(CData(:,2));
end
cylinder.BranchOrder = uint8(CData(:,1));
if max(CData(:,3)) <= 2^8
    cylinder.PositionInBranch = uint8(CData(:,3));
else
    cylinder.PositionInBranch = uint16(CData(:,3));
end
