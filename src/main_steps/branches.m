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

function branch = branches(cylinder)

% ---------------------------------------------------------------------
% BRANCHES.M        Determines the branching structure and computes branch
%                       attributes
%
% Version 3.0.0
% Latest update     2 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------

% Determines the branches (cylinders in a segment define a branch), their order
% and topological parent-child-relation. Branch number one is the trunk and
% its order is zero. Notice that branch number does not tell its age in the
% sense that branch number two would be the oldest branch and the number
% three the second oldest.
%
% Inputs:
% cylinder  Cylinders, structure array
%
% Outputs:
% branch    Branch structure array, contains fields:
%             Branch order, parent, volume, length, angle, height, azimuth
%             and diameter
% ---------------------------------------------------------------------

% Changes from version 2.1.0 to 3.0.0, 2 May 2022:
% 1) Changed the code such that the input "segment" and output "cylinder"
%    are not needed anymore, which simplified the code in many places.
%    Cylinder info is now computed in "cylinders" function.

% Changes from version 2.0.0 to 2.1.0, 25 Jan 2020:
% 1) Chanced the coding to simplify and shorten the code
% 2) Added branch area and zenith direction as new fields in the
%    branch-structure array
% 3) Removed the line were 'ChildCyls' and'CylsInSegment' fields are
%    removed from the cylinder-structure array

Rad = cylinder.radius;
Len = cylinder.length;
Axe = cylinder.axis;

%% Branches
nc = size(Rad,1);  % number of cylinder
ns = max(cylinder.branch);  % number of segments
BData = zeros(ns,9); % branch ord, dia, vol, are, len, ang, hei, azi, zen
ind = (1:1:nc)';
CiB = cell(ns,1);
for i = 1:ns
  C = ind(cylinder.branch == i);
  CiB{i} = C;
  if ~isempty(C)

    BData(i,1) = cylinder.BranchOrder(C(1)); % branch order
    BData(i,2) = 2*Rad(C(1)); % branch diameter
    BData(i,3) = 1000*pi*sum(Len(C).*Rad(C).^2); % branch volume
    BData(i,4) = 2*pi*sum(Len(C).*Rad(C)); % branch area
    BData(i,5) = sum(Len(C)); % branch length

    % if the first cylinder is added to fill a gap, then
    % use the second cylinder to compute the angle:
    if cylinder.added(C(1)) && length(C) > 1
      FC = C(2);  % first cyl in the branch
      PC = cylinder.parent(C(1)); % parent cylinder of the branch
    else
      FC = C(1);
      PC = cylinder.parent(FC);
    end
    if PC > 0
      BData(i,6) = 180/pi*acos(Axe(FC,:)*Axe(PC,:)'); % branch angle
    end

    BData(i,7) = cylinder.start(C(1),3)-cylinder.start(1,3); % branch height
    BData(i,8) = 180/pi*atan2(Axe(C(1),2),Axe(C(1),1)); % branch azimuth
    BData(i,9) = 180/pi*acos(Axe(C(1),3)); % branch zenith
  end
end
BData = single(BData);

%% Branching structure (topology, parent-child-relation)
branch.order = uint8(BData(:,1));
BPar = zeros(ns,1);
Chi = cell(nc,1);
for i = 1:nc
  c = ind(cylinder.parent == i);
  c = c(c ~= cylinder.extension(i));
  Chi{i} = c;
end
for i = 1:ns
  C = CiB{i};
  ChildCyls = unique(vertcat(Chi{C}));
  CB = unique(cylinder.branch(ChildCyls));  % Child branches
  BPar(CB) = i;
end
if ns <= 2^16
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
