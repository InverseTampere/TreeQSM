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

function [Partition,CubeCoord,Info,Cubes] = cubical_partition(P,EL,NE)

% ---------------------------------------------------------------------
% CUBICAL_PARTITION.M    Partitions the point cloud into cubes.
%
% Version 1.1.0
% Latest update     6 Oct 2021
%
% Copyright (C) 2015-2021 Pasi Raumonen
% ---------------------------------------------------------------------

% Inputs:
% P           Point cloud, (n_points x 3)-matrix
% EL          Length of the cube edges
% NE          Number of empty edge layers
%
% Outputs:
% Partition   Point cloud partitioned into cubical cells,
%                 (nx x ny x nz)-cell, where nx,ny,nz are the number
%                 of cubes in x,y,z-directions, respectively. If "Cubes"
%                 is outputed, then "Partition" is (n x 1)-cell, where each
%                 cell corresponds to a nonempty cube.
%
% CC          (n_points x 3)-matrix whose rows are the cube coordinates
%                 of each point: x,y,z-coordinates
% Info        The minimum coordinate values and number of cubes in each
%                 coordinate direction
% Cubes       (Optional) (nx x ny x nz)-matrix (array), each nonzero
%                 element indicates that its cube is nonempty and the
%                 number indicates which cell in "Partition" contains the
%                 points of the cube.
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 1.1.0, 6 Oct 2021:
% 1) Changed the determinationa EL and NE so that the while loop don't
%     continue endlessly in some cases

if nargin == 2
  NE = 3;
end

% The vertices of the big cube containing P
Min = double(min(P));
Max = double(max(P));

% Number of cubes with edge length "EdgeLength" in the sides
% of the big cube
N = double(ceil((Max-Min)/EL)+2*NE+1);
t = 0;
while t < 10 && 8*N(1)*N(2)*N(3) > 4e9
  t = t+1;
  EL = 1.1*EL;
  N = double(ceil((Max-Min)/EL)+2*NE+1);
end
if 8*N(1)*N(2)*N(3) > 4e9
  NE = 3;
  N = double(ceil((Max-Min)/EL)+2*NE+1);
end
Info = [Min N EL NE];

% Calculates the cube-coordinates of the points
CubeCoord = floor([P(:,1)-Min(1) P(:,2)-Min(2) P(:,3)-Min(3)]/EL)+NE+1;

% Sorts the points according a lexicographical order
LexOrd = [CubeCoord(:,1) CubeCoord(:,2)-1 CubeCoord(:,3)-1]*[1 N(1) N(1)*N(2)]';
CubeCoord = uint16(CubeCoord);
[LexOrd,SortOrd] = sort(LexOrd);
SortOrd = uint32(SortOrd);
LexOrd = uint32(LexOrd);

if nargout <= 3
  % Define "Partition"
  Partition = cell(N(1),N(2),N(3));
  np = size(P,1);     % number of points
  p = 1;              % The index of the point under comparison
  while p <= np
    t = 1;
    while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
      t = t+1;
    end
    q = SortOrd(p);
    Partition{CubeCoord(q,1),CubeCoord(q,2),CubeCoord(q,3)} = SortOrd(p:p+t-1);
    p = p+t;
  end

else
  nc = size(unique(LexOrd),1);

  % Define "Partition"
  Cubes = zeros(N(1),N(2),N(3),'uint32');
  Partition = cell(nc,1);
  np = size(P,1);     % number of points
  p = 1;              % The index of the point under comparison
  c = 0;
  while p <= np
    t = 1;
    while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
      t = t+1;
    end
    q = SortOrd(p);
    c = c+1;
    Partition{c,1} = SortOrd(p:p+t-1);
    Cubes(CubeCoord(q,1),CubeCoord(q,2),CubeCoord(q,3)) = c;
    p = p+t;
  end
end
