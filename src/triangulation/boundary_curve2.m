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

function Curve = boundary_curve2(P,Curve0,rball,dmax)

% ---------------------------------------------------------------------
% BOUNDARY_CURVE2.M      Determines the boundary curve based on the
%                           previously defined boundary curve.
%
% Version 1.0
% Latest update     16 Aug 2017
%
% Copyright (C) 2015-2017 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Inputs:
% P         Point cloud of the cross section
% Curve0    Seed points from previous cross section curve
% rball     Radius of the balls centered at seed points
% dmax      Maximum distance between concecutive curve points, if larger,
%               then create a new one between the points


%% Partition the point cloud into cubes
Min = double(min([P(:,1:2); Curve0(:,1:2)]));
Max = double(max([P(:,1:2); Curve0(:,1:2)]));
N = double(ceil((Max-Min)/rball)+5);
CC = floor([P(:,1)-Min(1) P(:,2)-Min(2)]/rball)+3; % cube coordinates of the section points
% Sorts the points according a lexicographical order
S = [CC(:,1) CC(:,2)-1]*[1 N(1)]';
[S,I] = sort(S);
% Define "partition"
np = size(P,1);
partition = cell(N(1),N(2));
p = 1;              % The index of the point under comparison
while p <= np
  t = 1;
  while (p+t <= np) && (S(p) == S(p+t))
    t = t+1;
  end
  q = I(p);
  partition{CC(q,1),CC(q,2)} = I(p:p+t-1);
  p = p+t;
end


%% Define segments using the previous points
CC = floor([Curve0(:,1)-Min(1) Curve0(:,2)-Min(2)]/rball)+3; % cube coordinates of the seed points
I = CC < 3;
CC(I) = 3;
nc = size(Curve0,1);  % number of sets
Dist = 1e8*ones(np,1);  % distance of point to the closest center
SoP = zeros(np,1);  % the segment the points belong to
Radius = rball^2;
for i = 1:nc
  points = partition(CC(i,1)-1:CC(i,1)+1,CC(i,2)-1:CC(i,2)+1);
  points = vertcat(points{:});
  V = [P(points,1)-Curve0(i,1) P(points,2)-Curve0(i,2)];
  dist = sum(V.*V,2);
  PointsInBall = dist < Radius;
  points = points(PointsInBall);
  dist = dist(PointsInBall);
  D = Dist(points);
  L = dist < D;
  I = points(L);
  Dist(I) = dist(L);
  SoP(I) = i;
end

%% Finalise the segments
% Number of points in each segment and index of each point in its segment
Num = zeros(nc,1);
IndPoints = zeros(np,1);
for i = 1:np
  if SoP(i) > 0
    Num(SoP(i)) = Num(SoP(i))+1;
    IndPoints(i) = Num(SoP(i));
  end
end
% Continue if enough non-emtpy segments
if nnz(Num) > 0.05*nc
  % Initialization of the "Seg"
  Seg = cell(nc,1);
  for i = 1:nc
    Seg{i} = zeros(Num(i),1);
  end
  % Define the "Seg"
  for i = 1:np
    if SoP(i) > 0
      Seg{SoP(i),1}(IndPoints(i),1) = i;
    end
  end

  %% Define the new curve points as the average of the segments
  Curve = zeros(nc,3);  % the new boundary curve
  for i = 1:nc
    S = Seg{i};
    if ~isempty(S)
      Curve(i,:) = mean(P(S,:),1);
      if norm(Curve(i,:)-Curve0(i,:)) > 1.25*dmax
        Curve(i,:) = Curve0(i,:);
      end
    else
      Curve(i,:) = Curve0(i,:);
    end
  end

  %% Add new points if too large distances
  V = Curve([2:end 1],:)-Curve(1:end,:);
  d = sum(V.*V,2);
  Large = d > dmax^2;
  m = nnz(Large);
  if m > 0
    Curve0 = zeros(nc+m,3);
    t = 0;
    for i = 1:nc
      if Large(i)
        t = t+1;
        Curve0(t,:) = Curve(i,:);
        t = t+1;
        Curve0(t,:) = Curve(i,:)+0.5*V(i,:);
      else
        t = t+1;
        Curve0(t,:) = Curve(i,:);
      end
    end
    Curve = Curve0;
  end

  %% Remove new points if too small distances
  nc = size(Curve,1);
  V = Curve([2:end 1],:)-Curve(1:end,:);
  d = sum(V.*V,2);
  Small = d < (0.333*dmax)^2;
  m = nnz(Small);
  if m > 0
    for i = 1:nc-1
      if Small(i) && Small(i+1)
        Small(i+1) = false;
      end
    end
    if ~Small(nc) && Small(1)
      Small(1) = false;
      Small(nc) = true;
    end
    Curve = Curve(~Small,:);
  end

else
  % If not enough new points, return the old curve
  Curve = Curve0;
end
