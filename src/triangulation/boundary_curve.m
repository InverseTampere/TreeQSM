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

function [Curve,Ind] = boundary_curve(P,Curve0,rball,dmax)

% ---------------------------------------------------------------------
% BOUNDARY_CURVE.M      Determines the boundary curve based on the
%                           previously defined boundary curve.
%
% Version 1.1.0
% Latest update     3 May 2022
%
% Copyright (C) 2015-2022 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Inputs:
% P         Point cloud of the cross section
% Curve0     Seed points from previous cross section curve
% rball     Radius of the balls centered at seed points
% dmax      Maximum distance between concecutive curve points, if larger,
%               then create a new one between the points
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 1.1.0, 3 May 2022:
% 1) Increased the cubical neighborhood in the generation of the segments

%% Partition the point cloud into cubes
Min = double(min([P(:,1:2); Curve0(:,1:2)]));
Max = double(max([P(:,1:2); Curve0(:,1:2)]));
N = double(ceil((Max-Min)/rball)+5);
% cube coordinates of the section points
CC = floor([P(:,1)-Min(1) P(:,2)-Min(2)]/rball)+3;
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
% cube coordinates of the seed points:
CC = floor([Curve0(:,1)-Min(1) Curve0(:,2)-Min(2)]/rball)+3;
I = CC < 3;
CC(I) = 3;
nc = size(Curve0,1);  % number of sets
Dist = 1e8*ones(np,1);  % distance of point to the closest center
SoP = zeros(np,1);  % the segment the points belong to
Radius = rball^2;
for i = 1:nc
  points = partition(CC(i,1)-2:CC(i,1)+2,CC(i,2)-2:CC(i,2)+2);
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
  Empty = false(nc,1);
  for i = 1:nc
    S = Seg{i};
    if ~isempty(S)
      Curve(i,:) = mean(P(S,:),1);
      if norm(Curve(i,:)-Curve0(i,:)) > 1.25*dmax
        Curve(i,:) = Curve0(i,:);
      end
    else
      Empty(i) = true;
    end
  end

  %% Interpolate for empty segments
  % For empty segments create points by interpolation from neighboring 
  % non-empty segments
  if any(Empty)
    for i = 1:nc
      if Empty(i)
        if i > 1 && i < nc
          k = 0;
          while i+k <= nc && Empty(i+k)
            k = k+1;
          end
          if i+k <= nc
            LineEle = Curve(i+k,:)-Curve(i-1,:);
          else
            LineEle = Curve(1,:)-Curve(i-1,:);
          end
          if k < 5
            for j = 1:k
              Curve(i+j-1,:) = Curve(i-1,:)+j/(k+1)*LineEle;
            end
          else
            Curve(i:i+k-1,:) = Curve0(i:i+k-1,:);
          end
        elseif i == 1
          a = 0;
          while Empty(end-a)
            a = a+1;
          end
          b = 1;
          while Empty(b)
            b = b+1;
          end
          LineEle = Curve(b,:)-Curve(nc-a,:);
          n = a+b-1;
          if n < 5
            for j = 1:a-1
              Curve(nc-a+1+j,:) = Curve(nc-a,:)+j/n*LineEle;
            end
            for j = 1:b-1
              Curve(j,:) = Curve(nc-a,:)+(j+a-1)/n*LineEle;
            end
          else
            Curve(nc-a+2:nc,1:2) = Curve0(nc-a+2:nc,1:2);
            Curve(nc-a+2:nc,3) = Curve0(nc-a+2:nc,3);
            Curve(1:b-1,1:2) = Curve0(1:b-1,1:2);
            Curve(1:b-1,3) = Curve0(1:b-1,3);
          end
        elseif i == nc
          LineEle = Curve(1,:)-Curve(nc-1,:);
          Curve(i,:) = Curve(nc-1,:)+0.5*LineEle;
        end
      end
    end
  end

  % Correct the height
  Curve(:,3) = min(Curve(:,3));

  % Check self-intersection
  [Intersect,IntersectLines] = check_self_intersection(Curve(:,1:2));

  % If self-intersection, try to modify the curve
  j = 1;
  while Intersect && j <= 5
    n = size(Curve,1);
    InterLines = (1:1:n)';
    NumberOfIntersections = cellfun('length',IntersectLines(:,1));
    I = NumberOfIntersections > 0;
    InterLines = InterLines(I);
    CrossLen = vertcat(IntersectLines{I,2});
    if length(CrossLen) == length(InterLines)
      LineEle = Curve([2:end 1],:)-Curve(1:end,:);
      d = sqrt(sum(LineEle.*LineEle,2));
      m = length(InterLines);
      for i = 1:2:m
        if InterLines(i) ~= n
          Curve(InterLines(i)+1,:) = Curve(InterLines(i),:)+...
            0.9*CrossLen(i)/d(InterLines(i))*LineEle(InterLines(i),:);
        else
          Curve(1,:) = Curve(InterLines(i),:)+...
            0.9*CrossLen(i)/d(InterLines(i))*LineEle(InterLines(i),:);
        end
      end
      [Intersect,IntersectLines] = check_self_intersection(Curve(:,1:2));
      j = j+1;
    else
      j = 6;
    end
  end

  %% Add new points if too large distances
  LineEle = Curve([2:end 1],:)-Curve(1:end,:);
  d = sum(LineEle.*LineEle,2);
  Large = d > dmax^2;
  m = nnz(Large);
  if m > 0
    Curve0 = zeros(nc+m,3);
    Ind = zeros(nc+m,2);
    t = 0;
    for i = 1:nc
      if Large(i)
        t = t+1;
        Curve0(t,:) = Curve(i,:);
        if i < nc
          Ind(t,:) = [i i+1];
        else
          Ind(t,:) = [i 1];
        end
        t = t+1;
        Curve0(t,:) = Curve(i,:)+0.5*LineEle(i,:);
        if i < nc
          Ind(t,:) = [i+1 0];
        else
          Ind(t,:) = [1 0];
        end
      else
        t = t+1;
        Curve0(t,:) = Curve(i,:);
        if i < nc
          Ind(t,:) = [i i+1];
        else
          Ind(t,:) = [i 1];
        end
      end
    end
    Curve = Curve0;

  else
    Ind = [(1:1:nc)' [(2:1:nc)'; 1]];
  end


  %% Remove new points if too small distances
  nc = size(Curve,1);
  LineEle = Curve([2:end 1],:)-Curve(1:end,:);
  d = sum(LineEle.*LineEle,2);
  Small = d < (0.333*dmax)^2;
  m = nnz(Small);
  if m > 0
    for i = 1:nc-1
      if ~Small(i) && Small(i+1)
        Ind(i,2) = -1;
      elseif Small(i) && Small(i+1)
        Small(i+1) = false;
      end
    end
    if ~Small(nc) && Small(1)
      Ind(nc,2) = -1;
      Ind(1,2) = -1;
      Small(1) = false;
      Small(nc) = true;
      I = Ind(:,2) > 0;
      Ind(2:end,1) = Ind(2:end,1)+1;
      Ind(I,2) = Ind(I,2)+1;

    end
    Ind = Ind(~Small,:);
    Curve = Curve(~Small,:);
  end

else
  % If not enough new points, return the old curve
  Ind = [(1:1:nc)' [(2:1:nc)'; 1]];
  Curve = Curve0;
end
