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

function Curve = initial_boundary_curve(P,TriaWidth)

% ---------------------------------------------------------------------
% INITIAL_BOUNDARY_CURVE.M      Determines the boundary curve adaptively.
%
% Version 1.0.1
% Latest update     26 Nov 2019
%
% Copyright (C) 2015-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 1.0.1, 26 Nov 2019:
% 1) Bug fix: Added "return" if the "Curve" is empty after it is first defined.

%% Define suitable center
% Use xy-data and even the z-coordinate to the top
Top = max(P(:,3));
P = [P(:,1:2) Top*ones(size(P,1),1)];

% Define the "center" of points as the mean
Center = mean(P);
Center0 = Center;

% If the center is outside or close to the boundary, define new center
i = 0;
A0 = 61;
ShortestDist = 0;
while ShortestDist < 0.075 && i < 100
  Center = Center0+[3*ShortestDist*randn(1,2) 0]; % Randomly move the center
  % Compute angles of points as seen from the center
  V = mat_vec_subtraction(P(:,1:2),Center(1:2));
  angle = 180/pi*atan2(V(:,2),V(:,1))+180;
  % % Check if the center is outside or near the boundary of the cross section
  A = false(70,1);
  a = ceil(angle/5);
  I = a > 0;
  A(a(I)) = true;
  if i == 0
    ShortestDist = 0.025;
  elseif nnz(A) < A0
    ShortestDist = 0.05;
  else
    PointDist = sqrt(sum(V.*V,2));
    [ShortestDist,FirstPoint] = min(PointDist);
  end
  i = i+1;
  if i == 100 && ShortestDist < 0.075
    i = 0;
    A0 = A0-2;
  end
end

%% Define first boundary curve based on the center
Curve = zeros(18,1); % the boundary curve, contains indexed of the point cloud rows
Curve(1) = FirstPoint; % start the curve from the point the closest the center
% Modify the angles so that first point has the angle 0
a0 = angle(FirstPoint);
I = angle < a0;
angle(I) = angle(I)+(360-a0);
angle(~I) = angle(~I)-a0;
% Select the rest of the points as the closest point in 15 deg sectors
% centered at 20 deg intervals
np = size(P,1);
Ind = (1:1:np)';
t = 0;
for i = 2:18
  J = angle > 12.5+20*(i-2) & angle < 27.5+20*(i-2);
  if ~any(J) % if no points, try 18 deg sector
    J = angle > 11+20*(i-2) & angle < 29+20*(i-2);
  end
  if any(J)
    % if sector has points, select the closest point as the curve point
    D = PointDist(J);
    ind = Ind(J);
    [~,J] = min(D);
    t = t+1;
    Curve(t) = ind(J);
  end
end
Curve = Curve(1:t);
if isempty(Curve)
  return
end
I = true(np,1);
I(Curve) = false;
Ind = Ind(I);


%% Adapt the initial curve to the data
V = P(Curve([(2:t)'; 1]),:)-P(Curve,:);
D = sqrt(sum(V(:,1:2).*V(:,1:2),2));
n = t;
n0 = 1;
% Continue adding new points as long as too long edges exists
while any(D > 1.25*TriaWidth) && n > n0
  N = [V(:,2) -V(:,1) V(:,3)];
  M = P(Curve,:)+0.5*V;

  Curve1 = Curve;
  t = 0;
  for i = 1:n
    if D(i) > 1.25*TriaWidth
      [d,~,hc] = distances_to_line(P(Curve1,:),N(i,:),M(i,:));
      I = hc > 0.01 & d < D(i)/2;
      if any(I)
        H = min(hc(I));
      else
        H = 1;
      end
      [d,~,h] = distances_to_line(P(Ind,:),N(i,:),M(i,:));
      I = d < D(i)/3 & h > -TriaWidth/2 & h < H;

      if any(I)
        ind = Ind(I);
        h = h(I);
        [h,J] = min(h);
        I = ind(J);

        t = t+1;
        if i < n
          Curve1 = [Curve1(1:t); I; Curve1(t+1:end)];
        else
          Curve1 = [Curve1(1:t); I];
        end
        J = Ind ~= I;
        Ind = Ind(J);
        t = t+1;

      else
        t = t+1;
      end
    else
      t = t+1;
    end
  end
  Curve = Curve1(1:t);

  n0 = n;
  n = size(Curve,1);
  V = P(Curve([(2:n)'; 1]),:)-P(Curve,:);
  D = sqrt(sum(V.*V,2));
end

%% Refine the curve for longer edges if far away points
n0 = n-1;
while n > n0
  N = [V(:,2) -V(:,1) V(:,3)];
  M = P(Curve,:)+0.5*V;

  Curve1 = Curve;
  t = 0;
  for i = 1:n
    if D(i) > 0.5*TriaWidth
      [d,~,hc] = distances_to_line(P(Curve1,:),N(i,:),M(i,:));
      I = hc > 0.01 & d < D(i)/2;
      if any(I)
        H = min(hc(I));
      else
        H = 1;
      end
      [d,~,h] = distances_to_line(P(Ind,:),N(i,:),M(i,:));
      I = d < D(i)/3 & h > -TriaWidth/3 & h < H;
      ind = Ind(I);
      h = h(I);
      [h,J] = min(h);

      if h > TriaWidth/10
        I = ind(J);
        t = t+1;
        if i < n
          Curve1 = [Curve1(1:t); I; Curve1(t+1:end)];
        else
          Curve1 = [Curve1(1:t); I];
        end
        J = Ind ~= I;
        Ind = Ind(J);
        t = t+1;

      else
        t = t+1;
      end
    else
      t = t+1;
    end

  end
  Curve = Curve1(1:t);

  n0 = n;
  n = size(Curve,1);
  V = P(Curve([(2:n)'; 1]),:)-P(Curve,:);
  D = sqrt(sum(V.*V,2));
end

%% Smooth the curve by defining the points by means of neighbors
Curve = P(Curve,:); % Change the curve from point indexes to coordinates
Curve = boundary_curve2(P,Curve,0.04,TriaWidth);
if isempty(Curve)
  return
end

%% Add points for too long edges
n = size(Curve,1);
V = Curve([(2:n)'; 1],:)-Curve;
D = sqrt(sum(V.*V,2));
Curve1 = Curve;
t = 0;
for i = 1:n
  if D(i) > TriaWidth
    m = floor(D(i)/TriaWidth);
    t = t+1;
    W = zeros(m,3);
    for j = 1:m
      W(j,:) = Curve(i,:)+j/(m+1)*V(i,:);
    end
    Curve1 = [Curve1(1:t,:); W; Curve1(t+1:end,:)];
    t = t+m ;
  else
    t = t+1;
  end
end
Curve = Curve1;
n = size(Curve,1);

%% Define the curve again by equalising the point distances along the curve
V = Curve([(2:n)'; 1],:)-Curve;
D = sqrt(sum(V.*V,2));
L = cumsum(D);
m = ceil(L(end)/TriaWidth);
TriaWidth = L(end)/m;
Curve1 = zeros(m,3);
Curve1(1,:) = Curve(1,:);
b = 1;
for i = 2:m
  while L(b) < (i-1)*TriaWidth
    b = b+1;
  end
  if b > 1
    a = ((i-1)*TriaWidth-L(b-1))/D(b);
    Curve1(i,:) = Curve(b,:)+a*V(b,:);
  else
    a = (L(b)-(i-1)*TriaWidth)/D(b);
    Curve1(i,:) = Curve(b,:)+a*V(b,:);
  end
end
Curve = Curve1;

Intersect = check_self_intersection(Curve(:,1:2));
if Intersect
  Curve = zeros(0,3);
end

