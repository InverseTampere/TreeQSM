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

function cover = cover_sets(P,inputs,RelSize)

% ---------------------------------------------------------------------
% COVER_SETS.M          Creates cover sets (surface patches) and their
%                       neighbor-relation for a point cloud
%
% Version 2.0.1
% Latest update     2 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------

% Covers the point cloud with small sets, which are along the surface,
% such that each point belongs at most one cover set; i.e. the cover is
% a partition of the point cloud.
%
% The cover is generated such that at first the point cloud is covered
% with balls with radius "BallRad". This first cover is such that
% 1) the minimum distance between the centers is "PatchDiam", and
% 2) the maximum distance from any point to nearest center is also "PatchDiam".
% Then the first cover of BallRad-balls is used to define a second cover:
% each BallRad-ball "A" defines corresponding cover set "B" in the second cover
% such that "B" contains those points of "A" that are nearer to the center of
% "A" than any other center of BallRad-balls. The BallRad-balls also define
% the neighbors for the second cover: Let CA and CB denote cover sets in
% the second cover, and BA and BB their BallRad-balls. Then CB is
% a neighbor of CA, and vice versa, if BA and CB intersect or
% BB and CA intersect.
%
% Inputs:
% P         Point cloud
% inputs    Input stucture, the following fields are needed:
%   PatchDiam1   Minimum distance between centers of cover sets; i.e. the
%                   minimum diameter of cover set in uniform covers. Does
%                   not need nor use the third optional input "RelSize".
%   PatchDiam2Min   Minimum diameter of cover sets for variable-size
%                     covers. Needed if "RelSize" is given as input.
%   PatchDiam2Max   Maximum diameter of cover sets for variable-size
%                     covers. Needed if "RelSize" is given as input.
% 	BallRad1    Radius of the balls used to generate the uniform cover. 
%                   These balls are also used to determine the neighbors
%   BallRad2    Maximum radius of the balls used to generate the 
%                   varibale-size cover. 
%   nmin1, nmin2    Minimum number of points in a BallRad1- and
%                       BallRad2-balls
% RelSize   Relative cover set size for each point
%
% Outputs:
% cover     Structure array containing the followin fields:
%   ball        Cover sets, (n_sets x 1)-cell
%   center      Center points of the cover sets, (n_sets x 1)-vector
%   neighbor    Neighboring cover sets of each cover set, (n_sets x 1)-cell

% Changes from version 2.0.0 to 2.0.1, 2 May 2022:
% 1) Added comments and changed some variable names
% 2) Enforced that input parameters are type double

if ~isa(P,'double')
  P = double(P);
end

%% Large balls and centers
np = size(P,1);
Ball = cell(np,1); % Large balls for generation of the cover sets and their neighbors
Cen = zeros(np,1,'uint32'); % the center points of the balls/cover sets
NotExa = true(np,1); % the points not yet examined
Dist = 1e8*ones(np,1);  % distance of point to the closest center
BoP = zeros(np,1,'uint32');  % the balls/cover sets the points belong
nb = 0;             % number of sets generated
if nargin == 2
  %% Same size cover sets everywhere
  BallRad = double(inputs.BallRad1);
  PatchDiamMax = double(inputs.PatchDiam1);
  nmin = double(inputs.nmin1);
  % Partition the point cloud into cubes for quick neighbor search
  [partition,CC] = cubical_partition(P,BallRad);

  % Generate the balls
  Radius = BallRad^2;
  MaxDist = PatchDiamMax^2;
  % random permutation of points, produces different covers for the same inputs:
  RandPerm = randperm(np); 
  for i = 1:np
    if NotExa(RandPerm(i))
      Q = RandPerm(i); % the center/seed point of the current cover set
      % Select the points in the cubical neighborhood of the seed:
      points = partition(CC(Q,1)-1:CC(Q,1)+1,CC(Q,2)-1:CC(Q,2)+1,CC(Q,3)-1:CC(Q,3)+1);
      points = vertcat(points{:});
      % Compute distances of the points to the seed:
      V = [P(points,1)-P(Q,1) P(points,2)-P(Q,2) P(points,3)-P(Q,3)];
      dist = sum(V.*V,2);
      % Select the points inside the ball:
      Inside = dist < Radius;
      if nnz(Inside) >= nmin
        ball = points(Inside); % the points forming the ball
        d = dist(Inside); % the distances of the ball's points
        core = (d < MaxDist); % the core points of the cover set
        NotExa(ball(core)) = false; % mark points as examined
        % define new ball:
        nb = nb+1; 
        Ball{nb} = ball;
        Cen(nb) = Q;
        % Select which points belong to this ball, i.e. are closer this
        % seed than previously tested seeds:
        D = Dist(ball); % the previous distances
        closer = d < D; % which points are closer to this seed
        ball = ball(closer); % define the ball
        % update the ball and distance information of the points
        Dist(ball) = d(closer); 
        BoP(ball) = nb; 
      end
    end
  end
else
  %% Use relative sizes (the size varies)
  % Partition the point cloud into cubes
  BallRad = double(inputs.BallRad2);
  PatchDiamMin = double(inputs.PatchDiam2Min);
  PatchDiamMax = double(inputs.PatchDiam2Max);
  nmin = double(inputs.nmin2);
  MRS = PatchDiamMin/PatchDiamMax;
  % minimum radius
  r = double(1.5*(double(min(RelSize))/256*(1-MRS)+MRS)*BallRad+1e-5); 
  NE = 1+ceil(BallRad/r);
  if NE > 4
    r = PatchDiamMax/4;
    NE = 1+ceil(BallRad/r);
  end
  [Partition,CC,~,Cubes] = cubical_partition(P,r,NE);

  I = RelSize == 0; % Don't use points with no size determined
  NotExa(I) = false;

  % Define random permutation of points (results in different covers for 
  % same input) so that first small sets are generated
  RandPerm = zeros(np,1,'uint32');
  I = RelSize <= 32;
  ind = uint32(1:1:np)';
  I = ind(I);
  t1 = length(I);
  RandPerm(1:1:t1) = I(randperm(t1));
  I = RelSize <= 128 & RelSize > 32;
  I = ind(I);
  t2 = length(I);
  RandPerm(t1+1:1:t1+t2) = I(randperm(t2));
  t2 = t2+t1;
  I = RelSize > 128;
  I = ind(I);
  t3 = length(I);
  RandPerm(t2+1:1:t2+t3) = I(randperm(t3));
  clearvars ind I

  Point = zeros(round(np/1000),1,'uint32');
  e = BallRad-PatchDiamMax;
  for i = 1:np
    if NotExa(RandPerm(i))
      Q = RandPerm(i); % the center/seed point of the current cover set
      % Compute the set size and the cubical neighborhood of the seed point:
      rs = double(RelSize(Q))/256*(1-MRS)+MRS; % relative radius
      MaxDist = PatchDiamMax*rs; % diameter of the cover set
      Radius = MaxDist+sqrt(rs)*e; % radius of the ball including the cover set
      N = ceil(Radius/r); % = number of cells needed to include the ball
      cubes = Cubes(CC(Q,1)-N:CC(Q,1)+N,CC(Q,2)-N:CC(Q,2)+N,CC(Q,3)-N:CC(Q,3)+N);
      I = cubes > 0;
      cubes = cubes(I); % Cubes forming the neighborhood
      Par = Partition(cubes); % cell-array of the points in the neighborhood
      % vertical catenation of the points from the cell-array
      S = cellfun('length',Par);
      stop = cumsum(S);
      start = [0; stop]+1;
      for k = 1:length(stop)
        Point(start(k):stop(k)) = Par{k};
      end
      points = Point(1:stop(k));
      % Compute the distance of the "points" to the seed:
      V = [P(points,1)-P(Q,1) P(points,2)-P(Q,2) P(points,3)-P(Q,3)];
      dist = sum(V.*V,2);
      % Select the points inside the ball:
      Inside = dist < Radius^2;
      if nnz(Inside) >= nmin
        ball = points(Inside); % the points forming the ball
        d = dist(Inside); % the distances of the ball's points
        core = (d < MaxDist^2); % the core points of the cover set
        NotExa(ball(core)) = false; % mark points as examined
        % define new ball:
        nb = nb+1; 
        Ball{nb} = ball;
        Cen(nb) = Q;
        % Select which points belong to this ball, i.e. are closer this
        % seed than previously tested seeds:
        D = Dist(ball); % the previous distances
        closer = d < D; % which points are closer to this seed
        ball = ball(closer); % define the ball
        % update the ball and distance information of the points
        Dist(ball) = d(closer); 
        BoP(ball) = nb; 
      end
    end
  end
end
Ball = Ball(1:nb,:);
Cen = Cen(1:nb);
clearvars RandPerm NotExa Dist

%% Cover sets
% Number of points in each ball and index of each point in its ball
Num = zeros(nb,1,'uint32');
Ind = zeros(np,1,'uint32');
for i = 1:np
  if BoP(i) > 0
    Num(BoP(i)) = Num(BoP(i))+1;
    Ind(i) = Num(BoP(i));
  end
end

% Initialization of the "PointsInSets"
PointsInSets = cell(nb,1);
for i = 1:nb
  PointsInSets{i} = zeros(Num(i),1,'uint32');
end

% Define the "PointsInSets"
for i = 1:np
  if BoP(i) > 0
    PointsInSets{BoP(i),1}(Ind(i)) = i;
  end
end

%% Neighbors
% Define neighbors. Sets A and B are neighbors if the large ball of A
% contains points of B. Notice that this is not a symmetric relation.
Nei = cell(nb,1);
Fal = false(nb,1);
for i = 1:nb
  B = Ball{i};        % the points in the big ball of cover set "i"
  I = (BoP(B) ~= i);
  N = B(I);           % the points of B not in the cover set "i"
  N = BoP(N);

  % select the unique elements of N:
  n = length(N);
  if n > 2
    Include = true(n,1);
    for j = 1:n
      if ~Fal(N(j))
        Fal(N(j)) = true;
      else
        Include(j) = false;
      end
    end
    Fal(N) = false;
    N = N(Include);
  elseif n == 2
    if N(1) == N(2)
      N = N(1);
    end
  end

  Nei{i} = uint32(N);
end

% Make the relation symmetric by adding, if needed, A as B's neighbor
% in the case B is A's neighbor
for i = 1:nb
  N = Nei{i};
  for j = 1:length(N)
    K = (Nei{N(j)} == i);
    if ~any(K)
      Nei{N(j)} = uint32([Nei{N(j)}; i]);
    end
  end
end

% Define output
cover.ball = PointsInSets;
cover.center = Cen;
cover.neighbor = Nei;

%% Display statistics
%disp(['    ',num2str(nb),' cover sets, points not covered: ',num2str(np-nnz(BoP))])