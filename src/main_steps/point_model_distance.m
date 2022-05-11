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

function pmdistance = point_model_distance(P,cylinder)

% ---------------------------------------------------------------------
% POINT_MODEL_DISTANCE.M    Computes the distances of the points to the 
%                               cylinder model
%
% Version 2.1.1
% Latest update     8 Oct 2021
%
% Copyright (C) 2015-2021 Pasi Raumonen
% ---------------------------------------------------------------------

% Changes from version 2.1.0 to 2.1.1, 8 Oct 2021:  
% 1) Changed the determinationa NE, the number of empty edge layers, so 
%     that is now limited in size, before it is given as input for 
%     cubical_partition function.

% Changes from version 2.0.0 to 2.1.0, 26 Nov 2019:  
% 1) Bug fix: Corrected the computation of the output at the end of the
%    code so that trees without branches are computed correctly.

% Cylinder data
Rad = cylinder.radius;
Len = cylinder.length;
Sta = cylinder.start;
Axe = cylinder.axis;
BOrd = cylinder.BranchOrder;

% Select randomly 25 % or max one million points for the distance comput.
np0 = size(P,1);
a = min(0.25*np0,1000000);
I = logical(round(0.5/(1-a/np0)*rand(np0,1)));
P = P(I,:);

% Partition the points into cubes 
L = 2*median(Len);
NE = max(3,min(10,ceil(max(Len)/L)))+3;
[Partition,~,Info] = cubical_partition(P,L,NE);
Min = Info(1:3);
EL = Info(7);
NE = Info(8);

% Calculates the cube-coordinates of the starting points
CC = floor([Sta(:,1)-Min(1) Sta(:,2)-Min(2) Sta(:,3)-Min(3)]/EL)+NE+1;

% Compute the number of cubes needed for each starting point
N = ceil(Len/L);

% Correct N so that cube indexes are not too small or large
I = CC(:,1) < N+1;
N(I) = CC(I,1)-1;
I = CC(:,2) < N+1;
N(I) = CC(I,2)-1;
I = CC(:,3) < N+1;
N(I) = CC(I,3)-1;
I = CC(:,1)+N+1 > Info(4);
N(I) = Info(4)-CC(I,1)-1;
I = CC(:,2)+N+1 > Info(5);
N(I) = Info(5)-CC(I,2)-1;
I = CC(:,3)+N+1 > Info(6);
N(I) = Info(6)-CC(I,3)-1;

% Calculate the distances to the cylinders
n = size(Rad,1);
np = size(P,1);
Dist = zeros(np,2); % Distance and the closest cylinder of each points
Dist(:,1) = 2; % Large distance initially
Points = zeros(ceil(np/10),1,'int32'); % Auxiliary variable
Data = cell(n,1);
for i = 1:n
  Par = Partition(CC(i,1)-N(i):CC(i,1)+N(i),CC(i,2)-N(i):CC(i,2)+N(i),...
    CC(i,3)-N(i):CC(i,3)+N(i));
  if N(i) > 1
    S = cellfun('length',Par);
    I = S > 0;
    S = S(I);
    Par = Par(I);
    stop = cumsum(S);
    start = [0; stop]+1;
    for k = 1:length(stop)
      Points(start(k):stop(k)) = Par{k}(:);
    end
    points = Points(1:stop(k));
  else
    points = vertcat(Par{:});
  end
  [d,~,h] = distances_to_line(P(points,:),Axe(i,:),Sta(i,:));
  d = abs(d-Rad(i));
  Data{i} = [d h double(points)];
  I = d < Dist(points,1);
  J = h >= 0;
  K = h <= Len(i);
  L = d < 0.5;
  M = I&J&K&L;
  points = points(M);
  Dist(points,1) = d(M);
  Dist(points,2) = i;
end

% Calculate the distances to the cylinders for points not yet calculated
% because they are not "on side of cylinder
for i = 1:n
  if ~isempty(Data{i})
    d = Data{i}(:,1);
    h = Data{i}(:,2);
    points = Data{i}(:,3);
    I = d < Dist(points,1);
    J = h >= -0.1 & h <= 0;
    K = h <= Len(i)+0.1 & h >= Len(i);
    L = d < 0.5;
    M = I&(J|K)&L;
    points = points(M);
    Dist(points,1) = d(M);
    Dist(points,2) = i;
  end
end

% Select only the shortest 95% of distances for each cylinder
N = zeros(n,1);
O = zeros(np,1);
for i = 1:np
  if Dist(i,2) > 0
    N(Dist(i,2)) = N(Dist(i,2))+1;
    O(i) = N(Dist(i,2));
  end
end
Cyl = cell(n,1);
for i = 1:n
  Cyl{i} = zeros(N(i),1);
end
for i = 1:np
  if Dist(i,2) > 0
    Cyl{Dist(i,2)}(O(i)) = i;
  end
end
DistCyl = zeros(n,1); % Average point distance to each cylinder
for i = 1:n
  I = Cyl{i};
  m = length(I);
  if m > 19 % select the smallest 95% of distances
    d = sort(Dist(I,1));
    DistCyl(i) = mean(d(1:floor(0.95*m)));
  elseif m > 0
    DistCyl(i) = mean(Dist(I,1));
  end
end

% Define the output
pmdistance.CylDist = single(DistCyl);
pmdistance.median = median(DistCyl(:,1));
pmdistance.mean = mean(DistCyl(:,1));
pmdistance.max = max(DistCyl(:,1));
pmdistance.std = std(DistCyl(:,1));

T = BOrd == 0;
B1 = BOrd == 1;
B2 = BOrd == 2;
B = DistCyl(~T,1);
T = DistCyl(T,1);
B1 = DistCyl(B1,1);
B2 = DistCyl(B2,1);

pmdistance.TrunkMedian = median(T);
pmdistance.TrunkMean = mean(T);
pmdistance.TrunkMax = max(T);
pmdistance.TrunkStd = std(T);

if ~isempty(B)
  pmdistance.BranchMedian = median(B);
  pmdistance.BranchMean = mean(B);
  pmdistance.BranchMax = max(B);
  pmdistance.BranchStd = std(B);
else
  pmdistance.BranchMedian = 0;
  pmdistance.BranchMean = 0;
  pmdistance.BranchMax = 0;
  pmdistance.BranchStd = 0;
end

if ~isempty(B1)
  pmdistance.Branch1Median = median(B1);
  pmdistance.Branch1Mean = mean(B1);
  pmdistance.Branch1Max = max(B1);
  pmdistance.Branch1Std = std(B1);
else
  pmdistance.Branch1Median = 0;
  pmdistance.Branch1Mean = 0;
  pmdistance.Branch1Max = 0;
  pmdistance.Branch1Std = 0;
end

if ~isempty(B2)
  pmdistance.Branch2Median = median(B2);
  pmdistance.Branch2Mean = mean(B2);
  pmdistance.Branch2Max = max(B2);
  pmdistance.Branch2Std = std(B2);
else
  pmdistance.Branch2Median = 0;
  pmdistance.Branch2Mean = 0;
  pmdistance.Branch2Max = 0;
  pmdistance.Branch2Std = 0;
end
