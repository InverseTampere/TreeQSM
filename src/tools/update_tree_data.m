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

function treedata = update_tree_data(QSM,cylinder,branch,inputs)

% ---------------------------------------------------------------------
% UPDATE_TREE_DATA.M       Updates the treedata structure, e.g. after
%                            simplification of QSM
%
% Version 1.0.0
% Latest update     4 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------

% Inputs:
% treedata    Treedata structure from "tree_data"
% cylinder    Cylinder structure from "cylinders"
% branch      Branch structure from "branches"
%
% Output:
% treedata     Tree data/attributes in a struct
% ---------------------------------------------------------------------


% Define some variables from cylinder:
treedata = QSM.treedata;
Rad = cylinder.radius;
Len = cylinder.length;
Sta = cylinder.start;
Axe = cylinder.axis;
nc = length(Rad);
ind = (1:1:nc)';
Trunk = cylinder.branch == 1; % Trunk cylinders

%% Tree attributes from cylinders
% Volumes, areas, lengths, branches
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
treedata.TotalLength = treedata.TrunkLength+treedata.BranchLength;
NB = length(branch.order)-1; % number of branches
treedata.NumberBranches = NB;
BO = max(branch.order); % maximum branch order
treedata.MaxBranchOrder = BO;
treedata.TrunkArea = 2*pi*sum(Rad(Trunk).*Len(Trunk));
treedata.BranchArea = 2*pi*sum(Rad(~Trunk).*Len(~Trunk));
treedata.TotalArea = 2*pi*sum(Rad.*Len);

%% Crown measures,Vertical profile and spreads
[treedata,spreads] = crown_measures(treedata,cylinder,branch);

%% Update triangulation information
if inputs.Tria
  treedata = update_triangulation(QSM,treedata,cylinder);
end

%% Tree Location
treedata.location = Sta(1,:);

%% Stem taper
R = Rad(Trunk);
n = length(R);
Taper = zeros(n+1,2);
Taper(1,2) = 2*R(1);
Taper(2:end,1) = cumsum(Len(Trunk));
Taper(2:end,2) = [2*R(2:end); 2*R(n)];
treedata.StemTaper = Taper';

%% Vertical profile and spreads
treedata.VerticalProfile = mean(spreads,2);
treedata.spreads = spreads;

%% CYLINDER DISTRIBUTIONS:
%% Wood part diameter distributions
% Volume, area and length of wood parts as functions of cylinder diameter
% (in 1cm diameter classes)
treedata = cylinder_distribution(treedata,Rad,Len,Axe,'Dia');

%% Wood part height distributions
% Volume, area and length of cylinders as a function of height
% (in 1 m height classes)
treedata = cylinder_height_distribution(treedata,Rad,Len,Sta,Axe,ind);

%% Wood part zenith direction distributions
% Volume, area and length of wood parts as functions of cylinder zenith
% direction (in 10 degree angle classes)
treedata = cylinder_distribution(treedata,Rad,Len,Axe,'Zen');

%% Wood part azimuth direction distributions
% Volume, area and length of wood parts as functions of cylinder zenith
% direction (in 10 degree angle classes)
treedata = cylinder_distribution(treedata,Rad,Len,Axe,'Azi');

%% BRANCH DISTRIBUTIONS:
%% Branch order distributions
% Volume, area, length and number of branches as a function of branch order
treedata = branch_order_distribution(treedata,branch);

%% Branch diameter distributions
% Volume, area, length and number of branches as a function of branch diameter
% (in 1cm diameter classes)
treedata = branch_distribution(treedata,branch,'Dia');

%% Branch height distribution
% Volume, area, length and number of branches as a function of branch height
% (in 1 meter classes) for all and 1st-order branches
treedata = branch_distribution(treedata,branch,'Hei');

%% Branch angle distribution
% Volume, area, length and number of branches as a function of branch angle
% (in 10 deg angle classes) for all and 1st-order branches
treedata = branch_distribution(treedata,branch,'Ang');

%% Branch azimuth distribution
% Volume, area, length and number of branches as a function of branch azimuth
% (in 22.5 deg angle classes) for all and 1st-order branches
treedata = branch_distribution(treedata,branch,'Azi');

%% Branch zenith distribution
% Volume, area, length and number of branches as a function of branch zenith
% (in 10 deg angle classes) for all and 1st-order branches
treedata = branch_distribution(treedata,branch,'Zen');

%% change into single-format
Names = fieldnames(treedata);
n = size(Names,1);
for i = 1:n
  treedata.(Names{i}) = single(treedata.(Names{i}));
end

if inputs.disp == 2
  %% Generate units for displaying the treedata
  Units = zeros(n,3);
  m = 23;
  for i = 1:n
    if ~inputs.Tria && strcmp(Names{i},'CrownVolumeAlpha')
      m = i;
    elseif inputs.Tria && strcmp(Names{i},'TriaTrunkLength')
      m = i;
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
    elseif strcmp(Names{i}(end-4:end),'aConv')
      Units(i,:) = 'm^2';
    elseif strcmp(Names{i}(end-5:end),'aAlpha')
      Units(i,:) = 'm^2';
    elseif strcmp(Names{i}(end-4:end),'eConv')
      Units(i,:) = 'm^3';
    elseif strcmp(Names{i}(end-5:end),'eAlpha')
      Units(i,:) = 'm^3';
    elseif strcmp(Names{i}(end-2:end),'Ave')
      Units(i,:) = 'm  ';
    elseif strcmp(Names{i}(end-2:end),'Max')
      Units(i,:) = 'm  ';
    end
  end
  %% Display treedata
  disp('------------')
  disp('  Tree attributes:')
  for i = 1:m
    v = change_precision(treedata.(Names{i}));
    if strcmp(Names{i},'DBHtri')
      disp('  -----')
      disp('  Tree attributes from triangulation:')
    end
    disp(['  ',Names{i},' = ',num2str(v),' ',Units(i,:)])
  end
  disp('  -----')
end

if inputs.plot > 1
  %% Plot distributions
  figure(6)
  subplot(2,4,1)
  plot(Taper(:,1),Taper(:,2),'-b')
  title('Stem taper')
  xlabel('Distance from base (m)')
  ylabel('Diameter (m)')
  axis tight
  grid on
  
  Q.treedata = treedata;
  subplot(2,4,2)
  plot_distribution(Q,6,0,0,'VolCylDia')
  
  subplot(2,4,3)
  plot_distribution(Q,6,0,0,'AreCylDia')
  
  subplot(2,4,4)
  plot_distribution(Q,6,0,0,'LenCylDia')
  
  subplot(2,4,5)
  plot_distribution(Q,6,0,0,'VolBranchOrd')
  
  subplot(2,4,6)
  plot_distribution(Q,6,0,0,'LenBranchOrd')
  
  subplot(2,4,7)
  plot_distribution(Q,6,0,0,'AreBranchOrd')
  
  subplot(2,4,8)
  plot_distribution(Q,6,0,0,'NumBranchOrd')
  
  figure(7)
  subplot(3,3,1)
  plot_distribution(Q,7,0,0,'VolCylHei')
  
  subplot(3,3,2)
  plot_distribution(Q,7,0,0,'AreCylHei')
  
  subplot(3,3,3)
  plot_distribution(Q,7,0,0,'LenCylHei')
  
  subplot(3,3,4)
  plot_distribution(Q,7,0,0,'VolCylZen')
  
  subplot(3,3,5)
  plot_distribution(Q,7,0,0,'AreCylZen')
  
  subplot(3,3,6)
  plot_distribution(Q,7,0,0,'LenCylZen')
  
  subplot(3,3,7)
  plot_distribution(Q,7,0,0,'VolCylAzi')
  
  subplot(3,3,8)
  plot_distribution(Q,7,0,0,'AreCylAzi')
  
  subplot(3,3,9)
  plot_distribution(Q,7,0,0,'LenCylAzi')
  
  figure(8)
  subplot(3,4,1)
  plot_distribution(Q,8,1,0,'VolBranchDia','VolBranch1Dia')
  
  subplot(3,4,2)
  plot_distribution(Q,8,1,0,'AreBranchDia','AreBranch1Dia')
  
  subplot(3,4,3)
  plot_distribution(Q,8,1,0,'LenBranchDia','LenBranch1Dia')
  
  subplot(3,4,4)
  plot_distribution(Q,8,1,0,'NumBranchDia','NumBranch1Dia')
  
  subplot(3,4,5)
  plot_distribution(Q,8,1,0,'VolBranchHei','VolBranch1Hei')
  
  subplot(3,4,6)
  plot_distribution(Q,8,1,0,'AreBranchHei','AreBranch1Hei')
  
  subplot(3,4,7)
  plot_distribution(Q,8,1,0,'LenBranchHei','LenBranch1Hei')
  
  subplot(3,4,8)
  plot_distribution(Q,8,1,0,'NumBranchHei','NumBranch1Hei')
  
  subplot(3,4,9)
  plot_distribution(Q,8,1,0,'VolBranchAng','VolBranch1Ang')
  
  subplot(3,4,10)
  plot_distribution(Q,8,1,0,'AreBranchAng','AreBranch1Ang')
  
  subplot(3,4,11)
  plot_distribution(Q,8,1,0,'LenBranchAng','LenBranch1Ang')
  
  subplot(3,4,12)
  plot_distribution(Q,8,1,0,'NumBranchAng','NumBranch1Ang')
  
  figure(9)
  subplot(2,4,1)
  plot_distribution(Q,9,1,0,'VolBranchZen','VolBranch1Zen')
  
  subplot(2,4,2)
  plot_distribution(Q,9,1,0,'AreBranchZen','AreBranch1Zen')
  
  subplot(2,4,3)
  plot_distribution(Q,9,1,0,'LenBranchZen','LenBranch1Zen')
  
  subplot(2,4,4)
  plot_distribution(Q,9,1,0,'NumBranchZen','NumBranch1Zen')
  
  subplot(2,4,5)
  plot_distribution(Q,9,1,0,'VolBranchAzi','VolBranch1Azi')
  
  subplot(2,4,6)
  plot_distribution(Q,9,1,0,'AreBranchAzi','AreBranch1Azi')
  
  subplot(2,4,7)
  plot_distribution(Q,9,1,0,'LenBranchAzi','LenBranch1Azi')
  
  subplot(2,4,8)
  plot_distribution(Q,9,1,0,'NumBranchAzi','NumBranch1Azi')
end

end % End of main function


function [treedata,spreads] = crown_measures(treedata,cylinder,branch)

%% Generate point clouds from the cylinder model
Axe = cylinder.axis;
Len = cylinder.length;
Sta = cylinder.start;
Tip = Sta+[Len.*Axe(:,1) Len.*Axe(:,2) Len.*Axe(:,3)]; % tips of the cylinders
nc = length(Len);
P = zeros(5*nc,3); % four mid points on the cylinder surface
t = 0;
for i = 1:nc
  [U,V] = orthonormal_vectors(Axe(i,:));
  U = cylinder.radius(i)*U;
  if cylinder.branch(i) == 1
    % For stem cylinders generate more points
    for k = 1:4
      M = Sta(i,:)+k*Len(i)/4*Axe(i,:);
      R = rotation_matrix(Axe(i,:),pi/12);
      for j = 1:12
        if j > 1
          U = R*U;
        end
        t = t+1;
        P(t,:) = M+U';
      end
    end
  else
    M = Sta(i,:)+Len(i)/2*Axe(i,:);
    R = rotation_matrix(Axe(i,:),pi/4);
    for j = 1:4
      if j > 1
        U = R*U;
      end
      t = t+1;
      P(t,:) = M+U';
    end
  end
end
P = P(1:t,:);
P = double([P; Sta; Tip]);
P = unique(P,'rows');

%% Vertical profiles (layer diameters/spreads), mean:
bot = min(P(:,3));
top = max(P(:,3));
Hei = top-bot;
if Hei > 10
  m = 20;
elseif Hei > 2
  m = 10;
else
  m = 5;
end
spreads = zeros(m,18);
for j = 1:m
  I = P(:,3) >= bot+(j-1)*Hei/m & P(:,3) < bot+j*Hei/m;
  X = unique(P(I,:),'rows');
  if size(X,1) > 5
    [K,A] = convhull(X(:,1),X(:,2));
    % compute center of gravity for the convex hull and use it as
    % center for computing average diameters
    n = length(K);
    x = X(K,1);
    y = X(K,2);
    CX = sum((x(1:n-1)+x(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
    CY = sum((y(1:n-1)+y(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
    
    V = mat_vec_subtraction(X(:,1:2),[CX CY]);
    ang = atan2(V(:,2),V(:,1))+pi;
    [ang,I] = sort(ang);
    L = sqrt(sum(V.*V,2));
    L = L(I);
    for i = 1:18
      I = ang >= (i-1)*pi/18 & ang < i*pi/18;
      if any(I)
        L1 = max(L(I));
      else
        L1 = 0;
      end
      J = ang >= (i-1)*pi/18+pi & ang < i*pi/18+pi;
      if any(J)
        L2 = max(L(J));
      else
        L2 = 0;
      end
      spreads(j,i) = L1+L2;
    end
  end
end

%% Crown diameters (spreads), mean and maximum:
X = unique(P(:,1:2),'rows');
[K,A] = convhull(X(:,1),X(:,2));
% compute center of gravity for the convex hull and use it as center for
% computing average diameters
n = length(K);
x = X(K,1);
y = X(K,2);
CX = sum((x(1:n-1)+x(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
CY = sum((y(1:n-1)+y(2:n)).*(x(1:n-1).*y(2:n)-x(2:n).*y(1:n-1)))/6/A;
V = mat_vec_subtraction(Tip(:,1:2),[CX CY]);
ang = atan2(V(:,2),V(:,1))+pi;
[ang,I] = sort(ang);
L = sqrt(sum(V.*V,2));
L = L(I);
S = zeros(18,1);
for i = 1:18
  I = ang >= (i-1)*pi/18 & ang < i*pi/18;
  if any(I)
    L1 = max(L(I));
  else
    L1 = 0;
  end
  J = ang >= (i-1)*pi/18+pi & ang < i*pi/18+pi;
  if any(J)
    L2 = max(L(J));
  else
    L2 = 0;
  end
  S(i) = L1+L2;
end
treedata.CrownDiamAve = mean(S);
MaxDiam = 0;
for i = 1:n
  V = mat_vec_subtraction([x y],[x(i) y(i)]);
  L = max(sqrt(sum(V.*V,2)));
  if L > MaxDiam
    MaxDiam = L;
  end
end
treedata.CrownDiamMax = L;

%% Crown areas from convex hull and alpha shape:
treedata.CrownAreaConv = A;
alp = max(0.5,treedata.CrownDiamAve/10);
shp = alphaShape(X(:,1),X(:,2),alp);
treedata.CrownAreaAlpha = shp.area;

%% Crown base
% Define first major branch as the branch whose diameter > min(0.05*dbh,5cm)
% and whose horizontal relative reach is more than the median reach of 1st-ord.
% branches (or at maximum 10). The reach is defined as the horizontal
% distance from the base to the tip divided by the dbh.
dbh = treedata.DBHcyl;
nb = length(branch.order);
HL = zeros(nb,1); % horizontal reach
branches1 = (1:1:nb)';
branches1 = branches1(branch.order == 1); % 1st-order branches
nb = length(branches1);
nc = size(Sta,1);
ind = (1:1:nc)';
for i = 1:nb
  C = ind(cylinder.branch == branches1(i));
  if ~isempty(C)
    base = Sta(C(1),:);
    C = C(end);
    tip = Sta(C,:)+Len(C)*Axe(C);
    V = tip(1:2)-base(1:2);
    HL(branches1(i)) = sqrt(V*V')/dbh*2;
  end
end
M = min(10,median(HL));

% Sort the branches according to the their heights
Hei = branch.height(branches1);
[Hei,SortOrd] = sort(Hei);
branches1 = branches1(SortOrd);

% Search the first/lowest branch:  
d = min(0.05,0.05*dbh);
b = 0;
if nb > 1
  i = 1;
  while i < nb
    i = i+1;
    if branch.diameter(branches1(i)) > d && HL(branches1(i)) > M
      b = branches1(i);
      i = nb+2;
    end
  end
  if i == nb+1 && nb > 1
    b = branches1(1);
  end
end

if b > 0
  % search all the children of the first major branch:
  nb = size(branch.parent,1);
  Ind = (1:1:nb)';
  chi = Ind(branch.parent == b);
  B = b;
  while ~isempty(chi)
    B = [B; chi];
    n = length(chi);
    C = cell(n,1);
    for i = 1:n
      C{i} = Ind(branch.parent == chi(i));
    end
    chi = vertcat(C{:});
  end
  
  % define crown base height from the ground:
  BaseHeight = max(Sta(:,3)); % Height of the crown base
  for i = 1:length(B)
    C = ind(cylinder.branch == B(i));
    ht = min(Tip(C,3));
    hb = min(Sta(C,3));
    h = min(hb,ht);
    if h < BaseHeight
      BaseHeight = h;
    end
  end
  treedata.CrownBaseHeight = BaseHeight-Sta(1,3);
  
  %% Crown length and ratio
  treedata.CrownLength = treedata.TreeHeight-treedata.CrownBaseHeight;
  treedata.CrownRatio = treedata.CrownLength/treedata.TreeHeight;
  
  %% Crown volume from convex hull and alpha shape:
  I = P(:,3) >= BaseHeight;
  X = P(I,:);
  [K,V] = convhull(X(:,1),X(:,2),X(:,3));
  treedata.CrownVolumeConv = V;
  alp = max(0.5,treedata.CrownDiamAve/5);
  shp = alphaShape(X(:,1),X(:,2),X(:,3),alp,'HoleThreshold',10000);
  treedata.CrownVolumeAlpha = shp.volume;

else 
  % No branches
  treedata.CrownBaseHeight = treedata.TreeHeight;
  treedata.CrownLength = 0;
  treedata.CrownRatio = 0;
  treedata.CrownVolumeConv = 0;
  treedata.CrownVolumeAlpha = 0;
end

end % End of function


function treedata = update_triangulation(QSM,treedata,cylinder)

% Update the mixed results:
if ~isempty(QSM.triangulation)
  CylInd = QSM.triangulation.cylind;
  Rad = cylinder.radius;
  Len = cylinder.length;
  % Determine the new stem cylinder that is about the location where the
  % triangulation stops:
  nc = length(Rad);
  ind = (1:1:nc)';
  ind = ind(cylinder.branch == 1); % cylinders in the stem
  S = QSM.cylinder.start(CylInd,:); % The place where the triangulation stops
  V = cylinder.start(ind,:)-S;
  d = sqrt(sum(V.*V,2));
  [d,I] = min(d);
  V = V(I,:);
  CylInd = ind(I); % The new cylinder closest to the correct place
  if d < 0.01
    TrunkVolMix = treedata.TrunkVolume-...
      1000*pi*sum(Rad(1:CylInd-1).^2.*Len(1:CylInd-1))+QSM.triangulation.volume;
    TrunkAreaMix = treedata.TrunkArea-...
      2*pi*sum(Rad(1:CylInd-1).*Len(1:CylInd-1))+QSM.triangulation.SideArea;
  else
    % Select the following cylinder
    h = V*cylinder.axis(CylInd,:)';
    if h < 0
      CylInd = CylInd+1;
      V = cylinder.start(CylInd,:)-S;
      h = V*cylinder.axis(CylInd,:)';
    end
    Len(CylInd-1) = Len(CylInd-1)-h;
    
    TrunkVolMix = treedata.TrunkVolume-...
      1000*pi*sum(Rad(1:CylInd-1).^2.*Len(1:CylInd-1))+QSM.triangulation.volume;
    TrunkAreaMix = treedata.TrunkArea-...
      2*pi*sum(Rad(1:CylInd-1).*Len(1:CylInd-1))+QSM.triangulation.SideArea;
  end
  treedata.MixTrunkVolume = TrunkVolMix;
  treedata.MixTotalVolume = TrunkVolMix+treedata.BranchVolume;
  treedata.MixTrunkArea = TrunkAreaMix;
  treedata.MixTotalArea = TrunkAreaMix+treedata.BranchArea;
end
end


function treedata = cylinder_distribution(treedata,Rad,Len,Axe,dist)
%% Wood part diameter, zenith and azimuth direction distributions
% Volume, area and length of wood parts as functions of cylinder
% diameter, zenith, and azimuth
if strcmp(dist,'Dia')
  Par = Rad;
  n = ceil(max(200*Rad));
  a = 0.005; % diameter in 1 cm classes
elseif strcmp(dist,'Zen')
  Par = 180/pi*acos(Axe(:,3));
  n = 18;
  a = 10; % zenith direction in 10 degree angle classes
elseif strcmp(dist,'Azi')
  Par = 180/pi*atan2(Axe(:,2),Axe(:,1))+180;
  n = 36;
  a = 10; % azimuth direction in 10 degree angle classes
end

CylDist = zeros(3,n);
for i = 1:n
  K = Par >= (i-1)*a & Par < i*a;
  CylDist(1,i) = 1000*pi*sum(Rad(K).^2.*Len(K)); % volumes in litres
  CylDist(2,i) = 2*pi*sum(Rad(K).*Len(K)); % areas in litres
  CylDist(3,i) = sum(Len(K)); % lengths in meters
end
treedata.(['VolCyl',dist]) = CylDist(1,:);
treedata.(['AreCyl',dist]) = CylDist(2,:);
treedata.(['LenCyl',dist]) = CylDist(3,:);
end


function treedata = cylinder_height_distribution(treedata,Rad,Len,Sta,Axe,ind)

%% Wood part height distributions
% Volume, area and length of cylinders as a function of height
% (in 1 m height classes)
MaxHei= ceil(treedata.TreeHeight);
treedata.VolCylHei = zeros(1,MaxHei);
treedata.AreCylHei = zeros(1,MaxHei);
treedata.LenCylHei = zeros(1,MaxHei);
End = Sta+[Len.*Axe(:,1) Len.*Axe(:,2) Len.*Axe(:,3)];
bot = min(Sta(:,3));
B = Sta(:,3)-bot;
T = End(:,3)-bot;
for j = 1:MaxHei
  I1 = B >= (j-2) & B < (j-1); % base below this bin
  J1 = B >= (j-1) & B < j; % base in this bin
  K1 = B >= j & B < (j+1); % base above this bin
  I2 = T >= (j-2) & T < (j-1); % top below this bin
  J2 = T >= (j-1) & T < j; % top in this bin
  K2 = T >= j & T < (j+1); % top above this bin
  C1 = ind(J1&J2); % base and top in this bin
  C2 = ind(J1&K2); % base in this bin, top above
  C3 = ind(J1&I2); % base in this bin, top below
  C4 = ind(I1&J2); % base in bin below, top in this
  C5 = ind(K1&J2); % base in bin above, top in this
  v1 = 1000*pi*sum(Rad(C1).^2.*Len(C1));
  a1 = 2*pi*sum(Rad(C1).*Len(C1));
  l1 = sum(Len(C1));
  r2 = (j-B(C2))./(T(C2)-B(C2)); % relative portion in this bin
  v2 = 1000*pi*sum(Rad(C2).^2.*Len(C2).*r2);
  a2 = 2*pi*sum(Rad(C2).*Len(C2).*r2);
  l2 = sum(Len(C2).*r2);
  r3 = (B(C3)-j+1)./(B(C3)-T(C3)); % relative portion in this bin
  v3 = 1000*pi*sum(Rad(C3).^2.*Len(C3).*r3);
  a3 = 2*pi*sum(Rad(C3).*Len(C3).*r3);
  l3 = sum(Len(C3).*r3);
  r4 = (T(C4)-j+1)./(T(C4)-B(C4)); % relative portion in this bin
  v4 = 1000*pi*sum(Rad(C4).^2.*Len(C4).*r4);
  a4 = 2*pi*sum(Rad(C4).*Len(C4).*r4);
  l4 = sum(Len(C4).*r4);
  r5 = (j-T(C5))./(B(C5)-T(C5)); % relative portion in this bin
  v5 = 1000*pi*sum(Rad(C5).^2.*Len(C5).*r5);
  a5 = 2*pi*sum(Rad(C5).*Len(C5).*r5);
  l5 = sum(Len(C5).*r5);
  treedata.VolCylHei(j) = v1+v2+v3+v4+v5;
  treedata.AreCylHei(j) = a1+a2+a3+a4+a5;
  treedata.LenCylHei(j) = l1+l2+l3+l4+l5;
end
end


function treedata = branch_distribution(treedata,branch,dist)
%% Branch diameter, height, angle, zenith and azimuth distributions
% Volume, area, length and number of branches as a function of branch
% diamater, height, angle, zenith and aximuth
BOrd = branch.order(2:end);
BVol = branch.volume(2:end);
BAre = branch.area(2:end);
BLen = branch.length(2:end);
if strcmp(dist,'Dia')
  Par = branch.diameter(2:end);
  n = ceil(max(100*Par));
  a = 0.005; % diameter in 1 cm classes
elseif strcmp(dist,'Hei')
  Par = branch.height(2:end);
  n = ceil(treedata.TreeHeight);
  a = 1; % height in 1 m classes
elseif strcmp(dist,'Ang')
  Par = branch.angle(2:end);
  n = 18;
  a = 10; % angle in 10 degree classes
elseif strcmp(dist,'Zen')
  Par = branch.zenith(2:end);
  n = 18;
  a = 10; % zenith direction in 10 degree angle classes
elseif strcmp(dist,'Azi')
  Par = branch.azimuth(2:end)+180;
  n = 36;
  a = 10; % azimuth direction in 10 degree angle classes
end

BranchDist = zeros(8,n);
for i = 1:n
  I = Par >= (i-1)*a & Par < i*a;
  BranchDist(1,i) = sum(BVol(I)); % volume (all branches)
  BranchDist(2,i) = sum(BVol(I & BOrd == 1)); % volume (1st-branches)
  BranchDist(3,i) = sum(BAre(I)); % area (all branches)
  BranchDist(4,i) = sum(BAre(I & BOrd == 1)); % area (1st-branches)
  BranchDist(5,i) = sum(BLen(I)); % length (all branches)
  BranchDist(6,i) = sum(BLen(I & BOrd == 1)); % length (1st-branches)
  BranchDist(7,i) = nnz(I); % number (all branches)
  BranchDist(8,i) = nnz(I & BOrd == 1); % number (1st-branches)
end
treedata.(['VolBranch',dist]) = BranchDist(1,:);
treedata.(['VolBranch1',dist]) = BranchDist(2,:);
treedata.(['AreBranch',dist]) = BranchDist(3,:);
treedata.(['AreBranch1',dist]) = BranchDist(4,:);
treedata.(['LenBranch',dist]) = BranchDist(5,:);
treedata.(['LenBranch1',dist]) = BranchDist(6,:);
treedata.(['NumBranch',dist]) = BranchDist(7,:);
treedata.(['NumBranch1',dist]) = BranchDist(8,:);
end


function treedata = branch_order_distribution(treedata,branch)
%% Branch order distributions
% Volume, area, length and number of branches as a function of branch order
BO = max(branch.order);
BranchOrdDist = zeros(BO,4);
for i = 1:max(1,BO)
  I = branch.order == i;
  BranchOrdDist(i,1) = sum(branch.volume(I)); % volumes
  BranchOrdDist(i,2) = sum(branch.area(I)); % areas
  BranchOrdDist(i,3) = sum(branch.length(I)); % lengths
  BranchOrdDist(i,4) = nnz(I); % number of ith-order branches
end
treedata.VolBranchOrd = BranchOrdDist(:,1)';
treedata.AreBranchOrd = BranchOrdDist(:,2)';
treedata.LenBranchOrd = BranchOrdDist(:,3)';
treedata.NumBranchOrd = BranchOrdDist(:,4)';
end
