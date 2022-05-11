function [SurfCov,Dis,CylVol,dis] = surface_coverage(P,Axis,Point,nl,ns,Dmin,Dmax)
 
% ---------------------------------------------------------------------
% SURFACE_COVERAGE.M   Computes point surface coverage measure
%
% Version       1.1.0
% Last update   7 Oct 2021
%
% Copyright (C) 2017-2021 Pasi Raumonen
% ---------------------------------------------------------------------
% Inputs:    
% Axis      Axis direction (1 x 3) 
% Point     Starting point of the cylinder (1 x 3)
% nl        Number of layers in the axis direction used for to partition
%               the cylinder surface into layer/sectors
% ns        Number of sectors used to partition the cylinder surface into 
%               layer/sectors
% Dmin      (Optional) Minimum point distance from the axis to be included
%               into SurfCov calculations
% Dmax      (Optional) Maximum point distance from the axis to be included
%               into SurfCov calculations
% 
% Output:  
% SurfCov   Number between 0 and 1 descring how big portion of the cylinder
%               surface is covered with points
% Dis       (Optional) Mean distances of the distances of the layer/sectors
% CylVol    (Optional) Volume of the cylinder estimated by the mean
%               distances of the layer/sectors as cylindrical sections
% dis       (Optional) Same as "Dis" but empty cells are interpolated
% ---------------------------------------------------------------------
% Computes surface coverage (number between 0 and 1) of points on cylinder 
% surface defined by "Axis" and "Point".

% Changes from version 1.0.0 to 1.1.0, 7 Oct 2021:
% 1) Added two possible inputs, minimum and maximum distance, 
%    Dmin and Dmax, which can be used to filter out points for the surface
%    coverage calculations
% 2) Computes the SurfCov estimate with four baseline directions used in
%    the sector determination and selects the largest value
% 3) Smalle changes to speed up computations

%% Compute the distances and heights of the points
[d,V,h] = distances_to_line(P,Axis,Point);
h = h-min(h);
Len = max(h);

%% (Optional) Filter out points based on the distance to the axis
if nargin >= 6
  Keep = d > Dmin;
  if nargin == 7
    Keep = Keep & d < Dmax;
  end
  V = V(Keep,:);
  h = h(Keep);
end

%% Compute SurfCov
% from 4 different baseline directions to determine the angles and select
% the maximum value
V0 = V;
[U,W] = orthonormal_vectors(Axis); % First planar axes
R = rotation_matrix(Axis,2*pi/ns/4); % Rotation matrix to rotate the axes
SurfCov = zeros(1,4);
for i = 1:4
  %% Rotate the axes
  if i > 1
    U = R*U;
    W = R*W;
  end
  
  %% Compute the angles (sectors) of the points
  V = V0*[U W];
  ang = atan2(V(:,2),V(:,1))+pi;
  
  %% Compute lexicographic order (sector,layer) of every point
  Layer = ceil(h/Len*nl);
  Layer(Layer <= 0) = 1;
  Layer(Layer > nl) = nl;
  Sector = ceil(ang/2/pi*ns);
  Sector(Sector <= 0) = 1;
  LexOrd = [Layer Sector-1]*[1 nl]';
  
  %% Compute SurfCov
  Cov = zeros(nl,ns);
  Cov(LexOrd) = 1;
  SurfCov(i) = nnz(Cov)/nl/ns;
end
SurfCov = max(SurfCov);


%% Compute volume estimate
if nargout > 1
  % Sort according to increasing lexicographic order
  [LexOrd,SortOrd] = sort(LexOrd);
  d = d(SortOrd);
  
  % Compute mean distance of the sector-layer intersections
  Dis = zeros(nl,ns); % mean distances
  np = length(LexOrd);     % number of points
  p = 1;
  while p <= np
    t = 1;
    while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
      t = t+1;
    end
    Dis(LexOrd(p)) = average(d(p:p+t-1));
    p = p+t;
  end
  
  if nargout > 2
    % Interpolate missing distances
    D = Dis;
    dis = Dis;
    Dinv = D((nl:-1:1)',:);
    D = [Dinv Dinv Dinv; D D D; Dinv Dinv Dinv];
    Zero = Dis == 0;
    RadMean = average(Dis(Dis > 0));
    for i = 1:nl
      for j = 1:ns
        if Zero(i,j)
          if nnz(D(i+nl-1:i+nl+1,j+ns-1:j+ns+1)) > 1
            d = D(i+nl-1:i+nl+1,j+ns-1:j+ns+1);
            dis(i,j) = average(d(d > 0));
          elseif nnz(D(i+nl-2:i+nl+2,j+ns-2:j+ns+2)) > 1
            d = D(i+nl-2:i+nl+2,j+ns-2:j+ns+2);
            dis(i,j) = average(d(d > 0));
          elseif nnz(D(i+nl-3:i+nl+3,j+ns-3:j+ns+3)) > 1
            d = D(i+nl-3:i+nl+3,j+ns-3:j+ns+3);
            dis(i,j) = average(d(d > 0));
          else
            dis(i,j) = RadMean;
          end
        end
      end
    end
    % Compute the volume estimate
    r = dis(:);
    CylVol = 1000*pi*sum(r.^2)/ns*Len/nl;
  end
end
