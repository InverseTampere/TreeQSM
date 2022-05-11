function [Pass,c] = surface_coverage_filtering(P,c,lh,ns)

% ---------------------------------------------------------------------
% SURFACE_COVERAGE_FILTERING.M    Filters a point cloud based on the 
%                                   assumption that it samples a cylinder
%
% Version 1.1.0
% Latest update     6 Oct 2021
%
% Copyright (C) 2017-2021 Pasi Raumonen
% ---------------------------------------------------------------------

% Filter a 3d-point cloud based on given cylinder (axis and radius) by 
% dividing the point cloud into "ns" equal-angle sectors and "lh"-height
% layers along the axis. For each sector-layer intersection (a region in 
% the cylinder surface) keep only the points closest to the axis. 

% Inputs:
% P             Point cloud, (n_points x 3)-matrix
% c             Cylinder, stucture array with fields "axis", "start",
%                 "length"
% lh            Height of the layers
% ns            Number of sectors
%
% Outputs:              
% Pass          Logical vector indicating which points pass the filtering
% c             Cylinder, stucture array with additional fields "radius",
%                 "SurfCov", "mad", "conv", "rel", estimated from the
%                 filtering
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 1.1.0, 6 Oct 2021:  
% 1) Small changes to make the code little faster
% 2) Change the radius estimation to make it much faster
 

% Compute the distances, heights and angles of the points
[d,V,h] = distances_to_line(P,c.axis,c.start);
h = h-min(h);
[U,W] = orthonormal_vectors(c.axis);
V = V*[U W];
ang = atan2(V(:,2),V(:,1))+pi;

% Sort based on lexicographic order of (sector,layer)
nl = ceil(c.length/lh);
Layer = ceil(h/c.length*nl);
Layer(Layer == 0) = 1;
Layer(Layer > nl) = nl;
Sector = ceil(ang/2/pi*ns);
Sector(Sector == 0) = 1;
LexOrd = [Layer Sector-1]*[1 nl]';
[LexOrd,SortOrd] = sort(LexOrd);
ds = d(SortOrd);

% Estimate the distances for each sector-layer intersection
Dis = zeros(nl,ns);
np = size(P,1);     % number of points
p = 1;
while p <= np
  t = 1;
  while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
    t = t+1;
  end
  D = min(ds(p:p+t-1));
  Dis(LexOrd(p)) = min(1.05*D,D+0.02);
  p = p+t;
end

% Compute the number of sectors based on the estimated radius
R = median(Dis(Dis > 0));
a = max(0.02,0.2*R);
ns = ceil(2*pi*R/a);
ns = min(36,max(ns,8));
nl = ceil(c.length/a);
nl = max(nl,3);

% Sort based on lexicographic order of (sector,layer)
Layer = ceil(h/c.length*nl);
Layer(Layer == 0) = 1;
Layer(Layer > nl) = nl;
Sector = ceil(ang/2/pi*ns);
Sector(Sector == 0) = 1;
LexOrd = [Layer Sector-1]*[1 nl]';
[LexOrd,SortOrd] = sort(LexOrd);
d = d(SortOrd);

% Filtering for each sector-layer intersection
Dis = zeros(nl,ns);
Pass = false(np,1);
p = 1; % index of point under processing
k = 0; % number of nonempty cells 
r = max(0.01,0.05*R); % cell diameter from the closest point
while p <= np
  t = 1;
  while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
    t = t+1;
  end
  ind = p:p+t-1;
  D = d(ind);
  Dmin = min(D);
  I = D <= Dmin+r;
  Pass(ind(I)) = true;
  Dis(LexOrd(p)) = min(1.05*Dmin,Dmin+0.02);
  p = p+t;
  k = k+1;
end
d = d(Pass);

% Sort the "Pass"-vector back to original point cloud order
n = length(SortOrd);
InvSortOrd = zeros(n,1);
InvSortOrd(SortOrd) = (1:1:n)';
Pass = Pass(InvSortOrd);

% Compute radius, SurfCov and mad
R = median(Dis(Dis > 0));
mad = sum(abs(d-R))/length(d);

c.radius = R;
c.SurfCov = k/nl/ns;
c.mad = mad;
c.conv = 1;
c.rel = 1;
