function [Pass,c] = surface_coverage_filtering(P,c,lh,ns)

% Filter a 3d-point cloud based on given cylinder (axis and radius) by 
% dividing the point cloud into "ns" equel-angle sectors and "nl" equal-length
% layers along the axis. For each sector-layer intersection (a region in 
% the cylinder surface) keep only the points closest to the centre. 

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
    D = ds(p:p+t-1);
    I = D <= min(D)+0.01;
    Dis(LexOrd(p)) = mean(D(I));
    p = p+t;
end

% Compute the number of sectors based on the estimated radius
R = median(Dis(Dis > 0));
a = max(0.02,0.2*R);
ns = ceil(2*pi*R/a);
ns = max(ns,10);
ns = min(ns,36);
nl = ceil(c.length/a);
nl = max(nl,4);

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
Part = cell(nl*ns,1);
Dis = zeros(nl,ns);
Pass = false(np,1);
p = 1;
k = 0;
r = max(0.01,0.05*R);
a = 0;
while p <= np
    t = 1;
    while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
        t = t+1;
    end
    a = a+1;
    Part{a} = p:p+t-1;
    D = d(p:p+t-1);
    I = D <= min(D)+r;
    ind = p:p+t-1;
    Pass(ind(I)) = true;
    Dis(LexOrd(p)) = mean(D(I));
    p = p+t;
    k = k+1;
end

n = length(SortOrd);
InvSortOrd = zeros(n,1);
InvSortOrd(SortOrd) = (1:1:n)';
d = d(Pass);
Pass = Pass(InvSortOrd);

% Compute radius, SurfCov and mad
R = median(Dis(Dis > 0));
mad = mean(abs(d-R));

c.radius = R;
c.SurfCov = k/nl/ns;
c.mad = mad;
c.conv = 1;
c.rel = 1;

% if R < 0.12
% P = P(SortOrd,:);
% plot_segs(P,Part,2,15)
% end