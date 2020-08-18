function [SurfCov,Dis,CylVol] = surface_coverage(P,Axis,Point,nl,ns)
 
% ---------------------------------------------------------------------
% SURFACE_COVERAGE.M   Computes point surface coverage measure
%
% Version 1.0.0
% Created     14 Jan 2020
%
% Copyright (C) 2017-2020 Pasi Raumonen
% ---------------------------------------------------------------------
% Inputs:    
% Axis      Axis direction (1 x 3) 
% Point     Starting point of the cylinder (1 x 3)
% nl        Number of layers in the axis direction used for to partition
%               the cylinder surface into layer/sectors
% ns        Number of sectors used to partition the cylinder surface into 
%               layer/sectors
% 
% Output:  
% SurfCov   Number between 0 and 1 descring how big portion of the cylinder
%               surface is covered with points
% Dis       (Optional) Mean distances of the distances of the layer/sectors
% CylVol    (Optional) Volume of the cylinder estimated by the mean
%               distances of the layer/sectors as cylindrical sections
% ---------------------------------------------------------------------
% Computes surface coverage (number between 0 and 1) of points on cylinder 
% surface defined by "Axis" and "Point".

%% Compute the distances and heights of the points
[d,V,h] = distances_to_line(P,Axis,Point);
h = h-min(h);
Len = max(h);

%% Compute the angles (sectors) of the points
[U,W] = orthonormal_vectors(Axis);
V = V*[U W];
ang = atan2(V(:,2),V(:,1))+pi;

%% Compute lexicographic order (sector,layer) of every point
Layer = ceil(h/Len*nl);
Layer(Layer <= 0) = 1;
Layer(Layer > nl) = nl;
Sector = ceil(ang/2/pi*ns);
Sector(Sector <= 0) = 1;
LexOrd = [Layer Sector-1]*[1 nl]';

%% Compute SurfCov
SurfCov = length(unique(LexOrd))/nl/ns;

%% Compute volume estimate
if nargout > 1
    % Sort according to increasing lexicographic order
    [LexOrd,SortOrd] = sort(LexOrd);
    d = d(SortOrd);
    
    % Compute mean distance of the sector-layer intersections
    Dis = zeros(nl,ns); % mean distances
    np = size(P,1);     % number of points
    p = 1;
    while p <= np
        t = 1;
        while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
            t = t+1;
        end
        Dis(LexOrd(p)) = mean(d(p:p+t-1));
        p = p+t;
    end
    
    if nargout > 2
        % Interpolate missing distances
        D = Dis;
        dis = Dis;
        Dinv = D((nl:-1:1)',:);
        D = [Dinv Dinv Dinv; D D D; Dinv Dinv Dinv];
        Zero = Dis == 0;
        RadMean = mean(Dis(Dis > 0));
        for i = 1:nl
            for j = 1:ns
                if Zero(i,j)
                    if nnz(D(i+nl-1:i+nl+1,j+ns-1:j+ns+1)) > 1
                        d = D(i+nl-1:i+nl+1,j+ns-1:j+ns+1);
                        dis(i,j) = mean(d(d > 0));
                    elseif nnz(D(i+nl-2:i+nl+2,j+ns-2:j+ns+2)) > 1
                        d = D(i+nl-2:i+nl+2,j+ns-2:j+ns+2);
                        dis(i,j) = mean(d(d > 0));
                    elseif nnz(D(i+nl-3:i+nl+3,j+ns-3:j+ns+3)) > 1
                        d = D(i+nl-3:i+nl+3,j+ns-3:j+ns+3);
                        dis(i,j) = mean(d(d > 0));
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
