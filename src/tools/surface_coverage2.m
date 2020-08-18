function SurfCov = surface_coverage2(Axis,Len,Vec,height,nl,ns)

% Computes surface coverage (number between 0 and 1) of points on cylinder 
% surface defined by "Axis" and "Len". "Vec" are the vectors connecting 
% points to the Axis and "height" are the heights of the points from 
% the base of the cylinder

[U,W] = orthonormal_vectors(Axis);
Vec = Vec*[U W];
ang = atan2(Vec(:,2),Vec(:,1))+pi;
I = ceil(height/Len*nl);
I(I == 0) = 1;
I(I > nl) = nl;
J = ceil(ang/2/pi*ns);
J(J == 0) = 1;
K = [I J-1]*[1 nl]';
SurfCov = length(unique(K))/nl/ns;