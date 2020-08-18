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

function cyl = least_squares_cylinder(P,cyl0,weight,Q)
% ---------------------------------------------------------------------
% LEAST_SQUARES_CYLINDER.M   Least-squares cylinder using Gauss-Newton.
%
% Version 1.3.0
% Latest update     14 July 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------
% Input    
% P         Point cloud
% cyl0      Initial estimates of the cylinder parameters
% weight    (Optional) Weights of the points for fitting
% Q         (Optional) Subset of "P" where the cylinder is intended
% 
% Output  
% cyl       Structure array containing the following fields:
%   radius      Radius of the cylinder
%   length      Length of the cylinder
%   start       Point on the axis at the bottom of the cylinder (1 x 3)
%   axis        Axis direction of the cylinder (1 x 3) 
%   mad         Mean absolute distance between points and cylinder surface
%   SurfCov     Relative cover of the cylinder's surface by the points 
%   SurfCovDis  Mean point distances to the axis from the layer-sectors
%   SurfCovVol  Estimate of cylinder's volume based on radii given by "SurfCovDis"
%   dist        Radial distances from the points to the cylinder (m x 1) 
%   conv        If conv = 1, the algorithm has converged 
%   rel         If rel = 1, the algorithm has reliable answer in terms of
%                   matrix inversion with a good enough condition number
% ---------------------------------------------------------------------

% Changes from version 1.2.0 to 1.3.0, 14 July 2020:  
% 1) Changed the input parameters of the cylinder to the struct format.
% 2) Added optional input for weights
% 3) Added optional input "Q", a subset of "P", the cylinder is intended
%    to be fitted in this subset but it is fitted to "P" to get better
%    estimate of the axis direction and radius

% Changes from version 1.1.0 to 1.2.0, 14 Jan 2020:  
% 1) Changed the outputs and optionally the inputs to the struct format.
% 2) Added new output, "mad", which is the mean absolute distance of the
%    points from the surface of the cylinder.
% 3) Added new output, "SurfCov", that measures how well the surface of the
%    cylinder is covered by the points.
% 4) Added new output, "SurfCovDis", which is a matrix of mean point distances 
%    from layer/sector-intersections to the axis.
% 5) Added new output, "SurfCovVol", which is an estimate of the cylinder's 
%    volume based on the radii in "SurfCovDis" and "cylindrical sectors".
% 6) Added new optional input "res" which gives the point resolution level
%    for computing SurfCov: the width and length of sectors/layers.

% Changes from version 1.0.0 to 1.1.0, 3 Oct 2019:  
% 1) Bug fix: --> "Point = Rot0'*([par(1) par(2) 0]')..."

% "Resolution level" for computing surface coverage
res = 0.03;

% Transform the data to close to standard position via a rotation 
% followed by a translation 
Rot0 = rotate_to_z_axis(cyl0.axis);
Point1 = Rot0*cyl0.start';
Pt = P*Rot0'-Point1';

% Initial estimates and tolerance information
par = [0 0 0 0 cyl0.radius]'; 

% Gauss-Newton algorithm to find estimate of rotation-translation 
% parameters that transform the data so that the best-fit circle is 
% one in standard position
if nargin == 2
    [par,dist,conv,rel] = nlssolver(par,Pt);
else
    [par,dist,conv,rel] = nlssolver(par,Pt,weight);
end
cyl.radius = single(par(5)); % radius

% Inverse transformation to find axis and point on axis 
% corresponding to original data
Rot = form_rotation_matrices(par(3:4));
Axis = Rot0'*Rot'*[0 0 1]'; % axis direction
Point = Rot0'*([par(1) par(2) 0]')+cyl0.start'; % axis point

% Compute the start, length and mad, translate the axis point to the 
% cylinder's bottom:
% If the fifth input (point cloud Q) is given, use it for the start, length, 
% mad, and SurfCov
if nargin == 4
    P = Q;
end
H = P*Axis;
hmin = min(H);
cyl.length = single(abs(max(H)-hmin));
hpoint = Axis'*Point;
Point = Point-(hpoint-hmin)*Axis; % axis point at the cylinder's bottom
cyl.start = single(Point');
cyl.axis = single(Axis');
cyl.mad = single(average(abs(dist))); % mean absolute distance
%cyl.dist = single(dist);
cyl.conv = conv;
cyl.rel = rel;

% Compute SurfCov, minimum 3*8 grid
if ~any(isnan(Axis)) && ~any(isnan(Point)) && rel && conv
    nl = ceil(cyl.length/res);
    nl = max(nl,3);
    ns = ceil(2*pi*cyl.radius/res);
    ns = max(ns,8);
    ns = min(36,ns);
    SurfCov = surface_coverage(P,Axis,Point,nl,ns);
    %[SurfCov,Dis,CylVol] = surface_coverage(P,Axis,Point,nl,ns);
    
    cyl.SurfCov = single(SurfCov);
    %cyl.SurfCovVol = single(CylVol);
    %cyl.SurfCovDis = single(Dis);
else
    cyl.SurfCov = single(0);
    %cyl.SurfCovVol = single(0);
    %cyl.SurfCovDis = single(0);
end
