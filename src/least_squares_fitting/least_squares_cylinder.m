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

function [Rad,Len,Point,Axis,dist,conv,rel] = least_squares_cylinder(P,Point0,Axis0,Rad0)
% ---------------------------------------------------------------------
% LEAST_SQUARES_CYLINDER.M   Least-squares cylinder using Gauss-Newton.
%
% Version 1.0
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------
% Input    
% P         Point cloud
% Point0    Initial estimate of the point on the axis (1 x 3) 
% Axis0     Initial estimate of the axis direction (1 x 3) 
% Rad0      Initial estimate of the cylinder radius
% 
% Output  
% Rad       Radius of the cylinder
% Len       Length of the cylinder
% Point     Point on the axis at the bottom of the cylinder (1 x 3)
% Axis      Axis direction of the cylinder (1 x 3) 
% dist      Radial distances from the points to the cylinder (m x 1) 
% conv      If conv = 1, the algorithm has converged 
% rel       If rel = 1, the algorithm has reliable answer in terms of
%               matrix inversion with a good enough condition number
% ---------------------------------------------------------------------


% Transform the data to close to standard position via a rotation 
% followed by a translation 
Rot0 = rotate_to_z_axis(Axis0);
Point1 = Rot0*Point0';
Pt = mat_vec_subtraction((P*Rot0'),Point1);

% Initial estimates and tolerance information
par = [0 0 0 0 Rad0]'; 

% Gauss-Newton algorithm to find estimate of roto-translation 
% parameters that transform the data so that the best-fit circle is 
% one in standard position
[par,dist,conv,rel] = nlssolver(par,Pt);

% Inverse transformation to find axis and point on axis 
% corresponding to original data 
Rad = par(5); % radius
Rot = form_rotation_matrices(par(3:4));
Axis = Rot0'*Rot'*[0 0 1]'; % axis direction
Point = Rot*([par(1) par(2) 0]')+Point0';

H = P*Axis;
hmin = min(H);
Len = abs(max(H)-hmin);
hpoint = Axis'*Point;
Point = Point-(hpoint-hmin)*Axis;
Point = Point';
Axis = Axis';
