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

function [dist,J] = func_grad_axis(P,par,weight)

% ---------------------------------------------------------------------
% FUNC_GRAD_CYLINDER.M   Function and gradient calculation for 
%                least-squares cylinder fit.
%
% Version 2.1.0
% Latest update     14 July 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Input 
% par       Cylinder parameters [x0 y0 alpha beta r]'
% P         Point cloud
% weight    (Optional) Weights for the points
% 
% Output
% dist      Signed distances of points to the cylinder surface:
%               dist(i) = sqrt(xh(i)^2 + yh(i)^2) - r, where 
%               [xh yh zh]' = Ry(beta) * Rx(alpha) * ([x y z]' - [x0 y0 0]')
% J         Jacobian matrix d dist(i)/d par(j).

% Changes from version 2.0.0 to 2.1.0, 14 July 2020:
% 1) Added optional input for weights of the points


% Five cylinder parameters: 
% Location = axis point intersects xy-plane: x0 and y0 (z0 == 0)
% Rotation angles around x and y axis = alpha and beta
% Radius = r
%
% Transformed points:
% Pt = [xh yx zh] = Ry(beta) * Rx(alpha) * (P - [x0 y0 0])
%
% "Plane points":
% Qt = Pt * I2 = [xh yh];
%
% Distance:
% D(x0,y0,alpha,beta,r) = sqrt( dot(Qt,Qt) )-r = sqrt( Qt*Qt' )-r
%
% Least squares = minimize Sum( D^2 ) over x0, y0, alpha, beta and r
%
% rt = sqrt( dot(Qt,Qt) )
% N = Qt/rt
%
% Jacobian for D with respect to x0, y0, alpha, beta:
% dD/di = dot( N,dQt/di ) = dot( Qt/rt,dQt/di )
%
% R = Ry*Rx
% dQt/dx0 = R*[-1 0 0]'
% dQt/dy0 = R*[0 -1 0]'
% dQt/dalpha = (P-[x0 y0 0])*DRx';
% dQt/dalpha = (P-[x0 y0 0])*DRy';
  
x0 = par(1);
y0 = par(2);
alpha = par(3);
beta = par(4);
r = par(5);

% Determine the rotation matrices and their derivatives
[R,DR1,DR2] = form_rotation_matrices([alpha beta]);

% Calculate the distances
Pt = (P-[x0 y0 0])*R';
xt = Pt(:,1);
yt = Pt(:,2);
rt = sqrt(xt.*xt + yt.*yt);
dist = rt-r; % Distances to the cylinder surface
if nargin == 3
  dist = weight.*dist; % Weighted distances
end

% form the Jacobian matrix
if nargout > 1
  
  N = [xt./rt yt./rt];
  m = size(P,1);
  J = zeros(m,2);
  
  A3 = (P-[x0 y0 0])*DR1';
  J(:,1) = sum(N(:,1:2).*A3(:,1:2),2);
  
  A4 = (P-[x0 y0 0])*DR2';
  J(:,2) = sum(N(:,1:2).*A4(:,1:2),2);
  
  if nargin == 3
    % Weighted Jacobian
    J = [weight.*J(:,1) weight.*J(:,2)];
  end
end
