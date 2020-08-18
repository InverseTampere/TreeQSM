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

function [dist,J] = func_grad_cylinder(par,P,weight)

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

m = size(P, 1);
x0 = par(1);
y0 = par(2);
alpha = par(3);
beta = par(4);
r = par(5);

% Determine the rotation matrices and their derivatives
[R,DR1,DR2] = form_rotation_matrices([alpha beta]);

% Calculate the distances
Pt = (P - ones(m,1)*[x0 y0 0])*R';
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
    J = zeros(m,5);
    
    A1 = (R*[-1 0 0]')';
    A = eye(2);
    A(1,1) = A1(1);  A(2,2) = A1(2);
    J(:,1) = sum(N(:,1:2)*A,2);
        
    A2 = (R*[0 -1 0]')';
    A(1,1) = A2(1);  A(2,2) = A2(2);
    J(:,2) = sum(N(:,1:2)*A,2);
    
    A3 = (P-repmat([x0 y0 0],[m 1]))*DR1';
    J(:,3) = sum(N(:,1:2).*A3(:,1:2),2);
        
    A4 = (P-repmat([x0 y0 0],[m 1]))*DR2';
    J(:,4) = sum(N(:,1:2).*A4(:,1:2),2);
    
    J(:,5) = -1*ones(m,1);
    if nargin == 3
        % Weighted Jacobian
        J = [weight.*J(:,1) weight.*J(:,2) weight.*J(:,3) ...
            weight.*J(:,4) weight.*J(:,5)];
    end
end
