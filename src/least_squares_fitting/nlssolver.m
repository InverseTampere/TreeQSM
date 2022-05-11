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

function [par,d,conv,rel] = nlssolver(par,P,weight)

% ---------------------------------------------------------------------
% NLSSOLVER.M   Nonlinear least squares solver for cylinders.
%
% Version 2.1.0
% Latest update     14 July 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Input 
% par       Intial estimates of the parameters
% P         Point cloud
%           
% Output 
% par       Optimised parameter values
% d         Distances of points to cylinder
% conv      True if fitting converged
% rel       True if condition number was not very bad, fit was reliable

% Changes from version 2.0.0 to 2.1.0, 14 July 2020:
% 1) Added optional input for weights of the points

maxiter = 50;
iter = 0;
conv = false;
rel = true;

if nargin == 2
    NoWeights = true;
else
    NoWeights = false;
end

%% Gauss-Newton iterations
while iter < maxiter && ~conv && rel
    
    %% Calculate the distances and Jacobian
    if NoWeights
        [d0, J] = func_grad_cylinder(par,P);
    else
        [d0, J] = func_grad_cylinder(par,P,weight);
    end
    
    %% Calculate update step
    SS0 = norm(d0); % Squared sum of the distances
    % solve for the system of equations: 
    % par(i+1) = par(i) - (J'J)^(-1)*J'd0(par(i))
    A = J'*J;
    b = J'*d0;
    warning off
    p = -A\b; % solve for the system of equations
    warning on
    par = par+p; % update the parameters
    
    %% Check reliability
    if rcond(-R) < 10000*eps
        rel = false;
    end
    
    %% Check convergence:
    % The distances with the new parameter values:
    if NoWeights
        d = func_grad_cylinder(par,P); 
    else
        d = func_grad_cylinder(par,P,weight); 
    end
    SS1 = norm(d); % Squared sum of the distances
    if abs(SS0-SS1) < 1e-4
        conv = true;
    end
    
    iter = iter + 1;
end
