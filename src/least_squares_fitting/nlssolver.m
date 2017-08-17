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

function [par,d,conv,rel] = nlssolver(par,P)

% ---------------------------------------------------------------------
% NLSSOLVER.M   Nonlinear least squares solver for cylinders.
%
% Version 2.0
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
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


maxiter = 50;
iter = 0;
conv = false;
rel = true;

% Gauss-Newton iterations
while iter < maxiter && ~conv && rel
    
    % Calculate the distances and Jacobian
    [d0, J] = func_grad_cylinder(par,P);
    
    % Calculate update step and gradient.
    SS0 = norm(d0); % Squared sum of the distances
    Ra = triu(qr([J, d0]));
    R = Ra(1:5,1:5);
    q = Ra(1:5,6);
    warning off
    p = -R\q; % solve for the system of equations
    warning on
    par = par+p; % update the parameters
    
    % Check reliability
    if rcond(-R) < 10000*eps
        rel = false;
    end
    
    % Check convergence
    d = func_grad_cylinder(par,P); % The distances with the new parameter values
    SS1 = norm(d); % Squared sum of the distances
    if abs(SS0-SS1) < 1e-4
        conv = true;
    end
    
    iter = iter + 1;
end
