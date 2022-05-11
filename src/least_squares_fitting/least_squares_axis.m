
function cyl = least_squares_axis(P,Axis,Point0,Rad0,weight)

% ---------------------------------------------------------------------
% LEAST_SQUARES_AXIS.M   Least-squares cylinder axis fitting using 
%                         Gauss-Newton when radius and point are given
%
% Version 1.0
% Latest update     1 Oct 2021
%
% Copyright (C) 2017-2021 Pasi Raumonen
% ---------------------------------------------------------------------
% Input    
% P         3d point cloud
% Axis0     Initial axis estimate (1 x 3)
% Point0    Initial estimate of axis point (1 x 3)
% Rad0      Initial estimate of the cylinder radius
% weight    (Optional) Weights for each point
% 
% Output
% cyl       Structure array with the following fields
%   axis      Cylinder axis (optimized here)
%   radius    Radius of the cylinder (from the input)
%   start     Axis point (from the input)
%   mad       Mean absolute distance of the points to the cylinder surface
%   SurfCov   Surface coverage, how much of the cylinder surface is covered 
%               with points
%   conv      If conv = 1, the algorithm has converged 
%   rel       If rel = 1, the algorithm has reliable answer in terms of
%               matrix inversion with a good enough condition number
% ---------------------------------------------------------------------


%% Initial estimates and other settings
res = 0.03; % "Resolution level" for computing surface coverage
par = [0 0]';
maxiter = 50; % maximum number of Gauss-Newton iteration
iter = 0; % number of iterations so far
conv = false; % converge of Gauss-Newton algorithm
rel = true; % are the results reliable, system matrix not badly conditioned
if nargin == 4
    weight = ones(size(P,1),1);
end
Rot0 = rotate_to_z_axis(Axis);
Pt = (P-Point0)*Rot0';

Par = [0 0 0 0 Rad0]';

%% Gauss-Newton iterations
while iter < maxiter && ~conv && rel
    
    % Calculate the distances and Jacobian
    [dist,J] = func_grad_axis(Pt,Par);
    
    % Calculate update step and gradient.
    SS0 = norm(dist); % Squared sum of the distances
    % solve for the system of equations: 
    % par(i+1) = par(i) - (J'J)^(-1)*J'd(par(i))
    A = J'*J;
    b = J'*dist;
    warning off
    p = -A\b; % solve for the system of equations
    warning on
    
    % Update
    par = par+p;
    
    % Check if the updated parameters lower the squared sum value
    Par = [0; 0; par; Rad0];
    dist = func_grad_axis(Pt,Par);
    SS1 = norm(dist);
    if SS1 > SS0
        % Update did not decreased the squared sum, use update with much
        % shorter update step
        par = par-0.95*p;
        Par = [0; 0; par; Rad0];
        dist = func_grad_axis(Pt,Par);
        SS1 = norm(dist);
    end
    
    % Check reliability
    rel = true;
    if rcond(A) < 10000*eps
        rel = false;
    end
    
    % Check convergence
    if abs(SS0-SS1) < 1e-5
        conv = true;
    end
    
    iter = iter+1;
end

%% Output
% Inverse transformation to find axis and point on axis 
% corresponding to original data
Rot = form_rotation_matrices(par);
Axis = Rot0'*Rot'*[0 0 1]'; % axis direction

% Compute the point distances to the axis
[dist,~,h] = distances_to_line(P,Axis,Point0); 
dist = dist-Rad0; % distances without weights
Len = max(h)-min(h);

% Compute mad (for points with maximum weights)
if nargin <= 4
  mad = mean(abs(dist)); % mean absolute distance to the circle
else
  I = weight == max(weight);
  mad = mean(abs(dist(I))); % mean absolute distance to the circle
end

% Compute SurfCov, minimum 3*8 grid
if ~any(isnan(par)) && rel && conv
  nl = ceil(Len/res);
  nl = max(nl,3);
  ns = ceil(2*pi*Rad0/res);
  ns = max(ns,8);
  ns = min(36,ns);
  SurfCov = single(surface_coverage(P,Axis,Point0,nl,ns,0.8*Rad0));
else
  SurfCov = single(0);
end


%% Define the output 
clear cir
cyl.radius = Rad0;
cyl.start = Point0;
cyl.axis = Axis';
cyl.mad = mad;
cyl.SurfCov = SurfCov;
cyl.conv = conv;
cyl.rel = rel;
