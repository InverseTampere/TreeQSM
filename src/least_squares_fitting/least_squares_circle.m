function cir = least_squares_circle(P,Point0,Rad0,weight)
% ---------------------------------------------------------------------
% LEAST_SQUARES_CIRCLE.M   Least-squares circle fitting using Gauss-Newton.
%
% Version 1.1.0
% Latest update     6 Oct 2021
%
% Copyright (C) 2017-2021 Pasi Raumonen
% ---------------------------------------------------------------------
% Input    
% P         2d point cloud
% Point0    Initial estimate of centre (1 x 2)
% Rad0      Initial estimate of the circle radius
% weight    Optional, weights for each point
% 
% Output  
% Rad       Radius of the cylinder
% Point     Centre point (1 x 2)
% ArcCov    Arc point coverage (%), how much of the circle arc is covered with points
% conv      If conv = 1, the algorithm has converged 
% rel       If rel = 1, the algorithm has reliable answer in terms of
%               matrix inversion with a good enough condition number
% ---------------------------------------------------------------------


%% Initial estimates and other settings
par = [Point0 Rad0]'; 
maxiter = 200; % maximum number of Gauss-Newton iteration
iter = 0; % number of iterations so far
conv = false; % converge of Gauss-Newton algorithm
rel = true; % are the reusults reliable in the sense that system matrix was not badly conditioned
if nargin == 3
  weight = ones(size(P,1),1);
end

%% Gauss-Newton iterations
while iter < maxiter && ~conv && rel
  
  % Calculate the distances and Jacobian
  [dist,J] = func_grad_circle(P,par,weight);
  
  % Calculate update step and gradient.
  SS0 = norm(dist); % Squared sum of the distances
  % solve for the system of equations: par(i+1) = par(i) - (J'J)^(-1)*J'd(par(i))
  A = J'*J;
  b = J'*dist;
  warning off
  p = -A\b; % solve for the system of equations
  warning on
  
  % Update
  par = par+p;
  
  % Check if the updated parameters lower the squared sum value
  dist = func_grad_circle(P,par,weight);
  SS1 = norm(dist);
  if SS1 > SS0
    % Update did not decreased the squared sum, use update with much
    % shorter update step
    par = par-0.95*p;
    dist = func_grad_circle(P,par,weight);
    SS1 = norm(dist);
  end
  
  % Check reliability
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
Rad = par(3);
Point = par(1:2);
U = P(:,1)-Point(1);
V = P(:,2)-Point(2);
dist = sqrt(U.*U+V.*V)-Rad;
if nargin <= 3
  mad = mean(abs(dist)); % mean absolute distance to the circle
else
  I = weight == max(weight);
  mad = mean(abs(dist(I))); % mean absolute distance to the circle
end
% Calculate ArcCov, how much of the circle arc is covered with points
if ~any(isnan(par))
  if nargin <= 3
    I = dist > -0.2*Rad;
  else
    I = dist > -0.2*Rad & weight == max(weight);
  end
  U = U(I,:);   V = V(I,:);
  ang = atan2(V,U)+pi;
  ang = ceil(ang/2/pi*100);
  ang(ang <= 0) = 1;
  Arc = false(100,1);
  Arc(ang) = true;
  ArcCov = nnz(Arc)/100;
else
  ArcCov = 0;
end

cir.radius = Rad;
cir.point = Point';
cir.mad = mad;
cir.ArcCov = ArcCov;
cir.conv = conv;
cir.rel = rel;