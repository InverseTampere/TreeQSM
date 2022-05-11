function cir = least_squares_circle_centre(P,Point0,Rad0)
% ---------------------------------------------------------------------
% LEAST_SQUARES_CIRCLE_CENTRE.M   Least-squares circle fitting such that
%                                   radius is given (fits the centre)
%
% Version 1.0.0
% Latest update     6 Oct 2021
%
% Copyright (C) 2017-2021 Pasi Raumonen
% ---------------------------------------------------------------------
% Input    
% P         2d point cloud
% Point0    Initial estimate of centre (1 x 2)
% Rad0      The circle radius
% weight    Optional, weights for each point
% 
% Output  
% cir     Structure array with the following fields
%   Rad       Radius of the cylinder
%   Point     Centre point (1 x 2)
%   ArcCov    Arc point coverage (%), how much of the circle arc is covered 
%               with points
%   conv      If conv = 1, the algorithm has converged 
%   rel       If rel = 1, the algorithm has reliable answer in terms of
%               matrix inversion with a good enough condition number
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 1.1.0, 6 Oct 2021:  
% 1) Streamlining code and some computations

%% Initial estimates and other settings
par = [Point0 Rad0]'; 
maxiter = 200; % maximum number of Gauss-Newton iteration
iter = 0; % number of iterations so far
conv = false; % converge of Gauss-Newton algorithm
rel = true; % the results reliable (system matrix was not badly conditioned)

%% Gauss-Newton iterations
while iter < maxiter && ~conv && rel
  
  % Calculate the distances and Jacobian
  [dist,J] = func_grad_circle_centre(P,par);
  
  % Calculate update step and gradient.
  SS0 = norm(dist); % Squared sum of the distances
  % solve for the system of equations: par(i+1) = par(i) - (J'J)^(-1)*J'd(par(i))
  A = J'*J;
  b = J'*dist;
  warning off
  p = -A\b; % solve for the system of equations
  warning on
  
  % Update
  par(1:2,1) = par(1:2,1)+p;
  
  % Check if the updated parameters lower the squared sum value
  dist = func_grad_circle_centre(P,par);
  SS1 = norm(dist);
  if SS1 > SS0
    % Update did not decreased the squared sum, use update with much
    % shorter update step
    par(1:2,1) = par(1:2,1)-0.95*p;
    dist = func_grad_circle_centre(P,par);
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
Point = par(1:2);
if conv && rel
  % Calculate ArcCov, how much of the circle arc is covered with points
  U = P(:,1)-par(1);
  V = P(:,2)-par(2);
  ang = atan2(V,U)+pi;
  I = false(100,1);
  ang = ceil(ang/2/pi*100);
  I(ang) = true;
  ArcCov = nnz(I)/100;
  % mean absolute distance to the circle
  d = sqrt(U.*U+V.*V)-Rad0;
  mad = mean(abs(d)); 
else
  mad = 0;
  ArcCov = 0;
end

cir.radius = Rad0;
cir.point = Point';
cir.mad = mad;
cir.ArcCov = ArcCov;
cir.conv = conv;
cir.rel = rel;