function [dist,J] = func_grad_circle(P,par,weight)

% ---------------------------------------------------------------------
% FUNC_GRAD_CIRCLE.M   Function and gradient calculation for 
%                           least-squares circle fit.
%
% Version 1.0
% Latest update     20 Oct 2017
%
% Copyright (C) 2017 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Input 
% P         Point cloud
% par       Circle parameters [x0 y0 r]'
% weight    Weights for the points. Weight the distances.
% 
% Output
% dist      Signed distances of points to the circle:
%               dist(i) = sqrt((xi-x0)^2 + (yi-y0)^2) - r, where 
%               
% J         Jacobian matrix d dist(i)/d par(j).


% Calculate the distances
Vx = P(:,1)-par(1);
Vy = P(:,2)-par(2);
rt = sqrt(Vx.*Vx + Vy.*Vy);
if nargin == 3
    dist = weight.*(rt-par(3)); % Weighted distances to the circle
else
    dist = rt-par(3); % Distances to the circle
end

% form the Jacobian matrix
if nargout > 1
    m = size(P, 1);
    J = zeros(m,3);
    J(:,1) = -Vx./rt;
    J(:,2) = -Vy./rt;
    J(:,3) = -1*ones(m,1);
    % apply the weights
    if nargin == 3
        J = [weight.*J(:,1) weight.*J(:,2) weight.*J(:,3)];
    end
end
