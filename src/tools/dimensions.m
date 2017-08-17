function [D,dir] = dimensions(points,varargin)

% Calculates the box-dimensions and dimension estimates of the point set 
% "points". Returns also the corresponding direction vectors.


if nargin == 2
    P = varargin{1};
    points = P(points,:);
elseif nargin == 3
    P = varargin{1};
    Bal = varargin{2};
    I = vertcat(Bal{points});
    points = P(I,:);
end

if size(points,2) == 3
    X = cov(points);
    [U,S,~] = svd(X);
    
    dp1 = points*U(:,1);
    dp2 = points*U(:,2);
    dp3 = points*U(:,3);
    
    D = [max(dp1)-min(dp1) max(dp2)-min(dp2) max(dp3)-min(dp3) ...
        (S(1,1)-S(2,2))/S(1,1) (S(2,2)-S(3,3))/S(1,1) S(3,3)/S(1,1)];
    
    dir = [U(:,1)' U(:,2)' U(:,3)'];
else
    X = cov(points);
    [U,S,~] = svd(X);
    
    dp1 = points*U(:,1);
    dp2 = points*U(:,2);
    
    D = [max(dp1)-min(dp1) max(dp2)-min(dp2) ...
        (S(1,1)-S(2,2))/S(1,1) S(2,2)/S(1,1)];
    
    dir = [U(:,1)' U(:,2)'];
end
