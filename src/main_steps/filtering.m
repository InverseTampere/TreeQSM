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

function Pass = filtering(P,r1,n1,d2,r2,n2,Scaling,AllPoints)

% ---------------------------------------------------------------------
% FILTERING.M       Filters noise from point clouds.
%
% Version 2.00
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Performs an initial filtering of the point cloud.
% 
% At first the possible NaNs are removed.
% 
% Secondly those points which belong to r1-balls that have at least
% n1 points them are included. Comprehensive filtering here means that
% the r1-ball is defined for every point whereas non-comprehensive means
% that only a cover of point cloud with r1-balls is defined
%
% Lastly small components (less than n2 r2-balls) are removed.
%
% For both filtering procedures the given minimum number (n1 and n2) 
% applies for the first three meters above the ground and after that it is
% modified for each meter by the change in the mean point density. That is,
% if the point density decreases, which usually is the case, then the
% minimum number is also decreased accordingly.
%
% Inputs:
% P         Point cloud
% r1        Radius of the balls used in the first filtering
% n1        Minimum number of points in the accepted balls of the first filtering
% r2        Radius of the balls used in the second filtering
% n2        Minimum number of balls in the components passing the second filtering
% Scaling   If true, the first filtering threshold "n1" is scaled along the 
%               height with average point density
% Comp      If true, does the first filtering process for every point
%
% Outputs:
% Pass      Logical vector indicating points passing the filtering


% Default firts filtering is not comprehensive and with scaling
if nargin == 6
    Scaling = false;
    AllPoints = false; 
elseif nargin == 7
    AllPoints = false;
end

% Only double precision data
if ~isa(P,'double')
    P = double(P);
end
% Only x,y,z-data
if size(P,2) > 3
    P = P(:,1:3);
end
np = size(P,1);
np0 = np;

% Remove possible NaNs
I = any(isnan(P),2);
if any(I)
    P = P(~I,:);
end


%% Partition the point cloud into cubes
[partition,CC] = cubical_partition(P,r1);

if ~AllPoints
    % Do the first filtering with a cover, not comprehensive filtering
    %% Generate a cover and determine the largest (number of points) ball for each point
    NotInspected = true(np,1);
    NumOfPoints = zeros(np,1);
    r1 = r1^2;
    for i = 1:np
        if NotInspected(i)
            points = partition(CC(i,1)-1:CC(i,1)+1,CC(i,2)-1:CC(i,2)+1,CC(i,3)-1:CC(i,3)+1);
            points = vertcat(points{:,:});
            cube = P(points,:);
            dist = (P(i,1)-cube(:,1)).^2+(P(i,2)-cube(:,2)).^2+(P(i,3)-cube(:,3)).^2;
            J = dist < r1;
            I = points(J);
            NotInspected(I) = false;
            N = NumOfPoints(I);
            n = length(I);
            K = N < n;
            NumOfPoints(I(K)) = n;
        end
    end
else
    % Do the first filtering comprehensively by defining the neighborhoods
    % for all points
    
    %% Generate the balls and determine the number of points for each point
    NumOfPoints = zeros(np,1);
    r1 = r1^2;
    for i = 1:np
        points = partition(CC(i,1)-1:CC(i,1)+1,CC(i,2)-1:CC(i,2)+1,CC(i,3)-1:CC(i,3)+1);
        points = vertcat(points{:,:});
        cube = P(points,:);
        dist = (P(i,1)-cube(:,1)).^2+(P(i,2)-cube(:,2)).^2+(P(i,3)-cube(:,3)).^2;
        J = dist < r1;
        NumOfPoints(i) = nnz(J);
    end
end
clearvars partition CC


%% First filtering
if Scaling
    % Use smaller treshold for upper parts of the tree
    % Change the treshold "n1" according the average points in the cover
    % sets for every meter
    hmin = min(P(:,3));
    hmax = max(P(:,3));
    H = ceil(hmax-hmin);
    D = zeros(H,1);
    J = false(np,1);
    Pass = false(np,1);
    A = (1:1:np)';
    for i = 1:H
        I = P(:,3) < hmin+i;
        K = I&~J;
        J = I;
        D(i) = ceil(mean(NumOfPoints(K)));
        if i <= 2
            M = A(K);
            N = NumOfPoints(K) >= n1;
            M = M(N);
            Pass(M) = true;
        else
            M = A(K);
            m = max(ceil(n1*D(i)/D(2)),ceil(n1/3));
            N = NumOfPoints(K) >= m;
            M = M(N);
            Pass(M) = true;
        end
    end
else
    Pass = NumOfPoints >= n1;
end
clearvars NumOfPoints

% Display filtering results
np = nnz(Pass);
nf = np0-np;
str = ['    All points: ',num2str(np0),', First filtering: ',num2str(nf),', Points left: ',num2str(np)];
disp(str)


%% Cover the point cloud with r2-balls for the second filtering
clear inputs
inputs.BallRad1 = r2;
inputs.PatchDiam1 = d2;
inputs.nmin1 = 0;
cover = cover_sets(P(Pass,1:3),inputs);


%% Determine the separate components
Components = connected_components(cover.neighbor,0,n2);


%% Second filtering
B = vertcat(Components{:});
I = vertcat(cover.ball{B});
J = false(np,1);
J(I) = true;
I = (1:1:np0)';
I = I(Pass);
I = I(~J);
Pass(I) = false;

% Display filtering results
nf = np-nnz(Pass);
npl = nnz(Pass);
str = ['    All points: ',num2str(np),', Second filtering: ',num2str(nf),', Points left: ',num2str(npl)];
disp(str)

nf = np0-npl;
str = ['    All points: ',num2str(np0),', All filtered points: ',num2str(nf),', Points left: ',num2str(npl)];
disp(str)


%% Plot the filtered and unfiltered point clouds
comparison_plot(P(Pass,:),P,1,5)
point_cloud_plotting(P(Pass,:),2,5)
