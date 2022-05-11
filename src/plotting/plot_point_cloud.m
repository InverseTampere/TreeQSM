function plot_point_cloud(P,fig,ms,col)

% Plots the given point cloud.
%
%   PLOT_POINT_CLOUD(P,FIG,MS,col) plots point cloud P in figure FIG using 
%   marker size MS and point color COL (string). P is a 2- or 3-column matrix 
%   where the first, second and third column gives the X-, Y-, and 
%   Z-coordinates of the points. 
%
%   PLOT_POINT_CLOUD(P,FIG) plots point cloud P in figure FIG using 
%   marker size 3 and color blue ('b').
%
%   PLOT_POINT_CLOUD(P) plots point cloud P in figure 1 using marker size 3
%   and color blue ('b').

if nargin == 1
    fig = 1;
    ms = 3;
    col = 'b';
elseif nargin == 2
    ms = 3;
    col = 'b';
elseif nargin == 3
    col = 'b';
end
if ms == 0
   ms = 3; 
end

col = ['.',col];

figure(fig)
if size(P,2) == 3
    plot3(P(:,1),P(:,2),P(:,3),col,'Markersize',ms)
elseif size(P,2) == 2
    plot(P(:,1),P(:,2),col,'Markersize',ms)
end
axis equal
    