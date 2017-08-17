function point_cloud_plotting(P,fig,ms,Bal,Sub)

% Plots the given point cloud "P". With additional inputs one can plot only
% those points that are included in the cover sets "Bal" or in the
% subcollection "Sub" of the cover sets. 
% "fig" and "ms" are the figure number and marker size.

if nargin == 2
    ms = 3;
elseif ms == 0
    ms = 3;
end

if nargin < 4
    figure(fig)
    if size(P,2) == 3
        plot3(P(:,1),P(:,2),P(:,3),'.b','Markersize',ms)
    elseif size(P,2) == 2
        plot(P(:,1),P(:,2),'.b','Markersize',ms)
    end
    axis equal
    
elseif nargin == 4
    I = vertcat(Bal{:});
    figure(fig)
    plot3(P(I,1),P(I,2),P(I,3),'.b','Markersize',ms)
    axis equal
    
else
    if iscell(Sub)
        S = vertcat(Sub{:});
        Sub = vertcat(S{:});
        I = vertcat(Bal{Sub});
        figure(fig)
        plot3(P(I,1),P(I,2),P(I,3),'.b','Markersize',ms)
        axis equal
    else
        I = vertcat(Bal{Sub});
        figure(fig)
        plot3(P(I,1),P(I,2),P(I,3),'.b','Markersize',ms)
        axis equal
    end
end