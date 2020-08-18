function plot_scatter(P,C,fig,ms)

% A scatter plot where the color of each 2d or 3d point is specified by a
% number. 
%
% Inputs:
% P     point cloud
% C     color data (vector, value for each point in P)
% fig   figure number
% ms    marker size

S = normalize(C);
figure(fig)
if size(P,2) == 3
    scatter3(P(:,1),P(:,2),P(:,3),ms*ones(size(P,1),1),C,'filled')
    %scatter3(P(:,1),P(:,2),P(:,3),ms*S.*ones(size(P,1),1),C,'filled')
elseif size(P,2) == 2
    scatter(P(:,1),P(:,2),ms*S.*ones(size(P,1),1),C,'filled')
end
axis equal
if size(C,2) == 1
    colormap(jet(25))
    %caxis([0 max(C)])
    if min(C) < max(C)
    caxis([min(C) max(C)])
    end
    colorbar
end
