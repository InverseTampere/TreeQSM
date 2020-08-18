function plot_comparison(P1,P2,fig,ms1,ms2)

% Plots two point clouds "P1" and "P2" so that those points of "P2" which are
% not in "P1" are plotted in red whereas the common points are plotted in
% blue. "fig" and "ms1" and "ms2" are the figure number and marker sizes.

if nargin == 3
    ms1 = 3;
    ms2 = 3;
elseif nargin == 4
    ms2 = 3;
end

if ms1 == 0
    ms1 = 3;
end
if ms2 == 0
    ms2 = 3;
end

P2 = setdiff(P2,P1,'rows');

figure(fig)
if size(P1,2) == 3
    plot3(P1(:,1),P1(:,2),P1(:,3),'.b','Markersize',ms1)
    hold on
    plot3(P2(:,1),P2(:,2),P2(:,3),'.r','Markersize',ms2)
elseif size(P1,2) == 2
    plot(P1(:,1),P1(:,2),'.b','Markersize',ms1)
    hold on
    plot(P2(:,1),P2(:,2),'.r','Markersize',ms2)
end
hold off
axis equal