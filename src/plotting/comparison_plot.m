function comparison_plot(P,R,fig,ms)

% Plots two point clouds "P" and "R" so that those points of "R" which are
% not in "P" are plotted in red whereas the common points are plotted in
% blue.
% "fig" and "ms" are the figure number and marker size.

if nargin == 3
    ms = 3;
elseif ms == 0
    ms = 3;
end

R = setdiff(R,P,'rows');

figure(fig)
plot3(P(:,1),P(:,2),P(:,3),'.b','Markersize',ms)
hold on
plot3(R(:,1),R(:,2),R(:,3),'.r','Markersize',ms)
hold off
axis equal
%set(gca,'Color',[0 0 0]);