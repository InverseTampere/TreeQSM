function plot_segments(P,Bal,fig,ms,seg1,seg2,seg3,seg4,seg5)

% Plots point cloud segments/subsets defined as subsets of cover sets.
% If the subsets intersect, then assiggnes the common points to the
% segments given first.
%
% Inputs
% P     Point cloud
% Bal   Cover sets. Bal = cover.ball
% fig   figure number
% seg1  Segment/subset 1, color blue
% seg2  (Optional) Segment/subset 2, color red
% seg3  (Optional) Segment/subset 3, color green
% seg4  (Optional) Segment/subset 4, color cyan
% seg5  (Optional) Segment/subset 5, color magenta


if nargin == 5
    S1 = unique(vertcat(Bal{seg1}));
    figure(fig)
    plot3(P(S1,1),P(S1,2),P(S1,3),'b.','Markersize',ms)
    axis equal
elseif nargin == 6
    S1 = unique(vertcat(Bal{seg1}));
    S2 = unique(vertcat(Bal{seg2}));
    S2 = setdiff(S2,S1);
    figure(fig)
    plot3(P(S1,1),P(S1,2),P(S1,3),'b.','Markersize',1.5*ms)
    hold on
    plot3(P(S2,1),P(S2,2),P(S2,3),'r.','Markersize',ms)
    axis equal
    hold off
elseif nargin == 7
    S1 = unique(vertcat(Bal{seg1}));
    S2 = unique(vertcat(Bal{seg2}));
    S3 = unique(vertcat(Bal{seg3}));
    S2 = setdiff(S2,S1);
    S3 = setdiff(S3,S1);
    S3 = setdiff(S3,S2);
    figure(fig)
    plot3(P(S1,1),P(S1,2),P(S1,3),'b.','Markersize',ms)
    hold on
    plot3(P(S2,1),P(S2,2),P(S2,3),'r.','Markersize',ms)
    plot3(P(S3,1),P(S3,2),P(S3,3),'g.','Markersize',ms)
    axis equal
    hold off
elseif nargin == 8
    S1 = unique(vertcat(Bal{seg1}));
    S2 = unique(vertcat(Bal{seg2}));
    S3 = unique(vertcat(Bal{seg3}));
    S4 = unique(vertcat(Bal{seg4}));
    S2 = setdiff(S2,S1);
    S3 = setdiff(S3,S1);
    S3 = setdiff(S3,S2);
    S4 = setdiff(S4,S1);
    S4 = setdiff(S4,S2);
    S4 = setdiff(S4,S3);
    figure(fig)
    plot3(P(S1,1),P(S1,2),P(S1,3),'b.','Markersize',ms)
    hold on
    plot3(P(S2,1),P(S2,2),P(S2,3),'r.','Markersize',ms)
    plot3(P(S3,1),P(S3,2),P(S3,3),'g.','Markersize',ms)
    plot3(P(S4,1),P(S4,2),P(S4,3),'c.','Markersize',ms)
    axis equal
    hold off
elseif nargin == 9
    S1 = unique(vertcat(Bal{seg1}));
    S2 = unique(vertcat(Bal{seg2}));
    S3 = unique(vertcat(Bal{seg3}));
    S4 = unique(vertcat(Bal{seg4}));
    S5 = unique(vertcat(Bal{seg5}));
    S2 = setdiff(S2,S1);
    S3 = setdiff(S3,S1);
    S3 = setdiff(S3,S2);
    S4 = setdiff(S4,S1);
    S4 = setdiff(S4,S2);
    S4 = setdiff(S4,S3);
    S5 = setdiff(S5,S1);
    S5 = setdiff(S5,S2);
    S5 = setdiff(S5,S3);
    S5 = setdiff(S5,S4);
    figure(fig)
    plot3(P(S1,1),P(S1,2),P(S1,3),'b.','Markersize',ms)
    hold on
    plot3(P(S2,1),P(S2,2),P(S2,3),'r.','Markersize',ms)
    plot3(P(S3,1),P(S3,2),P(S3,3),'g.','Markersize',ms)
    plot3(P(S4,1),P(S4,2),P(S4,3),'c.','Markersize',ms)
    plot3(P(S5,1),P(S5,2),P(S5,3),'m.','Markersize',ms)
    axis equal
    hold off
end