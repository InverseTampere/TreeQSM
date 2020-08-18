function plot_tree_structure(P,cover,segment,fig,ms,segind,BO)

% ---------------------------------------------------------------------
% PLOT_TREE_STRUCTURE.M       Plots branch-segmented point cloud with unique
%                               color for each branching order
%
% Version 1.1.0
% Latest update     13 July 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------
% 
% Blue = trunk, Green = 1st-order branches, Red = 2nd-order branches, etc.
% If segind = 1 and BO = 0, then plots the stem. If segind = 1 and BO = 1, 
% then plots the stem and the 1st-order branches. If segind = 1 and 
% BO >= maximum branching order or BO input is not given, then plots the 
% whole tree. If segind = 2 and BO is not given or it is high enough, then
% plots the branch whose index is 2 and all its sub-branches. 
%
% Inputs
% P         Point cloud
% cover     Cover sets structure
% Segs      Segments structure
% fig       Figure number
% ms        Marker size
% segind    Index of the segment where the plotting of tree structure
%                   starts. 
% BO        How many branching orders are plotted. 0 = stem, 1 = 1st order, etc
% 

% Changes from version 1.0.0 to 1.1.0, 13 July 2020:
% 1) Added option for choosing the coloring based either on branch order or
%    unique color for each branch

n = nargin;
if n < 7
    BO = 1000;
    if n < 6
        segind = 1;
        if n < 5
            ms = 1;
            if n == 3
                fig = 1;
            end
        end
    end
end

Bal = cover.ball;
Segs = segment.segments;
SChi = segment.ChildSegment;

col = [
	0.00  0.00  1.00
	0.00  0.50  0.00
	1.00  0.00  0.00
	0.00  0.75  0.75
	0.75  0.00  0.75
	0.75  0.75  0.00
	0.25  0.25  0.25
	0.75  0.25  0.25
	0.95  0.95  0.00
	0.25  0.25  0.75
	0.75  0.75  0.75
	0.00  1.00  0.00
	0.76  0.57  0.17
	0.54  0.63  0.22
	0.34  0.57  0.92
	1.00  0.10  0.60
	0.88  0.75  0.73
	0.10  0.49  0.47
	0.66  0.34  0.65
	0.99  0.41  0.23];
col = repmat(col,[1000,1]);

if iscell(Segs{1})
    n = max(size(Segs));
    Seg = cell(n,1);
    for i = 1:n
        m = size(Segs{i},1);
        S = zeros(0);
        for j = 1:m
            s = Segs{i}(j);
            s = s{:};
            S = [S; s];
        end
        Seg{i} = S;
    end
else
    Seg = Segs;
end

S = vertcat(Bal{Seg{segind}});
figure(fig)
plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(1,:),'Markersize',ms)
axis equal
%forb = S;
if BO > 0
    hold on
    c = SChi{segind};
    order = 1;
    while (order <= BO) && (~isempty(c))
        C = vertcat(Bal{vertcat(Seg{c})});
        %C = setdiff(C,forb);
        figure(fig)
        plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(order+1,:),'Markersize',ms)
        axis equal
        c = unique(vertcat(SChi{c}));
        order = order+1;
        %forb = union(forb,C);
    end
    hold off
end
