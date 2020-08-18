function plot_branch_segmentation(P,cover,segment,Color,fig,ms,segind,BO)

% ---------------------------------------------------------------------
% PLOT_BRANCH_SEGMENTATION.M   Plots branch-segmented point cloud, coloring
%                               based on branching order or branches
%
% Version 1.0.0
% Latest update     13 July 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------
% 
% If the coloring is based on branches (Color = 'branch'), then each segment 
% is colored with unique color. If the coloring is based on branching order 
% (Color = 'order'), then Blue = trunk, Green = 1st-order branches, 
% Red = 2nd-order branches, etc.
% 
% If segind = 1 and BO = 0, then plots the stem. If segind = 1 and BO = 1, 
% then plots the stem and the 1st-order branches. If segind = 1 and 
% BO >= maximum branching order or BO input is not given, then plots the 
% whole tree. If segind = 2 and BO is not given or it is high enough, then
% plots the branch whose index is 2 and all its sub-branches. 
%
% Inputs
% P         Point cloud
% cover     Cover sets structure
% segment   Segments structure
% Color     Color option, 'order' or 'branch'
% fig       Figure number
% ms        Marker size
% segind    Index of the segment where the plotting of tree structure starts. 
% BO        How many branching orders are plotted. 0 = stem, 1 = 1st order, etc

n = nargin;
if n < 8
    BO = 1000;
    if n < 7
        segind = 1;
        if n < 6
            ms = 1;
            if n < 5
                fig = 1;
                if n == 3
                    Color = 'order';
                end
            end
        end
    end
end

Bal = cover.ball;
Segs = segment.segments;
SChi = segment.ChildSegment;
SPar = segment.ParentSegment;
ns = max(size(Segs));

if iscell(Segs{1})
    Seg = cell(ns,1);
    for i = 1:ns
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

if strcmp(Color,'branch')
    Color = 1;
    % Color the segments with unique colors
    col = rand(ns,3);
    for i = 2:ns
        C = col(SPar(i),:);
        c = col(i,:);
        while sum(abs(C-c)) < 0.2
            c = rand(1,3);
        end
        col(i,:) = c;
    end
elseif strcmp(Color,'order')
    Color = 0;
    % Color the cylinders in branches based on the branch order
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
    col = repmat(col,[10,1]);
end

segments = segind;
C = SChi{segind};
b = 1;
order = 1;
while ~isempty(C) && b <= BO
    b = b+1;
    segments = [segments; C];
    order = [order; b*ones(length(C),1)];
    C = vertcat(SChi{C});
end

ns = length(segments);
figure(fig)
for i = 1:ns
    if i == 2
        hold on
    end
    S = vertcat(Bal{Seg{segments(i)}});
    if Color
    % Coloring based on branch 
    plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(segments(i),:),'Markersize',ms)
    else
    % Coloring based on branch order
    plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(order(i),:),'Markersize',ms)    
    end
end
hold off
axis equal
