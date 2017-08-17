function plot_tree_structure2(P,Bal,Segs,SChi,fig,ms,BO,segind)

% Plots the branch-segmented tree point cloud so that each branching order
% has its own color Blue = trunk, green = 1st-order branches, 
% red = 2nd-order branches, etc.
%
% Inputs
% P         Point cloud
% Bal       Cover sets, Bal = cover.bal
% Segs      Segments, Segs = segment.segments
% SChi      Child segments, SChi = segment.ChildSegment
% fig       Figure number
% ms        Marker size
% BO        How many branching orders are plotted. 0 = all orders
% segind    Index of the segment where the plotting of tree structure
%                   starts. If segnum = 1 and BO = 0, then plots the whole
%                   tree. If segnum = 1 and B0 = 2, then plots the stem and
%                   the 1st-order branches. If segnum = 2 and BO = 0, then 
%                   plots the branch whose index is 2 and all its sub-branches. 


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

if BO == 0
    BO = 1000;
end

S = vertcat(Bal{Seg{segind}});
figure(fig)
plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(1,:),'Markersize',ms)
axis equal
forb = S;
if BO > 1
    %pause
    hold on
    c = SChi{segind};
    i = 2;
    while (i <= BO) && (~isempty(c))
        C = vertcat(Bal{unique(vertcat(Seg{c}))});
        C = setdiff(C,forb);
        figure(fig)
        plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
        axis equal
        c = unique(vertcat(SChi{c}));
        i = i+1;
        forb = union(forb,C);
        if i <= BO
            %pause
        end
    end
    hold off
end
