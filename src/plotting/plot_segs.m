function plot_segs(P,comps,fig,ms,Bal)

% Plots the point cloud segments given in the cell array "comps".
% If 4 inputs, cells contain the point indexes. If 5 input, cells contain 
% the indexes of the cover sets given by "Bal".
% "fig" is the figure number and "ms" is the marker size.


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

n = max(size(comps));
if n < 100
    col = repmat(col,[ceil(n/20),1]);
else
    col = rand(n,3);
end

S = comps{1};
if iscell(S)
    n = size(comps,1);
    for i = 1:n
        S = comps{i};
        if ~isempty(S)
            S = vertcat(S{:});
            comps{i} = S;
        else
            comps{i} = zeros(0,1);
        end
    end
end


if nargin == 4
    
    % Plot the segments
    figure(fig)
    C = comps{1};
    plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(1,:),'Markersize',ms)
    hold on
    for i = 2:n
        C = comps{i};
        plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
    end
    axis equal
    hold off
    pause(0.1)
    
else
    
    np = size(P,1);
    D = false(np,1);
    C = unique(vertcat(Bal{comps{1}}));
    figure(fig)
    plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(1,:),'Markersize',ms)
    hold on
    for i = 2:n
        if ~isempty(comps{i})
            C = unique(vertcat(Bal{comps{i}}));
            I = D(C);
            C = C(~I);
            D(C) = true;
            plot3(P(C,1),P(C,2),P(C,3),'.','Color',col(i,:),'Markersize',ms)
        end
    end
    hold off
    axis equal
    pause(0.1)
end