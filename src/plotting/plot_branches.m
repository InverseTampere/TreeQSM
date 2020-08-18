function plot_branches(P,cover,segment,fig,ms,segind,BO)

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
SPar = segment.ParentSegment;

if iscell(Segs{1})
    ns = max(size(Segs));
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

segments = segind;
C = SChi{segind};
b = 0;
while ~isempty(C) && b <= BO
    b = b+1;
    segments = [segments; C];
    C = SChi{segind};
end

ns = length(segment);
figure(fig)
for i = 1:ns
    if i == 2
        hold on
    end
    S = vertcat(Bal{Seg{segments(i)}});
    plot3(P(S,1),P(S,2),P(S,3),'.','Color',col(segments(i),:),'Markersize',ms)
end
hold off
