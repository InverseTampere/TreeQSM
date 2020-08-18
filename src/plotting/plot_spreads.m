function plot_spreads(treedata,fig,lw,rel)

% Plots the spreads as a polar plot with different height layers presented
% with different colors. Inputs "fig" and "lw" define the figure number and
% the line width. Input Rel = 1 specifies relative spreads, i.e. the
% maximum spread is one, otherwise use the actual values.

if nargin == 2
    lw = 1;
    rel = 1;
elseif nargin == 3
    rel = 1;
end

spreads = treedata.spreads;
figure(fig)
n = size(spreads,1);
col = zeros(n,3);
col(:,1) = (0:1/n:(n-1)/n)';
col(:,3) = (1:-1/n:1/n)';
d = max(max(spreads));
D = [spreads(1,end) spreads(1,:)];
if rel
    polarplot(D/d,'-','Color',col(1,:),'Linewidth',lw)
else
    polarplot(D,'-','Color',col(1,:),'Linewidth',lw)
end
hold on
for i = 1:n
    D = [spreads(i,end) spreads(i,:)];
    if rel
        polarplot(D/d,'-','Color',col(i,:),'Linewidth',lw)
    else
        polarplot(D,'-','Color',col(i,:),'Linewidth',lw)
    end
end
hold off
if rel
    rlim([0 1])
else
    rlim([0 d])
end
