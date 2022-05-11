function plot_large_point_cloud(P,fig,ms,rel)

% Plots a random subset of a large point cloud. The user specifies the
% relative size of the subset (input "rel" given as in percentage points).
%
% Inputs:
% P     Point cloud
% fig   Figure number
% ms    Marker size
% rel   Subset size in percentage points (%). 
%           E.g. if rel = 12, then about 12 % poinst are plotted

rel = 0.5/(1-rel/100); % Compute a coeffiecient

I = logical(round(rel*rand(size(P,1),1)));
plot_point_cloud(P(I,:),fig,ms)