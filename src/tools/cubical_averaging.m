function DSP = cubical_averaging(P,CubeSize)

tic
% Downsamples the given point cloud by averaging points from each 
% cube of side length CubeSize. 

% The vertices of the big cube containing P
Min = double(min(P));
Max = double(max(P));

% Number of cubes with edge length "EdgeLength" in the sides 
% of the big cube
N = double(ceil((Max-Min)/CubeSize)+1);

CubeCoord = floor([P(:,1)-Min(1) P(:,2)-Min(2) P(:,3)-Min(3)]/CubeSize)+1;

% Sorts the points according a lexicographical order
LexOrd = [CubeCoord(:,1) CubeCoord(:,2)-1 CubeCoord(:,3)-1]*[1 N(1) N(1)*N(2)]';
[LexOrd,SortOrd] = sort(LexOrd);
nc = size(unique(LexOrd),1); % number of points in the downsampled point cloud
np = size(P,1);     % number of points
DSP = zeros(nc,3);  % Downsampled point cloud
p = 1;              % The index of the point under comparison
q = 0;
while p <= np
    t = 1;
    while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
        t = t+1;
    end
    q = q+1;
    DSP(q,:) = average(P(SortOrd(p:p+t-1),:));
    p = p+t;
end
toc

disp(['    Points before:  ',num2str(np)])
disp(['  Filtered points:  ',num2str(np-nc)])
disp(['      Points left:  ',num2str(nc)]);
