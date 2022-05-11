function Pass = cubical_downsampling(P,CubeSize)

% Downsamples the given point cloud by selecting one point from each 
% cube of side length "CubeSize". 

% The vertices of the big cube containing P
Min = double(min(P));
Max = double(max(P));

% Number of cubes with edge length "EdgeLength" in the sides 
% of the big cube
N = double(ceil((Max-Min)/CubeSize)+1);

% Process the data in 1e7-point blocks to consume much less memory 
np = size(P,1);
m = 1e7;
if np < m
    m = np;
end
nblocks = ceil(np/m); % number of blocks

% Downsample
R = cell(nblocks,1);
p = 1;
for i = 1:nblocks
    if i < nblocks
        % Compute the cube coordinates of the points
        C = floor([double(P(p:p+m-1,1))-Min(1) double(P(p:p+m-1,2))-Min(2)...
            double(P(p:p+m-1,3))-Min(3)]/CubeSize)+1;
        % Compute the lexicographical order of the cubes
        S = [C(:,1) C(:,2)-1 C(:,3)-1]*[1 N(1) N(1)*N(2)]';
        [S,I] = unique(S); % Select the unique cubes
        J = (p:1:p+m-1)';
        J = J(I);
        R{i} = [S J];
    else
        C = floor([double(P(p:end,1))-Min(1) double(P(p:end,2))-Min(2)...
            double(P(p:end,3))-Min(3)]/CubeSize)+1;
        S = [C(:,1) C(:,2)-1 C(:,3)-1]*[1 N(1) N(1)*N(2)]';
        [S,I] = unique(S);
        J = (p:1:np)';
        J = J(I);
        R{i} = [S J];
    end
    p = p+m;
end
% Select the unique cubes and their points
R = vertcat(R{:});
[~,I] = unique(R(:,1));
Pass = false(np,1);
Pass(R(I,2)) = true;
