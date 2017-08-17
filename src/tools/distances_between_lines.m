function [DistLines,DistOnRay,DistOnLines] = distances_between_lines(PointRay,DirRay,PointLines,DirLines)

% Calculates the distances between a ray and lines

% PointRay      A point of the ray
% DirRay        Unit direction vector of the line
% PointLines    One point of every line
% DirLines      Unit direction vectors of the lines 

PointLines = double(PointLines);
PointRay = double(PointRay);
DirLines = double(DirLines);
DirRay = double(DirRay);

% Calculate unit vectors N orthogonal to the ray and the lines
N = [DirRay(2)*DirLines(:,3)-DirRay(3)*DirLines(:,2) ...
    DirRay(3)*DirLines(:,1)-DirRay(1)*DirLines(:,3) ...
    DirRay(1)*DirLines(:,2)-DirRay(2)*DirLines(:,1)];
l = sqrt(sum(N.*N,2));
N = [1./l.*N(:,1) 1./l.*N(:,2) 1./l.*N(:,3)];

% Calculate the distances between the lines
A = -mat_vec_subtraction(PointLines,PointRay);
DistLines = sqrt(abs(sum(A.*N,2))); % distance between lines and the ray

% Calculate the distances on ray and on lines
b = DirLines*DirRay';
d = A*DirRay';
e = sum(A.*DirLines,2);
DistOnRay = (b.*e-d)./(1-b.^2); % Distances to PointRay from the closest points on the ray
DistOnLines = (e-b.*d)./(1-b.^2); % Distances to PointLines from the closest points on the lines
