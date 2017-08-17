function [d,V,h,B] = distances_to_line(Q,LineDirec,LinePoint)

% Calculates the distances of the points, given in the rows of the
% matrix Q, to the line defined by one of its point and its direction.

if size(LineDirec,1) == 1
    LineDirec = LineDirec';
end
LineDirec = LineDirec/norm(LineDirec);

A = mat_vec_subtraction(Q,LinePoint);
h = A*LineDirec;

B = repmat(LineDirec',length(Q(:,1)),1);
B = [h.*B(:,1) h.*B(:,2) h.*B(:,3)];
V = A-B;

d = sqrt(sum(V.*V,2));