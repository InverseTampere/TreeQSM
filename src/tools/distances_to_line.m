function [d,V,h,B] = distances_to_line(Q,LineDirec,LinePoint)

% Calculates the distances of the points, given in the rows of the
% matrix Q, to the line defined by one of its point and its direction.
% "LineDirec" must be a unit (1x3)-vector and LinePoint must be a (1x3)-vector.
% 
% Last update 8 Oct 2021

A = Q-LinePoint;
h = A*LineDirec';
B = h*LineDirec;
V = A-B;
d = sqrt(sum(V.*V,2));