function R = rotation_matrix(A,angle)

% Returns the rotation matrix for the given axis A and angle (in radians)

A = A/norm(A);
R = zeros(3,3);
c = cos(angle); 
s = sin(angle);
R(1,:) = [A(1)^2+(1-A(1)^2)*c  A(1)*A(2)*(1-c)-A(3)*s  A(1)*A(3)*(1-c)+A(2)*s];
R(2,:) = [A(1)*A(2)*(1-c)+A(3)*s  A(2)^2+(1-A(2)^2)*c  A(2)*A(3)*(1-c)-A(1)*s];
R(3,:) = [A(1)*A(3)*(1-c)-A(2)*s  A(2)*A(3)*(1-c)+A(1)*s  A(3)^2+(1-A(3)^2)*c];
