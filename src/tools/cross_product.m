function C = cross_product(A,B)

% Calculates the cross product C of the 3-vectors A and B

C = [A(2)*B(3)-A(3)*B(2);  A(3)*B(1)-A(1)*B(3);  A(1)*B(2)-A(2)*B(1)]; 