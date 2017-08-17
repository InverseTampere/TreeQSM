function C = dot_product(A,B)

% Computes the dot product of the corresponding rows of the matrices A and B

C = sum(A.*B,2);