function A = mat_vec_subtraction(A,v)

% Subtracts from each row of the matrix A the vector v.
% If A is (n x m)-matrix, then v needs to be m-vector.

for i = 1:size(A,2)
    A(:,i) = A(:,i)-v(i);
end