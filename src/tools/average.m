function A = average(X)

% Computes the average of columns of the matrix X

n = size(X,1);
if n > 1
    A = sum(X)/n;
else
    A = X;
end