function [A,L] = normalize(A)

% Normalize rows of the matrix

L = sqrt(sum(A.*A,2));
n = size(A,2);
for i = 1:n
    A(:,i) = A(:,i)./L;
end
