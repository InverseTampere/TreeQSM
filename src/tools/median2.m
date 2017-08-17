function y = median2(X)

% Computes the median of the given vector.

n = size(X,1);
if n > 1
    X = sort(X);
    m = floor(n/2);
    if 2*m == n
        y = (X(m)+X(m+1))/2;
    elseif m == 0
        y = (X(1)+X(2))/2;
    else
        y = X(m+1);
    end
else
    y = X;
end