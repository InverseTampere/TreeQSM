function [v,mean_res,sigmah,residual] = optimal_parallel_vector(V)

% For a given set of unit vectors (the rows of the matrix "V"),
% returns a unit vector ("v") that is the most parallel to them all
% in the sense that the sum of squared dot products of v with the
% vectors of V is maximized.
 
A = V'*V;
[U,~,~] = svd(A);
v = U(:,1)';

if nargout > 1
    residual = abs(V*v');
    mean_res = mean(residual);
    sigmah = std(residual);
end