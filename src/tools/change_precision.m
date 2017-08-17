function v = change_precision(v)

% Decrease the number of nonzero decimals in the vector v according to the
% exponent of the number for displaying and writing.

n = length(v);
for i = 1:n
    if abs(v(i)) >= 1e3
        v(i) = round(v(i));
    elseif abs(v(i)) >= 1e2
        v(i) = round(10*v(i))/10;
    elseif abs(v(i)) >= 1e1
        v(i) = round(100*v(i))/100;
    elseif abs(v(i)) >= 1e0
        v(i) = round(1000*v(i))/1000;
    elseif abs(v(i)) >= 1e-1
        v(i) = round(10000*v(i))/10000;
    else
        v(i) = round(100000*v(i))/100000;
    end
end