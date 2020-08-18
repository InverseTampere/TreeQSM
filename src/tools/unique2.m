function SetUni = unique2(Set)


n = length(Set);
if n > 0
    Set = sort(Set);
    d = Set(2:n)-Set(1:n-1);
    A = Set(2:n);
    I = d > 0;
    SetUni = [Set(1); A(I)];
else
    SetUni = Set;
end