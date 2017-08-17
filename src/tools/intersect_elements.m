function Set = intersect_elements(Set1,Set2,False1,False2)

% Determines the intersection of Set1 and Set2.

Set = unique_elements([Set1; Set2],False1);
False1(Set1) = true;
False2(Set2) = true;
I = False1(Set)&False2(Set);
Set = Set(I);
