function Set1 = set_difference(Set1,Set2,False)

% Performs the set difference so that the common elements of Set1 and Set2
% are removed from Set1, which is the output. Uses logical vector whose
% length must be up to the maximum element of the sets.

False(Set2) = true;
I = False(Set1);
Set1 = Set1(~I);