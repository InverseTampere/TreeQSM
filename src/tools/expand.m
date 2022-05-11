function C = expand(Nei,C,n,Forb)

% Expands the given subset "C" of cover sets "n" times with their neighbors, 
% and optionally, prevents the expansion into "Forb" sets. "C" is a vector 
% and "Forb" can be a number vector or a logical vector.

if nargin == 3
  for i = 1:n
    C = union(C,vertcat(Nei{C}));
  end
  if size(C,2) > 1
    C = C';
  end
else
  if islogical(Forb)
    for i = 1:n
      C = union(C,vertcat(Nei{C}));
      I = Forb(C);
      C = C(~I);
    end
  else
    for i = 1:n
      C = union(C,vertcat(Nei{C}));
      C = setdiff(C,Forb);
    end
  end
  if size(C,2) > 1
    C = C';
  end
end