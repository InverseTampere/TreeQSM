function [Vector,IndElements] = verticalcat(CellArray)

% Vertical concatenation of the given cell-array into a vector.

CellSize = cellfun('length',CellArray); % determine the size of each cell
nc = max(size(CellArray)); % number of cells
IndElements = ones(nc,2); % indexes for elements in each cell
IndElements(:,2) = cumsum(CellSize);
IndElements(2:end,1) = IndElements(2:end,1)+IndElements(1:end-1,2);
Vector = zeros(sum(CellSize),1); % concatenation of the cell-array into a vector
for j = 1:nc
    Vector(IndElements(j,1):IndElements(j,2)) = CellArray{j};
end