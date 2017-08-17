% This file is part of TREEQSM.
% 
% TREEQSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TREEQSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.

function [Components,CompSize] = connected_components(Nei,Sub,MinSize,Fal)

% ---------------------------------------------------------------------
% CONNECTED_COMPONENTS.M      Determines the connected components of cover
%                                   sets using their neighbour-relation
%
% Version 1.1
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Determines connected components of the subset of cover sets defined
% by "Sub" such that each component has at least "MinSize"
% number of cover sets.
%
% Inputs:
% Nei       Neighboring cover sets of each cover set, (n_sets x 1)-cell
% Sub       Subset whose components are determined,
%               length(Sub) < 2 means no subset and thus the whole point cloud
%               "Sub" may be also a vector of cover set indexes in the subset
%               or a logical (n_sets)-vector, where n_sets is the number of
%               all cover sets
% MinSize   Minimum number of cover sets in an acceptable component
% Fal       Logical false vector for the cover sets
%
% Outputs:
% Components    Connected components, (n_comp x 1)-cell
% CompSize      Number of sets in the components, (n_comp x 1)-vector

if length(Sub) <= 3 && ~islogical(Sub) && Sub(1) > 0
    % Very small subset, i.e. at most 3 cover sets
    n = length(Sub);
    if n == 1
        Components = cell(1,1);
        Components{1} = uint32(Sub);
        CompSize = 1;
    elseif n == 2
        I = Nei{Sub(1)} == Sub(2);
        if any(I)
            Components = cell(1,1);
            Components{1} = uint32((Sub));
            CompSize = 1;
        else
            Components = cell(2,1);
            Components{1} = uint32(Sub(1));
            Components{2} = uint32(Sub(2));
            CompSize = [1 1];
        end
    elseif n == 3
        I = Nei{Sub(1)} == Sub(2);
        J = Nei{Sub(1)} == Sub(3);
        K = Nei{Sub(2)} == Sub(3);
        if any(I)+any(J)+any(K) >= 2
            Components = cell(1,1);
            Components{1} = uint32(Sub);
            CompSize = 1;
        elseif any(I)
            Components = cell(2,1);
            Components{1} = uint32(Sub(1:2));
            Components{2} = uint32(Sub(3));
            CompSize = [2 1];
        elseif any(J)
            Components = cell(2,1);
            Components{1} = uint32(Sub([1 3]));
            Components{2} = uint32(Sub(2));
            CompSize = [2 1];
        elseif any(K)
            Components = cell(2,1);
            Components{1} = uint32(Sub(2:3));
            Components{2} = uint32(Sub(1));
            CompSize = [2 1];
        else
            Components = cell(3,1);
            Components{1} = uint32(Sub(1));
            Components{2} = uint32(Sub(2));
            Components{3} = uint32(Sub(3));
            CompSize = [1 1 1];
        end
    end
    
elseif any(Sub) || (length(Sub) == 1 && Sub(1) == 0)
    nb = size(Nei,1);
    if nargin == 3
        Fal = false(nb,1);
    end
    if length(Sub) == 1 && Sub == 0
        % All the cover sets
        ns = nb;
        if nargin == 3
            Sub = true(nb,1);
        else
            Sub = ~Fal;
        end
    elseif ~islogical(Sub)
        % Subset of cover sets
        ns = length(Sub);
        if nargin == 3
            sub = false(nb,1);
        else
            sub = Fal;
        end
        sub(Sub) = true;
        Sub = sub;
    else
        % Subset of cover sets
        ns = nnz(Sub);
    end
    
    Components = cell(ns,1);
    CompSize = zeros(ns,1,'uint32');
    nc = 0;      % number of components found
    m = 1;
    while ~Sub(m)
        m = m+1;
    end
    i = 0;
    Comp = zeros(ns,1,'uint32');
    while i < ns
        Add = Nei{m};
        I = Sub(Add);
        Add = Add(I);
        a = length(Add);
        Comp(1) = m;
        Sub(m) = false;
        t = 1;
        while a > 0
            Comp(t+1:t+a) = Add;
            Sub(Add) = false;
            t = t+a;
            Add = vertcat(Nei{Add});
            I = Sub(Add);
            Add = Add(I);
            % select the unique elements of Add:
            n = length(Add);
            if n > 2
                I = true(n,1);
                for j = 1:n
                    if ~Fal(Add(j))
                        Fal(Add(j)) = true;
                    else
                        I(j) = false;
                    end
                end
                Fal(Add) = false;
                Add = Add(I);
            elseif n == 2
                if Add(1) == Add(2)
                    Add = Add(1);
                end
            end
            a = length(Add);
        end
        i = i+t;
        if t >= MinSize
            nc = nc+1;
            Components{nc} = uint32(Comp(1:t));
            CompSize(nc) = t;
        end
        if i < ns
            while m <= nb && Sub(m) == false
                m = m+1;
            end
        end
    end
    Components = Components(1:nc);
    CompSize = CompSize(1:nc);
else
    Components = cell(0,1);
    CompSize = 0;
end