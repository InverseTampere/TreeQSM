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

function [Intersect,IntersectLines] = check_self_intersection(Curve)

% The function takes in a curve (the coordinates of the vertices, in the
% right order) and checks if the curve intersects itself
%
% Outputs:
% Intersect         Logical value indicating if the curve self-intersects
% IntersectLines    Cell array containing for each line element which are
%                       the intersecting elements and how far away along
%                       the line the intersection point is


if ~isempty(Curve)
  dim = size(Curve,2); % two or three dimensional curve
  n = size(Curve,1); % number of points in the curve
  V = Curve([(2:n)'; 1],:)-Curve; % line elements forming the curve
  L = sqrt(sum(V.*V,2)); % the lengths of the line elements
  i = 1; % the line element under inspection
  Ind = (1:1:n)'; % indexes of the line elements
  if dim == 2 % 2d curves
    % directions (unit vectors) of the line elements:
    DirLines = [1./L.*V(:,1) 1./L.*V(:,2)]; 
    Intersect = false;
    if nargout == 1 % check only if the curve intersects
      while i <= n-1 && ~Intersect
        % Select the line elements that can intersect element i
        if i > 1
          I = Ind > i+1 | Ind < i-1;
        else
          I = Ind > i+1 & Ind < n;
        end
        ind = Ind(I)';
        for j = ind
          % Solve for the crossing points of every line element
          A = [DirLines(j,:)' -DirLines(i,:)'];
          b = Curve(i,:)'-Curve(j,:)';
          Ainv = 1/(A(1,1)*A(2,2)-A(1,2)*A(2,1))*[A(2,2) -A(1,2); -A(2,1) A(1,1)];
          x = Ainv*b; % signed length along the line elements to the crossing
          if x(1) >= 0 && x(1) <= L(j) && x(2) >= 0 && x(2) <= L(i)
            Intersect = true;
          end
        end
        i = i+1; % study the next line element
      end
    else % determine also all intersection points (line elements)
      IntersectLines = cell(n,2);
      for i = 1:n-1
        % Select the line elements that can intersect element i
        if i > 1
          I = Ind > i+1 | Ind < i-1;
        else
          I = Ind > i+1 & Ind < n;
        end
        ind = Ind(I)';
        for j = ind
          % Solve for the crossing points of every line element
          A = [DirLines(j,:)' -DirLines(i,:)'];
          b = Curve(i,:)'-Curve(j,:)';
          Ainv = 1/(A(1,1)*A(2,2)-A(1,2)*A(2,1))*[A(2,2) -A(1,2); -A(2,1) A(1,1)];
          x = Ainv*b;
          if x(1) >= 0 && x(1) <= L(j) && x(2) >= 0 && x(2) <= L(i)
            Intersect = true;
            % which line elements cross element i:
            IntersectLines{i,1} = [IntersectLines{i,1}; j]; 
            % which line elements cross element j:
            IntersectLines{j,1} = [IntersectLines{j,1}; i]; 
            % distances along element i to intersection points:
            IntersectLines{i,2} = [IntersectLines{i,2}; x(1)]; 
            % distances along element j to intersection points:
            IntersectLines{j,2} = [IntersectLines{j,2}; x(2)]; 
          end
        end
      end
      % remove possible multiple values
      for i = 1:n
        IntersectLines{i,1} = unique(IntersectLines{i,1});
        IntersectLines{i,2} = min(IntersectLines{i,2});
      end
    end

  elseif dim == 3 % 3d curves
    % directions (unit vectors) of the line elements
    DirLines = [1./L.*V(:,1) 1./L.*V(:,2) 1./L.*V(:,3)];
    Intersect = false;
    if nargout == 1 % check only if the curve intersects
      while i <= n-1
        % Select the line elements that can intersect element i
        if i > 1
          I = Ind > i+1 | Ind < i-1;
        else
          I = Ind > i+1 & Ind < n;
        end
        % Solve for possible intersection points
        [~,DistOnRay,DistOnLines] = distances_between_lines(...
          Curve(i,:),DirLines(i,:),Curve(I,:),DirLines(I,:));
        if any(DistOnRay >= 0 & DistOnRay <= L(i) &...
            DistOnLines > 0 & DistOnLines <= L(I))
          Intersect = true;
          i = n;
        else
          i = i+1; % study the next line element
        end
      end
    else % determine also all intersection points (line elements)
      IntersectLines = cell(n,2);
      for i = 1:n-1
        % Select the line elements that can intersect element i
        if i > 1
          I = Ind > i+1 | Ind < i-1;
        else
          I = Ind > i+1 & Ind < n;
        end
        % Solve for possible intersection points
        [D,DistOnRay,DistOnLines] = distances_between_lines(...
          Curve(i,:),DirLines(i,:),Curve(I,:),DirLines(I,:));
        if any(DistOnRay >= 0 & DistOnRay <= L(i) & ...
            DistOnLines > 0 & DistOnLines <= L(I))
          Intersect = true;
          J = DistOnRay >= 0 & DistOnRay <= L(i) & ...
            DistOnLines > 0 & DistOnLines <= L(I);
          ind = Ind(I);
          ind = ind(J);
          DistOnLines = DistOnLines(J);
          IntersectLines{i,1} = ind;
          IntersectLines{i,2} = DistOnRay(J);
          % Record the elements intersecting
          for j = 1:length(ind)
            IntersectLines{ind(j),1} = [IntersectLines{ind(j),1}; i];
            IntersectLines{ind(j),2} = [IntersectLines{ind(j),2}; DistOnLines(j)];
          end
        end
      end
      % remove possible multiple values
      for i = 1:n
        IntersectLines{i} = unique(IntersectLines{i});
        IntersectLines{i,2} = min(IntersectLines{i,2});
      end
    end
  end
else % Empty curve
  Intersect = false;
  IntersectLines = cell(1,1);
end
