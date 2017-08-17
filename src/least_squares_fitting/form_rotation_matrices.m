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

function [R,dR1,dR2] = form_rotation_matrices(theta)
  
% --------------------------------------------------------------------------
% FORM_ROTATION_MATRICES.M      Forms rotation matrices R = R2*R1 and its
%                                   derivatives
%
% Input 
% theta    Plane rotation angles (t1, t2) 
%
% Output 
% R        Rotation matrix
% R1       Plane rotation [1 0 0; 0 c1 -s1; 0 s1 c1]
% R2       Plane rotation [c2 0 s2; 0 1 0; -s2 0 c2]

c = cos(theta);
s = sin(theta);

R1 = [1 0 0; 0 c(1) -s(1); 0 s(1) c(1)];
R = R1;

R2 = [c(2) 0 s(2); 0 1 0; -s(2) 0 c(2)];
R = R2*R;

if nargout > 1
    dR1 = [0 0 0; 0 -R1(3,2) -R1(2,2); 0 R1(2,2) -R1(3,2)];
end

if nargout > 2
    dR2 = [-R2(1,3) 0 R2(1,1); 0 0 0; -R2(1,1) 0 -R2(1,3)];
end