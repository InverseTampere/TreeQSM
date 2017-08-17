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

function save_model_text(QSM,string)

% Save the cylinder, branch, and treedata structures in text-formats (.txt) 
% into /result-folder

% !!! Notice that only part of the treedata is saved in the text-file and
% the treedata.txt is different than in previous versions.
% Every user can change this code easily to define what is saved into 
% their text-files.

cylinder = QSM.cylinder;
branch = QSM.branch;
treedata = QSM.treedata;

%% Form cylinder data, branch data and tree data
% Use less decimals
Rad = round(10000*cylinder.radius)/10000; % radius in meters
Len = round(10000*cylinder.length)/10000; % length in meters
Sta = round(10000*cylinder.start)/10000; % starting point in meters
Axe = round(10000*cylinder.axis)/10000; % axis in meters
CPar = single(cylinder.parent);
CExt = single(cylinder.extension);
Added = single(cylinder.added);
Rad0 = round(10000*cylinder.UnmodRadius)/10000;
B = single(cylinder.branch);
BO = single(cylinder.BranchOrder);
PIB = single(cylinder.PositionInBranch);
CylData = [Rad Len Sta Axe CPar CExt B BO PIB Added Rad0];

BOrd = single(branch.order);
BPar = single(branch.parent);
BVol = round(10000*branch.volume)/10000; % volume in litres
BLen = round(10000*branch.length)/10000; % length in meters
BAng = branch.angle; % angle in degrees
BHei = branch.height; % height in metres
BAzi = branch.azimuth; % azimuth in degrees
BDia = round(10000*branch.diameter)/10000; % diameter in meters
BranchData = [BOrd BPar BVol BLen BAng BHei BAzi BDia];

Names = fieldnames(treedata);
n = 1;
while ~strcmp(Names{n}(1:3),'DBH')
    n = n+1;
end
if size(Names,1) <= 18
    % No triangulation model
    n = n+1;
else
    % Triangulation model
    n = n+6;
end
TreeData = zeros(n,1); 
% TreeData contains TotalVolume, TrunkVolume, BranchVolume, TreeHeight, 
% TrunkLength, BranchLength, NumberBranches, MaxBranchOrder, TotalArea, DBHqsm, DBHcyl.
% If triangulation model, then TreeData contains the following additional
% data: DBHtri, TriaTrunkVolume, MixTrunkVolume, MixTotalVolume, TriaTrunkLength 
for i = 1:n
    TreeData(i) = treedata.(Names{i,:});
end
% % Add branch order data up to sixth order as was used in some earlier versions
% % in the treedata.txt
% for i = 1:min(6,treedata.MaxBranchOrder)
%     TreeData(n+i) = treedata.NumberBranchOrder(i);
% end
% n = n+6;
% for i = 1:min(6,treedata.MaxBranchOrder)
%     TreeData(n+i) = treedata.VolumeBranchOrder(i);
% end
% n = n+6;
% for i = 1:min(6,treedata.MaxBranchOrder)
%     TreeData(n+i) = treedata.LengthBranchOrder(i);
% end
% if treedata.MaxBranchOrder < 6
%     n = n+6;
%     TreeData(n) = 0;
% end

TreeData = change_precision(TreeData); % use less decimals

%% Save the data as text-files
str = ['results/cyl_data_',string,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%g\t', 1, size(CylData,2)-1) '%g\n'], CylData.');
fclose(fid);

str = ['results/branch_data_',string,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%g\t', 1, size(BranchData,2)-1) '%g\n'], BranchData.');
fclose(fid);

str = ['results/tree_data_',string,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%g\t', 1, size(TreeData,2)-1) '%g\n'], TreeData.');
fclose(fid);
