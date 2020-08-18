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

function save_model_text(QSM,savename)

% ---------------------------------------------------------------------
% SAVE_MODEL_TEXT.M       Saves QSM (cylinder, branch, treedata) into text
%                           files
%
% Version 1.1.0
% Latest update     17 Aug 2020
%
% Copyright (C) 2013-2020 Pasi Raumonen
% ---------------------------------------------------------------------

% Save the cylinder, branch, and treedata structures in text-formats (.txt) 
% into /result-folder with the input "savename" defining the file names:
% 'cylinder_',savename,'.txt'
% 'branch_',savename,'.txt'
% 'treedata_',savename,'.txt'
% !!! Notice that only part of the treedata, the single number tree 
% attributes are saved in the text-file.
% Every user can change this code easily to define what is saved into 
% their text-files.

% Changes from version 1.0.0 to 1.1.0, 17 Aug 2020:
% 1) Added the new fields of cylinder, branch and treedata structures
% 2) Added header names to the files
% 3) Changed the names of the files to be saved
% 4) Changed the name of second input from "string" to "savename"
% 5) Changed the rounding of some parameters and attributes

cylinder = QSM.cylinder;
branch = QSM.branch;
treedata = QSM.treedata;

%% Form cylinder data, branch data and tree data
% Use less decimals
Rad = round(10000*cylinder.radius)/10000; % radius (m)
Len = round(10000*cylinder.length)/10000; % length (m)
Sta = round(10000*cylinder.start)/10000; % starting point (m)
Axe = round(10000*cylinder.axis)/10000; % axis (m)
CPar = single(cylinder.parent); % parent cylinder
CExt = single(cylinder.extension); % extension cylinder
Added = single(cylinder.added); % is cylinder added to fil a gap
Rad0 = round(10000*cylinder.UnmodRadius)/10000; % unmodified radius (m)
B = single(cylinder.branch); % branch index of the cylinder
BO = single(cylinder.BranchOrder); % branch order of the branch
PIB = single(cylinder.PositionInBranch); % position of the cyl. in the branch
Mad = single(round(10000*cylinder.mad)/10000); % mean abso. distance (m)
SC = single(round(10000*cylinder.SurfCov)/10000); % surface coverage
CylData = [Rad Len Sta Axe CPar CExt B BO PIB Mad SC Added Rad0];
NamesC = ['radius (m)',"length (m)","start_point","axis_direction",...
  "parent","extension","branch","branch_order","position_in_branch",...
  "mad","SurfCov","added","UnmodRadius (m)"];

BOrd = single(branch.order); % branch order
BPar = single(branch.parent); % parent branch
BDia = round(10000*branch.diameter)/10000; % diameter (m)
BVol = round(10000*branch.volume)/10000; % volume (L)
BAre = round(10000*branch.area)/10000; % area (m^2)
BLen = round(1000*branch.length)/1000; % length (m)
BAng = round(10*branch.angle)/10; % angle (deg)
BHei = round(1000*branch.height)/1000; % height (m)
BAzi = round(10*branch.azimuth)/10; % azimuth (deg)
BZen = round(10*branch.zenith)/10; % zenith (deg)
BranchData = [BOrd BPar BDia BVol BAre BLen BHei BAng BAzi BZen];
NamesB = ["order","parent","diameter (m)","volume (L)","area (m^2)",...
  "length (m)","height (m)","angle (deg)","azimuth (deg)","zenith (deg)"];

% Extract the field names of treedata
Names = fieldnames(treedata);
n = 1;
while ~strcmp(Names{n},'location')
    n = n+1;
end
n = n-1;
Names = Names(1:n);

TreeData = zeros(n,1); 
% TreeData contains TotalVolume, TrunkVolume, BranchVolume, etc
for i = 1:n
    TreeData(i) = treedata.(Names{i,:});
end
TreeData = change_precision(TreeData); % use less decimals
NamesD = string(Names);

%% Save the data as text-files
str = ['results/cylinder_',savename,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%s\t', 1, size(NamesC,2)-1) '%s\n'], NamesC.');
fprintf(fid, [repmat('%g\t', 1, size(CylData,2)-1) '%g\n'], CylData.');
fclose(fid);

str = ['results/branch_',savename,'.txt'];
fid = fopen(str, 'wt');
fprintf(fid, [repmat('%s\t', 1, size(NamesB,2)-1) '%s\n'], NamesB.');
fprintf(fid, [repmat('%g\t', 1, size(BranchData,2)-1) '%g\n'], BranchData.');
fclose(fid);

str = ['results/treedata_',savename,'.txt'];
fid = fopen(str, 'wt');
NamesD(:,2) = TreeData;
fprintf(fid,'%s\t %g\n',NamesD.');
fclose(fid);
