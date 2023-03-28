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


function QSM = treeqsm(P,inputs)

% ---------------------------------------------------------------------
% TREEQSM.M     Reconstructs quantitative structure tree models from point 
%                   clouds containing a tree.
%
% Version 2.4.1
% Latest update     2 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------
%
% INPUTS:
%
% P                 (Filtered) point cloud, (m_points x 3)-matrix, the rows
%                       give the coordinates of the points.
%
% inputs            Structure field defining reconstruction parameters.
%                       Created with the "create_input.m" script. Contains 
%                       the following main fields:
%   PatchDiam1        Patch size of the first uniform-size cover
%
%   PatchDiam2Min     Minimum patch size of the cover sets in the second cover
%
%   PatchDiam2Max     Maximum cover set size in the stem's base in the 
%                       second cover
%
%   BallRad1          Ball size used for the first cover generation
%
%   BallRad2          Maximum ball radius used for the second cover generation
%
%   nmin1             Minimum number of points in BallRad1-balls, 
%                       default value is 3.
%
%   nmin2             Minimum number of points in BallRad2-balls, 
%                       default value is 1.
%
%   OnlyTree          If "1", the point cloud contains only points from the 
%                       tree and the trunk's base is defined as the lowest 
%                       part of the point cloud. Default value is "1". 
%
%   Tria              If "1", tries to make triangulation for the stem up 
%                       to first main branch. Default value is "0". 
%
%   Dist              If "1", compute the point-model distances. 
%                       Default value is "1".
%
%   MinCylRad         Minimum cylinder radius, used particularly in the 
%                       taper corrections
%
%   ParentCor         If "1", child branch cylinders radii are always 
%                       smaller than the parent branche's cylinder radii
%
%   TaperCor          If "1", use partially linear (stem) and parabola 
%                       (branches) taper corrections
%
%   GrowthVolCor      If "1", use growth volume correction introduced 
%                       by Jan Hackenberg
%
%   GrowthVolFac      fac-parameter of the growth volume approach, 
%                       defines upper and lower bound
%
%   name              Name string for saving output files and name for the
%                       model in the output object
% 
%   tree              Numerical id/index given to the tree
% 
%   model             Model number of the tree, e.g. with the same inputs
%
%   savemat           If "1", saves the output struct QSM as a matlab-file
%                       into \result folder 
%
%   savetxt           If "1", saves the models in .txt-files into 
%                       \result folder 
%
%   plot              Defines what is plotted during the reconstruction:
%                       2 = same as below plus distributions
%                       1 = plots the segmented point cloud and QSMs
%                       0 = plots nothing
%
%   disp              Defines what is displayed during the reconstruction:
%                       2 = same as below plus times and tree attributes; 
%                       1 = display name, parameters and fit metrics;
%                       0 = display only the name
% ---------------------------------------------------------------------
% OUTPUT:
%
% QSM           Structure array with the following fields:
%               cylinder        Cylinder data  
%               branch          Branch data
%               treedata        Tree attributes  
%               rundata         Information about the modelling run
%               pmdistances     Point-to-model distance statistics
%               triangulation   Triangulation of the stem (if inputs.Tria = 1)
% ---------------------------------------------------------------------

% cylinder (structure-array) contains the following fields:
% radius
% length
% start         xyz-coordinates of the starting point
% axis          xyz-component of the cylinder axis
% parent        index (in this file) of the parent cylinder
% extension     index (in this file) of the extension cylinder
% added         is cylinder added after normal cylinder fitting (= 1 if added)
% UnmodRadius   unmodified radius of the cylinder
% branch        branch (index in the branch structure array) of the cylinder
% BranchOrder   branch order of the branch the cylinder belongs
% PositionInBranch	running number of the cylinder in the branch it belongs
%
% branch (structure-array) contains the following fields:
% order     branch order (0 for trunk, 1 for branches originating from 
%               the trunk, etc.)
% parent	index (in this file) of the parent branch
% volume	volume (L) of the branch (sum of the volumes of the cylinders 
%               forming the branch)
% length	length (m) of the branch (sum of the lengths of the cylinders)
% angle     branching angle (deg) (angle between the branch and its parent 
%               at the branching point)
% height    height (m) of the base of the branch
% azimuth   azimuth (deg) of the branch at the base 
% diameter  diameter (m) of the branch at the base
%
% treedata (structure-array) contains the following fields:
% TotalVolume
% TrunkVolume
% BranchVolume
% TreeHeight
% TrunkLength
% BranchLength
% NumberBranches    Total number of branches
% MaxBranchOrder 
% TotalArea 
% DBHqsm        From the cylinder of the QSM at the right heigth
% DBHcyl        From the cylinder fitted to the section 1.1-1.5m
% location      (x,y,z)-coordinates of the base of the tree
% StemTaper     Stem taper function/curve from the QSM
% VolumeCylDiam     Distribution of the total volume in diameter classes
% LengthCylDiam     Distribution of the total length in diameter classes
% VolumeBranchOrder     Branch volume per branching order
% LengthBranchOrder     Branch length per branching order
% NumberBranchOrder     Number of branches per branching order

% treedata from mixed model (cylinders and triangulation) contains also 
% the following fields:
% DBHtri            Computed from triangulation model
% TriaTrunkVolume   Triangulated trunk volume (up to first branch)
% MixTrunkVolume    Mixed trunk volume, bottom (triang.) + top (cylinders)
% MixTotalVolume    Mixed total volume, mixed trunk volume + branch volume
% TriaTrunkLength   Triangulated trunk length
%
% pmdistances (structure-array) contains the following fields (and others):
% CylDists  Average point-model distance for each cylinder
% median    median of CylDist for all, stem, 1branch, 2branch cylinder
% mean      mean of CylDist for all, stem, 1branch, 2branch cylinder
% max       max of CylDist for all, stem, 1branch, 2branch cylinder
% std       standard dev. of CylDist for all, stem, 1branch, 2branch cylinder
% 
% rundata (structure-array) contains the following fields:
% inputs    The input parameters in a structure-array
% time      Computation times for each step
% date      Starting and stopping dates (year,month,day,hour,minute,second) 
%             of the computation
% 
% triangulation (structure-array) contains the following fields:
% vert      Vertices (xyz-coordinates) of the triangulation
% facet     Facet information
% fvd       Color information for plotting the model
% volume    Volume enclosed by the triangulation
% bottom    Z-coordinate of the bottom plane of the triangulation
% top       Z-coordinate of the top plane of the triangulation
% triah     Height of the triangles
% triah     Width of the triangles
% cylind    Cylinder index in the stem where the triangulation stops
% ---------------------------------------------------------------------

% Changes from version 2.4.0 to 2.4.1, 2 May 2022:  
% Minor update. New filtering options, new code ("define_input") for 
% selecting automatically PatchDiam and BallRad parameter values for 
% the optimization process, added sensitivity estimates of the results, 
% new smoother plotting of QSMs, corrected some bugs, rewrote some 
% functions (e.g. "branches").
% Particular changes in treeqsm.m file:
% 1) Deleted the remove of the field "ChildCyls" and "CylsInSegment".

% Changes from version 2.3.2 to 2.4.0, 17 Aug 2020:  
% First major update. Cylinder fitting process and the taper correction 
% has changed. The fitting is adaptive and no more “lcyl” and “FilRad” 
% parameters. Treedata has many new outputs: Branch and cylinder 
% distributions; surface areas; crown dimensions. More robust triangulation 
% of stem. Branch, cylinder and triangulation structures have new fields. 
% More optimisation metrics, more plots of the results and more plotting 
% functions.
% Particular changes in treeqsm.m file:
% 1) Removed the for-loops for lcyl and FilRad.
% 2) Changes what is displayed about the quality of QSMs 
%    (point-model-distances and surface coverage) during reconstruction
% 3) Added version number to rundata
% 4) Added remove of the field "ChildCyls" and "CylsInSegment" of "cylinder"
%    from "branches" to "treeqsm".

% Changes from version 2.3.1 to 2.3.2, 2 Dec 2019:  
% Small changes in the subfunction to allow trees without branches

% Changes from version 2.3.0 to 2.3.1, 8 Oct 2019:  
% 1) Some changes in the subfunctions, particularly in "cylinders" and 
%    "tree_sets"
% 2) Changed how "treeqsm" displays things during the running of the
%    function


%% Code starts -->
Time = zeros(11,1); % Save computation times for modelling steps
Date = zeros(2,6); % Starting and stopping dates of the computation
Date(1,:) = clock;
% Names of the steps to display
name = ['Cover sets      ';
        'Tree sets       ';
        'Initial segments';
        'Final segments  ';
        'Cylinders       ';
        'Branch & data   ';
        'Distances       '];
 
if inputs.disp > 0
  disp('---------------')
  disp(['  ',inputs.name,', Tree = ',num2str(inputs.tree),...
    ', Model = ',num2str(inputs.model)])
end

% Input parameters
PatchDiam1 = inputs.PatchDiam1;
PatchDiam2Min = inputs.PatchDiam2Min;
PatchDiam2Max = inputs.PatchDiam2Max;
BallRad1 = inputs.BallRad1; 
BallRad2 = inputs.BallRad2; 
nd = length(PatchDiam1);
ni = length(PatchDiam2Min);
na = length(PatchDiam2Max);

if inputs.disp == 2
  % Display parameter values
  disp(['  PatchDiam1 = ',num2str(PatchDiam1)])
  disp(['  BallRad1 = ',num2str(BallRad1)])
  disp(['  PatchDiam2Min = ',num2str(PatchDiam2Min)])
  disp(['  PatchDiam2Max = ',num2str(PatchDiam2Max)])
  disp(['  BallRad2 = ',num2str(BallRad2)])
  disp(['  Tria = ',num2str(inputs.Tria),...
      ', OnlyTree = ',num2str(inputs.OnlyTree)])
  disp('Progress:')
end

P=transpose(P)

%% Make the point cloud into proper form
% only 3-dimensional data
if size(P,2) > 3
    P = P(:,1:3);
end
% Only double precision data
if ~isa(P,'double')
    P = double(P);
end

%% Initialize the output file
QSM = struct('cylinder',{},'branch',{},'treedata',{},'rundata',{},...
    'pmdistance',{},'triangulation',{});

%% Reconstruct QSMs
nmodel = 0;
for h = 1:nd
  tic
  Inputs = inputs;
  Inputs.PatchDiam1 = PatchDiam1(h);
  Inputs.BallRad1 = BallRad1(h);
  if nd > 1 && inputs.disp >= 1
    disp('  -----------------')
    disp(['  PatchDiam1 = ',num2str(PatchDiam1(h))]);
    disp('  -----------------')
  end
  
  %% Generate cover sets
  cover1 = cover_sets(P,Inputs);
  Time(1) = toc;
  if inputs.disp == 2
    display_time(Time(1),Time(1),name(1,:),1)
  end
  
  %% Determine tree sets and update neighbors
  [cover1,Base,Forb] = tree_sets(P,cover1,Inputs);
  Time(2) = toc-Time(1);
  if inputs.disp == 2
    display_time(Time(2),sum(Time(1:2)),name(2,:),1)
  end
  
  %% Determine initial segments
  segment1 = segments(cover1,Base,Forb);
  Time(3) = toc-sum(Time(1:2));
  if inputs.disp == 2
    display_time(Time(3),sum(Time(1:3)),name(3,:),1)
  end
  
  %% Correct segments
  % Don't remove small segments and add the modified base to the segment
  segment1 = correct_segments(P,cover1,segment1,Inputs,0,1,1);
  Time(4) = toc-sum(Time(1:3));
  if inputs.disp == 2
    display_time(Time(4),sum(Time(1:4)),name(4,:),1)
  end
  
  for i = 1:na
    % Modify inputs
    Inputs.PatchDiam2Max = PatchDiam2Max(i);
    Inputs.BallRad2 = BallRad2(i);
    if na > 1 && inputs.disp >= 1
      disp('    -----------------')
      disp(['    PatchDiam2Max = ',num2str(PatchDiam2Max(i))]);
      disp('    -----------------')
    end
    for j = 1:ni
      tic
      % Modify inputs
      Inputs.PatchDiam2Min = PatchDiam2Min(j);
      if ni > 1 && inputs.disp >= 1
        disp('      -----------------')
        disp(['      PatchDiam2Min = ',num2str(PatchDiam2Min(j))]);
        disp('      -----------------')
      end
      
      %% Generate new cover sets
      % Determine relative size of new cover sets and use only tree points
      RS = relative_size(P,cover1,segment1);
      
      % Generate new cover
      cover2 = cover_sets(P,Inputs,RS);
      Time(5) = toc;
      if inputs.disp == 2
          display_time(Time(5),sum(Time(1:5)),name(1,:),1)
      end
      
      %% Determine tree sets and update neighbors
      [cover2,Base,Forb] = tree_sets(P,cover2,Inputs,segment1);
      Time(6) = toc-Time(5);
      if inputs.disp == 2
        display_time(Time(6),sum(Time(1:6)),name(2,:),1)
      end
      
      %% Determine segments
      segment2 = segments(cover2,Base,Forb);
      Time(7) = toc-sum(Time(5:6));
      if inputs.disp == 2
        display_time(Time(7),sum(Time(1:7)),name(3,:),1)
      end
      
      %% Correct segments
      % Remove small segments and the extended bases.
      segment2 = correct_segments(P,cover2,segment2,Inputs,1,1,0);
      Time(8) = toc-sum(Time(5:7));
      if inputs.disp == 2
        display_time(Time(8),sum(Time(1:8)),name(4,:),1)
      end
      
      %% Define cylinders
      cylinder = cylinders(P,cover2,segment2,Inputs);
      Time(9) = toc;
      if inputs.disp == 2
        display_time(Time(9),sum(Time(1:9)),name(5,:),1)
      end
      
      if ~isempty(cylinder.radius)
        %% Determine the branches
        branch = branches(cylinder);
        
        %% Compute (and display) model attributes
        T = segment2.segments{1};
        T = vertcat(T{:});
        T = vertcat(cover2.ball{T});
        trunk = P(T,:); % point cloud of the trunk
        % Compute attributes and distibutions from the cylinder model
        % and possibly some from a triangulation
        [treedata,triangulation] = tree_data(cylinder,branch,trunk,inputs);
        Time(10) = toc-Time(9);
        if inputs.disp == 2
          display_time(Time(10),sum(Time(1:10)),name(6,:),1)
        end
        
        %% Compute point model distances
        if inputs.Dist
          pmdis = point_model_distance(P,cylinder);
          
          % Display the mean point-model distances and surface coverages
          % for stem, branch, 1branc and 2branch cylinders
          if inputs.disp >= 1
            D = [pmdis.TrunkMean pmdis.BranchMean ...
                pmdis.Branch1Mean pmdis.Branch2Mean];
            D = round(10000*D)/10;
            
            T = cylinder.branch == 1;
            B1 = cylinder.BranchOrder == 1;
            B2 = cylinder.BranchOrder == 2;
            SC = 100*cylinder.SurfCov;
            S = [mean(SC(T)) mean(SC(~T)) mean(SC(B1)) mean(SC(B2))];
            S = round(10*S)/10;
            
            disp('  ----------')
            str = ['  PatchDiam1 = ',num2str(PatchDiam1(h)), ...
                ', PatchDiam2Max = ',num2str(PatchDiam2Max(i)), ...
                ', PatchDiam2Min = ',num2str(PatchDiam2Min(j))];
            disp(str)
            str = ['  Distances and surface coverages for ',...
                'trunk, branch, 1branch, 2branch:'];
            disp(str)
            str = ['  Average cylinder-point distance:  '...
                num2str(D(1)),'  ',num2str(D(2)),'  ',...
                num2str(D(3)),'  ',num2str(D(4)),' mm'];
            disp(str)
            str = ['  Average surface coverage:  '...
                num2str(S(1)),'  ',num2str(S(2)),'  ',...
                num2str(S(3)),'  ',num2str(S(4)),' %'];
            disp(str)
            disp('  ----------')
          end
          Time(11) = toc-sum(Time(9:10));
          if inputs.disp == 2
            display_time(Time(11),sum(Time(1:11)),name(7,:),1)
          end
        end
        
        %% Reconstruct the output "QSM"
        Date(2,:) = clock;
        Time(12) = sum(Time(1:11));
        clear qsm
        qsm = struct('cylinder',{},'branch',{},'treedata',{},'rundata',{},...
          'pmdistance',{},'triangulation',{});
        qsm(1).cylinder = cylinder;
        qsm(1).branch = branch;
        qsm(1).treedata = treedata;
        qsm(1).rundata.inputs = Inputs;
        qsm(1).rundata.time = single(Time);
        qsm(1).rundata.date = single(Date);
        qsm(1).rundata.version = '2.4.1';
        if inputs.Dist
          qsm(1).pmdistance = pmdis;
        end
        if inputs.Tria
          qsm(1).triangulation = triangulation;
        end
        nmodel = nmodel+1;
        QSM(nmodel) = qsm;
        
        %% Save the output into results-folder
        % matlab-format (.mat)
        if inputs.savemat
          str = [inputs.name,'_t',num2str(inputs.tree),'_m',...
            num2str(inputs.model)];
          save(['results/QSM_',str],'QSM')
        end
        % text-format (.txt)
        if inputs.savetxt
          if nd > 1 || na > 1 || ni > 1
            str = [inputs.name,'_t',num2str(inputs.tree),'_m',...
              num2str(inputs.model)];
            if nd > 1
              str = [str,'_D',num2str(PatchDiam1(h))];
            end
            if na > 1
              str = [str,'_DA',num2str(PatchDiam2Max(i))];
            end
            if ni > 1
              str = [str,'_DI',num2str(PatchDiam2Min(j))];
            end
          else
            str = [inputs.name,'_t',num2str(inputs.tree),'_m',...
              num2str(inputs.model)];
          end
          save_model_text(qsm,str)
        end

        %% Plot models and segmentations
        if inputs.plot >= 1
          if inputs.Tria
            plot_models_segmentations(P,cover2,segment2,cylinder,trunk,...
                triangulation)
          else
            plot_models_segmentations(P,cover2,segment2,cylinder)
          end
          if nd > 1 || na > 1 || ni > 1
            pause
          end
        end
      end
    end
  end
end
