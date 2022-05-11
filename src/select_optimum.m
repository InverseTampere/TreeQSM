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

function [TreeData,OptModels,OptInputs,OptQSM] = ...
  select_optimum(QSMs,Metric,savename)

% ---------------------------------------------------------------------
% SELECT_OPTIMUM.M       Selects optimum models based on point-cylinder model
%                           distances or standard deviations of attributes
%
% Version 1.4.0 
% Latest update     2 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Works for single or multiple tree cases where the input QSMs contains
% multiple models for the same tree with different inputs and multiple runs
% with the same inputs. Allows the user to select from 34 different metrics
% for the optimization. These include average point-model distances from
% all, trunk, branch, 1st-order branch and 2nd-order branch cylinders plus
% some combinations where e.g. "mean trunk and mean branch" or "mean trunk
% and mean 1st-order branch" point-model distances are added together.
% Similarly for the maximum point-model distances and the sums of mean and
% the maximum distances.
%   The difference between "all" and "trunk and branch" is that "all"
% is the average of all cylinder distances which usually emphasizes
% branch cylinder as there usually much more those, whereas "trunk and branch"
% gives equal weight for trunk and branch cylinders.
%   The other options for metric are based on minimizing the standard deviations
% of volumes (total, trunk, branch, trunk+branch which have equal emphasis
% between trunk and branches), lengths (trunk, branches) or total number of
% branches. Here the idea is that if the variance (standard deviation) of
% some attribute between models with the same inputs is small, then it
% indicates some kind of robustness which might indicate that the inputs
% are close to optimal.
%   The optimal single model out of the models with the optimal inputs is
% selected based on the minimum mean point-model-distance.
%
% Inputs:
% QSMs      Contain all the QSMs, possibly from multiple trees
% Metric    Optional input, Metric to be minimized:
%           CYLINDER-DISTANCE METRICS:
%           'all_mean_dis' = mean distances from (mdf) all cylinders, DEFAULT option
%           'trunk_mean_dis' = mdf trunk cylinders,
%           'branch_mean_dis' = mdf all branch cylinders,
%           '1branch_mean_dis' = mdf 1st-order branch cylinders,
%           '2branch_mean_dis' = mdf 2nd-order branch cylinders,
%           'trunk+branch_mean_dis' = mdf trunk + mdf branch cylinders,
%           'trunk+1branch_mean_dis' = mdf trunk + mdf 1st-ord branch cyls,
%           'trunk+1branch+2branch_mean_dis' = above + mdf 2nd-ord branch cyls
%           '1branch+2branch_mean_dis' = mdf 1branch cyls + mdf 2branch cyls
%           'all_max_dis' = maximum distances from (mdf) all cylinders
%           'trunk_max_dis' = mdf trunk cylinders,
%           'branch_max_dis' = mdf all branch cylinders,
%           '1branch_max_dis' = mdf 1st-order branch cylinders,
%           '2branch_max_dis' = mdf 2nd-order branch cylinders,
%           'trunk+branch_max_dis' = mdf trunk + mdf branch cylinders,
%           'trunk+1branch_max_dis' = mdf trunk + mdf 1st-ord branch cyls,
%           'trunk+1branch+2branch_max_dis' = above + mdf 2nd-ord branch cyls.
%           '1branch+2branch_max_dis' = mdf 1branch cyls + mdf 2branch cyls
%           'all_mean+max_dis' = mean + maximum distances from (m+mdf) all cylinders
%           'trunk_mean+max_dis' = (m+mdf) trunk cylinders,
%           'branch_mean+max_dis' = (m+mdf) all branch cylinders,
%           '1branch_mean+max_dis' = (m+mdf) 1st-order branch cylinders,
%           '2branch_mean+max_dis' = (m+mdf) 2nd-order branch cylinders,
%           'trunk+branch_mean+max_dis' = (m+mdf) trunk + (m+mdf) branch cylinders,
%           'trunk+1branch_mean+max_dis' = (m+mdf) trunk + (m+mdf) 1branch cyls,
%           'trunk+1branch+2branch_mean+max_dis' = above + (m+mdf) 2branch cyls.
%           '1branch+2branch_mean+max_dis' = (m+mdf) 1branch cyls + (m+mdf) 2branch cyls
%           STANDARD DEVIATION METRICS:
%           'tot_vol_std' = standard deviation of total volume
%           'trunk_vol_std' = standard deviation of trunk volume
%           'branch_vol_std' = standard deviation of branch volume
%           'trunk+branch_vol_std' = standard deviation of trunk plus branch volume
%           'tot_are_std' = standard deviation of total area
%           'trunk_are_std' = standard deviation of trunk area
%           'branch_are_std' = standard deviation of branch area
%           'trunk+branch_are_std' = standard deviation of trunk plus branch area
%           'trunk_len_std' = standard deviation of trunk length
%           'branch_len_std' = standard deviation of branch length
%           'branch_num_std' = standard deviation of number of branches
%           BRANCH-ORDER DISTRIBUTION METRICS:
%           'branch_vol_ord3_mean' = mean difference in volume of 1-3 branch orders
%           'branch_are_ord3_mean' = mean difference in area of 1-3 branch orders
%           'branch_len_ord3_mean' = mean difference in length of 1-3 branch orders
%           'branch_num_ord3_mean' = mean difference in number of 1-3 branch orders
%           'branch_vol_ord3_max' = max difference in volume of 1-3 branch orders
%           'branch_are_ord3_max' = max difference in area of 1-3 branch orders
%           'branch_len_ord3_max' = max difference in length of 1-3 branch orders
%           'branch_num_ord3_max' = max difference in number of 1-3 branch orders
%           'branch_vol_ord6_mean' = mean difference in volume of 1-6 branch orders
%           'branch_are_ord6_mean' = mean difference in area of 1-6 branch orders
%           'branch_len_ord6_mean' = mean difference in length of 1-6 branch orders
%           'branch_num_ord6_mean' = mean difference in number of 1-6 branch orders
%           'branch_vol_ord6_max' = max difference in volume of 1-6 branch orders
%           'branch_are_ord6_max' = max difference in area of 1-6 branch orders
%           'branch_len_ord6_max' = max difference in length of 1-6 branch orders
%           'branch_num_ord6_max' = max difference in number of 1-6 branch orders
%           CYLINDER DISTRIBUTION METRICS:
%           'cyl_vol_dia10_mean') = mean diff. in volume of 1-10cm diam cyl classes
%           'cyl_are_dia10_mean') = mean diff. in area of 1-10cm diam cyl classes
%           'cyl_len_dia10_mean') = mean diff. in length of 1-10cm diam cyl classes
%           'cyl_vol_dia10_max') = max diff. in volume of 1-10cm diam cyl classes
%           'cyl_are_dia10_max') = max diff. in area of 1-10cm diam cyl classes
%           'cyl_len_dia10_max') = max diff. in length of 1-10cm diam cyl classes
%           'cyl_vol_dia20_mean') = mean diff. in volume of 1-20cm diam cyl classes
%           'cyl_are_dia20_mean') = mean diff. in area of 1-20cm diam cyl classes
%           'cyl_len_dia20_mean') = mean diff. in length of 1-20cm diam cyl classes
%           'cyl_vol_dia20_max') = max diff. in volume of 1-20cm diam cyl classes
%           'cyl_are_dia20_max') = max diff. in area of 1-20cm diam cyl classes
%           'cyl_len_dia20_max') = max diff. in length of 1-20cm diam cyl classes
%           'cyl_vol_zen_mean') = mean diff. in volume of cyl zenith distribution
%           'cyl_are_zen_mean') = mean diff. in area of cyl zenith distribution
%           'cyl_len_zen_mean') = mean diff. in length of cyl zenith distribution
%           'cyl_vol_zen_max') = max diff. in volume of cyl zenith distribution
%           'cyl_are_zen_max') = max diff. in area of cyl zenith distribution
%           'cyl_len_zen_max') = max diff. in length of cyl zenith distribution
%           SURFACE COVERAGE METRICS:
%               metric to be minimized is 1-mean(surface_coverage) or 1-min(SC)
%           'all_mean_surf' = mean surface coverage from (msc) all cylinders
%           'trunk_mean_surf' = msc trunk cylinders,
%           'branch_mean_surf' = msc all branch cylinders,
%           '1branch_mean_surf' = msc 1st-order branch cylinders,
%           '2branch_mean_surf' = msc 2nd-order branch cylinders,
%           'trunk+branch_mean_surf' = msc trunk + msc branch cylinders,
%           'trunk+1branch_mean_surf' = msc trunk + msc 1st-ord branch cyls,
%           'trunk+1branch+2branch_mean_surf' = above + msc 2nd-ord branch cyls
%           '1branch+2branch_mean_surf' = msc 1branch cyls + msc 2branch cyls
%           'all_min_surf' = minimum surface coverage from (msc) all cylinders
%           'trunk_min_surf' = msc trunk cylinders,
%           'branch_min_surf' = msc all branch cylinders,
%           '1branch_min_surf' = msc 1st-order branch cylinders,
%           '2branch_min_surf' = msc 2nd-order branch cylinders,
%           'trunk+branch_min_surf' = msc trunk + msc branch cylinders,
%           'trunk+1branch_min_surf' = msc trunk + msc 1st-ord branch cyls,
%           'trunk+1branch+2branch_min_surf' = above + msc 2nd-ord branch cyls.
%           '1branch+2branch_min_surf' = msc 1branch cyls + msc 2branch cyls
% savename      Optional input, name string specifying the name of the saved file
%               containing the outputs
%
% Outputs:
% TreeData      Similar structure array as the "treedata" in QSMs but now each
%               attribute contains the mean and std computed from the models
%               with the optimal inputs. Also contains the sensitivities
%               for the inputs PatchDiam1, PatchDiam2Min, PatchDiam2Max.
%               Thus for single number attributes (e.g. TotalVolume) there
%               are five numbers [mean std sensi_PD1 sensi_PD2Min sensi_PD2Max]
% OptModels     Indexes of the models with the optimal inputs (column 1) and
%               the index of the optimal single model (column 2) in "QSMs" 
%               for each tree
% OptInputs     The optimal input parameters for each tree
% OptQSMs       The single best QSM for each tree, OptQSMs = QSMs(OptModel);
% ---------------------------------------------------------------------


% Changes from version 1.3.1 to 1.4.0, 2 May 2022:
% 1) Added estimation of (relative) sensitivity of the single number
%    attributes in TreeData for the inputs PatchDiam1, PatchDiam2Min,
%    PatchDiam2Max. Now TreeData contains also these values as the columns
%    3 to 5.
% 2) Corrected a small bug in the subfunction "collect_data" (assignment
%    of values for "CylSurfCov(i,:)"). The bug caused error for QSMs whose
%    maximum branch order is less than 2.
% 3) Bug fix for 3 lines (caused error for some cases and for other cases
%    the optimal single model was wrongly selected):
%    [~,T] = min(dist(ind,best));  -->  [~,T] = min(Data.CylDist(ind,best));

% Changes from version 1.2.0 to 1.3.0, 4 Aug 2020:
% 1) Removed two inputs ("lcyl" and "FilRad") from the inputs to be
%    optimised. This corresponds to changes in the cylinder fitting.
% 2) Added more choices for the optimisation criteria or cost
%    functions ("metric") that are minimised. There is now 91 metrics and
%    the new ones include surface coverage based metrics.

% Changes from version 1.1.1 to 1.2.0, 4 Feb 2020:
% 1) Major change in the structure: subfunctions
% 2) Added more choices for the optimisation criteria or cost
%    functions ("metric") that are minimised. There is now 73 metrics and in
%    particular the new ones include some area related metrics and branch
%    and cylinder distribution based metrics.

% Changes from version 1.1.0 to 1.1.1, 26 Nov 2019:
% 1) Added the "name" of the point cloud from the inputs.name to the output
%    TreeData as a field. Also now displays the name together with the tree
%    number.
% 2) TreeData contains now correctly fields ("location", "StemTaper",
%    "VolumeBranchOrder", etc) from the Optimal QSMs.

% Changes from version 1.0.0 to 1.1.0, 08 Oct 2019:
% 1) Added the posibility to select the optimisation criteria or cost
%    function ("metric") that is minimised from 34 different options.
%    Previously only one option was used. The used metric is also included
%    in "OptInputs" output as one of the fields.
% 2) Added OptQSM as one of the outputs

%% Select the metric based on the input
if nargin > 1
  [met,Metric] = select_metric(Metric);
else
  met = 1;
  Metric = 'all_mean_dis';
end


% The metric for selecting the optimal single model from the models with
% the optimal inputs is the mean point-model-distance.
best = 1;

%% Collect data
% Find the first non-empty model
i = 1;
while isempty(QSMs(i).cylinder)
  i = i+1;
end
% Determine how many single-number attributes there are in treedata
names = fieldnames(QSMs(i).treedata);
n = 1;
while numel(QSMs(i).treedata.(names{n})) == 1
  n = n+1;
end
n = n-1;
Names = names(1:n);
L = max(cellfun('length',Names))+1;
for i = 1:n
  name = Names{i};
  name(L) = ' ';
  Names{i} = name;
end

% Collect data:
[treedata,inputs,TreeId,Data] = collect_data(QSMs,names,n);

% Trees and their unique IDs
TreeIds = unique(TreeId(:,1));
nt = length(TreeIds); % number of trees

DataM = zeros(n,nt);
DataS = zeros(n,nt); % Standard deviation of tree data for each tree
DataM2 = DataM;     DataM3 = DataM;
DataS2 = DataS;     DataS3 = DataS;

OptIn = zeros(nt,9); % Optimal input values
OptDist = zeros(nt,9); % Smallest metric values


% average treedata and inputs for each tree-input-combination:
TreeDataAll = zeros(nt,5*5*5,n);
Inputs = zeros(nt,5*5*5,3);

IndAll = (1:1:size(TreeId,1))';

% Indexes of the optimal single models in QSMs:
OptModel = zeros(nt,3);
% The indexes of models in QSMs with the optimal inputs (col 1)
% and the indexes of the optimal single models (col 2):
OptModels = cell(nt,2);

NInputs = zeros(nt,1);

%% Process each tree separately
for tree = 1:nt
  % Select the models for the tree
  Models = TreeId(:,1) == TreeIds(tree);

  %% Determine the input parameter values
  InputParComb = unique(inputs(Models,:),'rows'); % Input parameter combinations
  IV = cell(3,1);
  N = zeros(3,1);
  for i = 1:3
    I = unique(InputParComb(:,i));
    IV{i} = I;
    N(i) = length(I);
  end

  %% Determine metric-value for each input
  % (average over number of models with the same inputs)
  input = cell(1,N(1)*N(2)*N(3));
  distM = zeros(1,N(1)*N(2)*N(3)); % average distances or volume stds
  b = 0;
  for d = 1:N(1) % PatchDiam1
    J = abs(inputs(:,1)-IV{1}(d)) < 0.0001;
    for a = 1:N(2) % PatchDiam2Min
      K = abs(inputs(:,2)-IV{2}(a)) < 0.0001;
      for i = 1:N(3) % PatchDiam2Max
        L = abs(inputs(:,3)-IV{3}(i)) < 0.0001;

        % Select models for the tree with the same inputs:
        T = Models & J & K & L;
        b = b+1;
        input{b} = [d a i];

        % Compute the metric value;
        D = compute_metric_value(met,T,treedata,Data);
        distM(b) = D;

        % Collect the data and inputs
        TreeDataAll(tree,b,:) = mean(treedata(:,T),2);
        Inputs(tree,b,:) = [IV{1}(d) IV{2}(a) IV{3}(i)];
      end
    end
  end

  %% Determine the optimal inputs and models
  ninputs = prod(N);
  NInputs(tree) = ninputs;
  [d,J] = sort(distM);
  O = input{J(1)};
  OptIn(tree,1:3) = [IV{1}(O(1)) IV{2}(O(2)) IV{3}(O(3))];
  OptDist(tree,1) = d(1);
  if ninputs > 1
    O = input{J(2)};
    OptIn(tree,4:6) = [IV{1}(O(1)) IV{2}(O(2)) IV{3}(O(3))];
    OptDist(tree,2) = d(2);
    if ninputs > 2
      O = input{J(3)};
      OptIn(tree,7:9) = [IV{1}(O(1)) IV{2}(O(2)) IV{3}(O(3))];
      OptDist(tree,3) = d(3);
    end
  end

  %% Mean of tree data for each tree computed from the optimal models:
  % Select the optimal models for each tree: In the case of multiple models
  % with same inputs, select the one model with the optimal inputs that
  % has the minimum metric value.
  J = abs(inputs(:,1)-OptIn(tree,1)) < 0.0001;
  K = abs(inputs(:,2)-OptIn(tree,2)) < 0.0001;
  L = abs(inputs(:,3)-OptIn(tree,3)) < 0.0001;
  T = Models & J & K & L;
  ind = IndAll(T);
  [~,T] = min(Data.CylDist(ind,best));
  OptModel(tree,1) = ind(T);
  OptModels{tree,1} = ind;
  OptModels{tree,2} = ind(T);
  DataM(:,tree) = mean(treedata(:,ind),2);
  DataS(:,tree) = std(treedata(:,ind),[],2);
  if ninputs > 1
    J = abs(inputs(:,1)-OptIn(tree,4)) < 0.0001;
    K = abs(inputs(:,2)-OptIn(tree,5)) < 0.0001;
    L = abs(inputs(:,3)-OptIn(tree,6)) < 0.0001;
    T = Models & J & K & L;
    ind = IndAll(T);
    [~,T] = min(Data.CylDist(ind,best));
    OptModel(tree,2) = ind(T);
    DataM2(:,tree) = mean(treedata(:,ind),2);
    DataS2(:,tree) = std(treedata(:,ind),[],2);
    if ninputs > 2
      J = abs(inputs(:,1)-OptIn(tree,7)) < 0.0001;
      K = abs(inputs(:,2)-OptIn(tree,8)) < 0.0001;
      L = abs(inputs(:,3)-OptIn(tree,9)) < 0.0001;
      T = Models & J & K & L;
      ind = IndAll(T);
      [~,T] = min(Data.CylDist(ind,best));
      OptModel(tree,3) = ind(T);
      DataM3(:,tree) = mean(treedata(:,ind),2);
      DataS3(:,tree) = std(treedata(:,ind),[],2);
    end
  end

  % Decrease the number on non-zero decimals
  DataM(:,tree) = change_precision(DataM(:,tree));
  DataS(:,tree) = change_precision(DataS(:,tree));
  if ninputs > 1
    DataM2(:,tree) = change_precision(DataM2(:,tree));
    DataS2(:,tree) = change_precision(DataS2(:,tree));
    if ninputs > 2
      DataM3(:,tree) = change_precision(DataM3(:,tree));
      DataS3(:,tree) = change_precision(DataS3(:,tree));
    end
  end

  % Define the output "OptInputs"
  OptM = IndAll(OptModel(tree,1));
  OptInputs(tree) = QSMs(OptM).rundata.inputs;
  if ninputs > 1
    OptM2 = IndAll(OptModel(tree,2));
    OI2(tree) = QSMs(OptM2).rundata.inputs;
    if ninputs > 2
      OptM3 = IndAll(OptModel(tree,3));
      OI3(tree) = QSMs(OptM3).rundata.inputs;
    end
  end

end
N = max(NInputs);
TreeDataAll = TreeDataAll(:,1:N,:);
Inputs = Inputs(:,1:N,:);

% Compute Coefficient of variation for the data
OptModel = IndAll(OptModel(:,1));
OptQSM = QSMs(OptModel);
DataCV = DataS./DataM*100; % Coefficient of variation
if ninputs > 1
  DataCV2 = DataS2./DataM2*100; % Coefficient of variation
  if ninputs > 2
    DataCV3 = DataS3./DataM3*100; % Coefficient of variation
  end
end
% Decrease the number on non-zero decimals
for j = 1:nt
  DataCV(:,j) = change_precision(DataCV(:,j));
  if ninputs > 1
    DataCV2(:,j) = change_precision(DataCV2(:,j));
    if ninputs > 2
      DataCV3(:,j) = change_precision(DataCV3(:,j));
    end
  end
end

%% Display some data about optimal models
% Display optimal inputs, model and attributes for each tree
for t = 1:nt
  disp('-------------------------------')
  disp(['  Tree: ',num2str(OptInputs(t).tree),', ',OptInputs(t).name])
  if NInputs(t) == 1
    disp(['    Metric: ',Metric])
    disp(['    Metric value:  ',num2str(1000*OptDist(t,1))])
    disp(['    Optimal inputs:  PatchDiam1 = ',...
      num2str(OptInputs(t).PatchDiam1)])
    disp(['                  PatchDiam2Min = ',...
      num2str(OptInputs(t).PatchDiam2Min)])
    disp(['                  PatchDiam2Max = ',...
      num2str(OptInputs(t).PatchDiam2Max)])
    disp(['    Optimal model: ',num2str(OptModel(t))])
    sec = num2str(round(QSMs(OptModel(t)).rundata.time(end)));
    disp(['    Reconstruction time for the optimal model: ',sec,' seconds'])
    disp('    Attributes (mean, std, CV(%)):')
    for i = 1:n
      str = (['      ',Names{i},': ',num2str([...
        DataM(i,t) DataS(i,t) DataCV(i,t)])]);
      disp(str)
    end
  elseif NInputs(t) == 2
    disp('    The best two cases:')
    disp(['    Metric: ',Metric])
    disp(['    Metric values:  ',num2str(OptDist(t,1:2))])
    disp(['            inputs:  PatchDiam1 = ',...
      num2str([OptInputs(t).PatchDiam1 OI2(t).PatchDiam1])])
    disp(['                  PatchDiam2Min = ',...
      num2str([OptInputs(t).PatchDiam2Min OI2(t).PatchDiam2Min])])
    disp(['                  PatchDiam2Max = ',...
      num2str([OptInputs(t).PatchDiam2Max OI2(t).PatchDiam2Max])])
    disp(['    Optimal model: ',num2str(OptModel(t))])
    sec = num2str(round(QSMs(OptModel(t)).rundata.time(end)));
    disp(['    Reconstruction time for the optimal model: ',sec,' seconds'])
    disp('    Attributes (mean, std, CV(%), second best mean):')
    for i = 1:n
      str = (['      ',Names{i},':  ',num2str([DataM(i,t) ...
        DataS(i,t) DataCV(i,t) DataM2(i,t)])]);
      disp(str)
    end
  elseif NInputs(t) > 2
    disp('    The best three cases:')
    disp(['    Metric: ',Metric])
    disp(['    Metric values:  ',num2str(OptDist(t,1:3))])
    disp(['            inputs:  PatchDiam1 = ',num2str([...
      OptInputs(t).PatchDiam1 OI2(t).PatchDiam1 OI3(t).PatchDiam1])])
    disp(['                  PatchDiam2Min = ',num2str([...
      OptInputs(t).PatchDiam2Min OI2(t).PatchDiam2Min OI3(t).PatchDiam2Min])])
    disp(['                  PatchDiam2Max = ',num2str([...
      OptInputs(t).PatchDiam2Max OI2(t).PatchDiam2Max OI3(t).PatchDiam2Max])])
    disp(['    Optimal model: ',num2str(OptModel(t))])
    sec = num2str(round(QSMs(OptModel(t)).rundata.time(end)));
    disp(['    Reconstruction time for the optimal model: ',sec,' seconds'])
    str = ['    Attributes (mean, std, CV(%),',...
      ' second best mean, third best mean, sensitivity):'];
    disp(str)
    for i = 1:n
      sensi = max(abs([DataM(i,t)-DataM2(i,t)...
        DataM(i,t)-DataM3(i,t)])/DataM(i,t));
      sensi2 = 100*sensi;
      sensi = 100*sensi/DataCV(i,t);
      sensi2 = change_precision(sensi2);
      sensi = change_precision(sensi);
      str = (['      ',Names{i},':  ',num2str([DataM(i,t) DataS(i,t) ...
        DataCV(i,t) DataM2(i,t) DataM3(i,t) sensi sensi2])]);
      disp(str)
    end
  end
  disp('------')
end

%% Compute the sensitivity of the tree attributes relative to PatchDiam-parameters
Sensi = sensitivity_analysis(TreeDataAll,TreeId,Inputs,OptIn,NInputs);

%% Generate TreeData sructure for optimal models
clear TreeData
TreeData = vertcat(OptQSM(:).treedata);
for t = 1:nt
  for i = 1:n
    TreeData(t).(names{i}) = [DataM(i,t) DataS(i,t) squeeze(Sensi(t,i,:))'];
  end
  TreeData(t).name = OptInputs(t).name;
end

%% Add the metric for the "OptInputs"
for i = 1:nt
  OptInputs(i).metric = Metric;
end

%% Save results
if nargin == 3
  str = ['results/OptimalQSMs_',savename];
  save(str,'TreeData','OptModels','OptInputs','OptQSM')

  str = ['results/tree_data_',savename,'.txt'];
  fid = fopen(str, 'wt');
  fprintf(fid, [repmat('%g\t', 1, size(DataM,2)-1) '%g\n'], DataM.');
  fclose(fid);
end

% End of main function
end


function [treedata,inputs,TreeId,Data] = collect_data(...
  QSMs,names,Nattri)

Nmod = max(size(QSMs)); % number of models
treedata = zeros(Nattri,Nmod); % Collect all tree attributes from all models
inputs = zeros(Nmod,3); % collect the inputs from all models
% ([PatchDiam1 PatchDiam2Min PatchDiam2Max])
CylDist = zeros(Nmod,10); % collect the distances from all models
CylSurfCov = zeros(Nmod,10); % collect the surface coverages from all models
s = 6; % maximum branch order
OrdDis = zeros(Nmod,4*s); % collect the distributions from all the models
r = 20; % maximum cylinder diameter
CylDiaDis = zeros(Nmod,3*r);
CylZenDis = zeros(Nmod,54);
TreeId = zeros(Nmod,2); % collectd the tree and model indexes from all models
Keep = true(Nmod,1); % Non-empty models

for i = 1:Nmod
  if ~isempty(QSMs(i).cylinder)
    % Collect input-parameter values and tree IDs:
    p = QSMs(i).rundata.inputs;
    inputs(i,:) = [p.PatchDiam1 p.PatchDiam2Min p.PatchDiam2Max];
    TreeId(i,:) = [p.tree p.model];

    % Collect cylinder-point distances: mean of all cylinders,
    % mean of trunk, branch, 1st- and 2nd-order branch cylinders.
    % And the maximum of the previous:
    D = QSMs(i).pmdistance;
    CylDist(i,:) =  [D.mean  D.TrunkMean  D.BranchMean  D.Branch1Mean ...
      D.Branch2Mean D.max  D.TrunkMax  D.BranchMax  D.Branch1Max ...
      D.Branch2Max];

    % Collect surface coverages: mean of all cylinders,
    % mean of trunk, branch, 1st- and 2nd-order branch cylinders.
    % And the minimum of the previous:
    D = QSMs(i).cylinder.SurfCov;
    T = QSMs(i).cylinder.branch == 1;
    B1 = QSMs(i).cylinder.BranchOrder == 1;
    B2 = QSMs(i).cylinder.BranchOrder == 2;
    if ~any(B1)
      CylSurfCov(i,:) =  [mean(D)  mean(D(T))  0 0 0 ...
         min(D) min(D(T)) 0 0 0];
    elseif ~any(B2)
      CylSurfCov(i,:) =  [mean(D)  mean(D(T))  mean(D(~T))  mean(D(B1)) ...
        0 min(D) min(D(T)) min(D(~T)) min(D(B1)) 0];
    else
      CylSurfCov(i,:) =  [mean(D)  mean(D(T))  mean(D(~T))  mean(D(B1)) ...
        mean(D(B2)) min(D) min(D(T)) min(D(~T)) min(D(B1)) min(D(B2))];
    end

    % Collect branch-order distributions:
    d = QSMs(i).treedata.VolBranchOrd;
    nd = length(d);
    if nd > 0
      a = min(nd,s);
      OrdDis(i,1:a) = d(1:a);
      OrdDis(i,s+1:s+a) = QSMs(i).treedata.AreBranchOrd(1:a);
      OrdDis(i,2*s+1:2*s+a) = QSMs(i).treedata.LenBranchOrd(1:a);
      OrdDis(i,3*s+1:3*s+a) = QSMs(i).treedata.NumBranchOrd(1:a);
    end

    % Collect cylinder diameter distributions:
    d = QSMs(i).treedata.VolCylDia;
    nd = length(d);
    if nd > 0
      a = min(nd,r);
      CylDiaDis(i,1:a) = d(1:a);
      CylDiaDis(i,r+1:r+a) = QSMs(i).treedata.AreCylDia(1:a);
      CylDiaDis(i,2*r+1:2*r+a) = QSMs(i).treedata.LenCylDia(1:a);
    end

    % Collect cylinder zenith direction distributions:
    d = QSMs(i).treedata.VolCylZen;
    if ~isempty(d)
      CylZenDis(i,1:18) = d;
      CylZenDis(i,19:36) = QSMs(i).treedata.AreCylZen;
      CylZenDis(i,37:54) = QSMs(i).treedata.LenCylZen;
    end

    % Collect the treedata values from each model
    for j = 1:Nattri
      treedata(j,i) = QSMs(i).treedata.(names{j});
    end

  else
    Keep(i) = false;
  end
end
treedata = treedata(:,Keep);
inputs = inputs(Keep,:);
TreeId = TreeId(Keep,:);
clear Data
Data.CylDist = CylDist(Keep,:);
Data.CylSurfCov = CylSurfCov(Keep,:);
Data.BranchOrdDis = OrdDis(Keep,:);
Data.CylDiaDis = CylDiaDis(Keep,:);
Data.CylZenDis = CylZenDis(Keep,:);

% End of function
end


function [met,Metric] = select_metric(Metric)

% Mean distance metrics:
if strcmp(Metric,'all_mean_dis')
  met = 1;
elseif strcmp(Metric,'trunk_mean_dis')
  met = 2;
elseif strcmp(Metric,'branch_mean_dis')
  met = 3;
elseif strcmp(Metric,'1branch_mean_dis')
  met = 4;
elseif strcmp(Metric,'2branch_mean_dis')
  met = 5;
elseif strcmp(Metric,'trunk+branch_mean_dis')
  met = 6;
elseif strcmp(Metric,'trunk+1branch_mean_dis')
  met = 7;
elseif strcmp(Metric,'trunk+1branch+2branch_mean_dis')
  met = 8;
elseif strcmp(Metric,'1branch+2branch_mean_dis')
  met = 9;

  % Maximum distance metrics:
elseif strcmp(Metric,'all_max_dis')
  met = 10;
elseif strcmp(Metric,'trunk_max_dis')
  met = 11;
elseif strcmp(Metric,'branch_max_dis')
  met = 12;
elseif strcmp(Metric,'1branch_max_dis')
  met = 13;
elseif strcmp(Metric,'2branch_max_dis')
  met = 14;
elseif strcmp(Metric,'trunk+branch_max_dis')
  met = 15;
elseif strcmp(Metric,'trunk+1branch_max_dis')
  met = 16;
elseif strcmp(Metric,'trunk+1branch+2branch_max_dis')
  met = 17;
elseif strcmp(Metric,'1branch+2branch_max_dis')
  met = 18;

  % Mean plus Maximum distance metrics:
elseif strcmp(Metric,'all_mean+max_dis')
  met = 19;
elseif strcmp(Metric,'trunk_mean+max_dis')
  met = 20;
elseif strcmp(Metric,'branch_mean+max_dis')
  met = 21;
elseif strcmp(Metric,'1branch_mean+max_dis')
  met = 22;
elseif strcmp(Metric,'2branch_mean+max_dis')
  met = 23;
elseif strcmp(Metric,'trunk+branch_mean+max_dis')
  met = 24;
elseif strcmp(Metric,'trunk+1branch_mean+max_dis')
  met = 25;
elseif strcmp(Metric,'trunk+1branch+2branch_mean+max_dis')
  met = 26;
elseif strcmp(Metric,'1branch+2branch_mean+max_dis')
  met = 27;

  % Standard deviation metrics:
elseif strcmp(Metric,'tot_vol_std')
  met = 28;
elseif strcmp(Metric,'trunk_vol_std')
  met = 29;
elseif strcmp(Metric,'branch_vol_std')
  met = 30;
elseif strcmp(Metric,'trunk+branch_vol_std')
  met = 31;
elseif strcmp(Metric,'tot_are_std')
  met = 32;
elseif strcmp(Metric,'trunk_are_std')
  met = 33;
elseif strcmp(Metric,'branch_are_std')
  met = 34;
elseif strcmp(Metric,'trunk+branch_are_std')
  met = 35;
elseif strcmp(Metric,'trunk_len_std')
  met = 36;
elseif strcmp(Metric,'trunk+branch_len_std')
  met = 37;
elseif strcmp(Metric,'branch_len_std')
  met = 38;
elseif strcmp(Metric,'branch_num_std')
  met = 39;

  % Branch order distribution metrics:
elseif strcmp(Metric,'branch_vol_ord3_mean')
  met = 40;
elseif strcmp(Metric,'branch_are_ord3_mean')
  met = 41;
elseif strcmp(Metric,'branch_len_ord3_mean')
  met = 42;
elseif strcmp(Metric,'branch_num_ord3_mean')
  met = 43;
elseif strcmp(Metric,'branch_vol_ord3_max')
  met = 44;
elseif strcmp(Metric,'branch_are_ord3_max')
  met = 45;
elseif strcmp(Metric,'branch_len_ord3_max')
  met = 46;
elseif strcmp(Metric,'branch_num_ord3_max')
  met = 47;
elseif strcmp(Metric,'branch_vol_ord6_mean')
  met = 48;
elseif strcmp(Metric,'branch_are_ord6_mean')
  met = 49;
elseif strcmp(Metric,'branch_len_ord6_mean')
  met = 50;
elseif strcmp(Metric,'branch_num_ord6_mean')
  met = 51;
elseif strcmp(Metric,'branch_vol_ord6_max')
  met = 52;
elseif strcmp(Metric,'branch_are_ord6_max')
  met = 53;
elseif strcmp(Metric,'branch_len_ord6_max')
  met = 54;
elseif strcmp(Metric,'branch_num_ord6_max')
  met = 55;

  % Cylinder distribution metrics:
elseif strcmp(Metric,'cyl_vol_dia10_mean')
  met = 56;
elseif strcmp(Metric,'cyl_are_dia10_mean')
  met = 57;
elseif strcmp(Metric,'cyl_len_dia10_mean')
  met = 58;
elseif strcmp(Metric,'cyl_vol_dia10_max')
  met = 59;
elseif strcmp(Metric,'cyl_are_dia10_max')
  met = 60;
elseif strcmp(Metric,'cyl_len_dia10_max')
  met = 61;
elseif strcmp(Metric,'cyl_vol_dia20_mean')
  met = 62;
elseif strcmp(Metric,'cyl_are_dia20_mean')
  met = 63;
elseif strcmp(Metric,'cyl_len_dia20_mean')
  met = 64;
elseif strcmp(Metric,'cyl_vol_dia20_max')
  met = 65;
elseif strcmp(Metric,'cyl_are_dia20_max')
  met = 66;
elseif strcmp(Metric,'cyl_len_dia20_max')
  met = 67;
elseif strcmp(Metric,'cyl_vol_zen_mean')
  met = 68;
elseif strcmp(Metric,'cyl_are_zen_mean')
  met = 69;
elseif strcmp(Metric,'cyl_len_zen_mean')
  met = 70;
elseif strcmp(Metric,'cyl_vol_zen_max')
  met = 71;
elseif strcmp(Metric,'cyl_are_zen_max')
  met = 72;
elseif strcmp(Metric,'cyl_len_zen_max')
  met = 73;

  % Mean surface coverage metrics:
elseif strcmp(Metric,'all_mean_surf')
  met = 74;
elseif strcmp(Metric,'trunk_mean_surf')
  met = 75;
elseif strcmp(Metric,'branch_mean_surf')
  met = 76;
elseif strcmp(Metric,'1branch_mean_surf')
  met = 77;
elseif strcmp(Metric,'2branch_mean_surf')
  met = 78;
elseif strcmp(Metric,'trunk+branch_mean_surf')
  met = 79;
elseif strcmp(Metric,'trunk+1branch_mean_surf')
  met = 80;
elseif strcmp(Metric,'trunk+1branch+2branch_mean_surf')
  met = 81;
elseif strcmp(Metric,'1branch+2branch_mean_surf')
  met = 82;

  % Minimum surface coverage metrics:
elseif strcmp(Metric,'all_min_surf')
  met = 83;
elseif strcmp(Metric,'trunk_min_surf')
  met = 84;
elseif strcmp(Metric,'branch_min_surf')
  met = 85;
elseif strcmp(Metric,'1branch_min_surf')
  met = 86;
elseif strcmp(Metric,'2branch_min_surf')
  met = 87;
elseif strcmp(Metric,'trunk+branch_min_surf')
  met = 88;
elseif strcmp(Metric,'trunk+1branch_min_surf')
  met = 89;
elseif strcmp(Metric,'trunk+1branch+2branch_min_surf')
  met = 90;
elseif strcmp(Metric,'1branch+2branch_min_surf')
  met = 91;

  % Not given in right form, take the default option
else
  met = 1;
  Metric = 'all_mean_dis';
end
% End of function
end


function D = compute_metric_value(met,T,treedata,Data)


if met <= 27 % cylinder distance metrics:
  D = mean(Data.CylDist(T,:),1);
  D(6:10) = 0.5*D(6:10); % Half the maximum values
end

if met < 10 % mean cylinder distance metrics:
  if met == 1 % all_mean_dis
    D = D(1);
  elseif met == 2 % trunk_mean_dis
    D = D(2);
  elseif met == 3 % branch_mean_dis
    D = D(3);
  elseif met == 4 % 1branch_mean_dis
    D = D(4);
  elseif met == 5 % 2branch_mean_dis
    D = D(5);
  elseif met == 6 % trunk+branch_mean_dis
    D = D(2)+D(3);
  elseif met == 7 % trunk+1branch_mean_dis
    D = D(2)+D(4);
  elseif met == 8 % trunk+1branch+2branch_mean_dis
    D = D(2)+D(4)+D(5);
  elseif met == 9 % 1branch+2branch_mean_dis
    D = D(4)+D(5);
  end

elseif met < 19 % maximum cylinder distance metrics:
  if met == 10 % all_max_dis
    D = D(6);
  elseif met == 11 % trunk_max_dis
    D = D(7);
  elseif met == 12 % branch_max_dis
    D = D(8);
  elseif met == 13 % 1branch_max_dis
    D = D(9);
  elseif met == 14 % 2branch_max_dis
    D = D(10);
  elseif met == 15 % trunk+branch_max_dis
    D = D(7)+D(8);
  elseif met == 16 % trunk+1branch_max_dis
    D = D(7)+D(9);
  elseif met == 17 % trunk+1branch+2branch_max_dis
    D = D(7)+D(9)+D(10);
  elseif met == 18 % 1branch+2branch_max_dis
    D = D(9)+D(10);
  end

elseif met < 28 % Mean plus maximum cylinder distance metrics:
  if met == 19 % all_mean+max_dis
    D = D(1)+D(6);
  elseif met == 20 % trunk_mean+max_dis
    D = D(2)+D(7);
  elseif met == 21 % branch_mean+max_dis
    D = D(3)+D(8);
  elseif met == 22 % 1branch_mean+max_dis
    D = D(4)+D(9);
  elseif met == 23 % 2branch_mean+max_dis
    D = D(5)+D(10);
  elseif met == 24 % trunk+branch_mean+max_dis
    D = D(2)+D(3)+D(7)+D(8);
  elseif met == 25 % trunk+1branch_mean+max_dis
    D = D(2)+D(4)+D(7)+D(9);
  elseif met == 26 % trunk+1branch+2branch_mean+max_dis
    D = D(2)+D(4)+D(5)+D(7)+D(9)+D(10);
  elseif met == 27 % 1branch+2branch_mean+max_dis
    D = D(4)+D(5)+D(9)+D(10);
  end

elseif met < 39 % Standard deviation metrics:
  if met == 28 % tot_vol_std
    D = std(treedata(1,T));
  elseif met == 29 % trunk_vol_std
    D = std(treedata(2,T));
  elseif met == 30 % branch_vol_std
    D = std(treedata(3,T));
  elseif met == 31 % trunk+branch_vol_std
    D = std(treedata(2,T))+std(treedata(3,T));
  elseif met == 32 % tot_are_std
    D = std(treedata(12,T));
  elseif met == 33 % trunk_are_std
    D = std(treedata(10,T));
  elseif met == 34 % branch_are_std
    D = std(treedata(11,T));
  elseif met == 35 % trunk+branch_are_std
    D = std(treedata(10,T))+std(treedata(11,T));
  elseif met == 36 % trunk_len_std
    D = std(treedata(5,T));
  elseif met == 37 % branch_len_std
    D = std(treedata(6,T));
  elseif met == 38 % trunk+branch_len_std
    D = std(treedata(5,T))+std(treedata(6,T));
  elseif met == 39 % branch_num_std
    D = std(treedata(8,T));
  end

elseif met < 56 % Branch order metrics:
  dis = max(Data.BranchOrdDis(T,:),[],1)-min(Data.BranchOrdDis(T,:),[],1);
  M = mean(Data.BranchOrdDis(T,:),1);
  I = M > 0;
  dis(I) = dis(I)./M(I);
  if met == 40 % branch_vol_ord3_mean
    D = mean(dis(1:3));
  elseif met == 41 % branch_are_ord3_mean
    D = mean(dis(7:9));
  elseif met == 42 % branch_len_ord3_mean
    D = mean(dis(13:15));
  elseif met == 43 % branch_num_ord3_mean
    D = mean(dis(19:21));
  elseif met == 44 % branch_vol_ord3_max
    D = max(dis(1:3));
  elseif met == 45 % branch_are_ord3_max
    D = max(dis(7:9));
  elseif met == 46 % branch_len_ord3_max
    D = max(dis(13:15));
  elseif met == 47 % branch_vol_ord3_max
    D = max(dis(19:21));
  elseif met == 48 % branch_vol_ord6_mean
    D = mean(dis(1:6));
  elseif met == 49 % branch_are_ord6_mean
    D = mean(dis(7:12));
  elseif met == 50 % branch_len_ord6_mean
    D = mean(dis(13:18));
  elseif met == 51 % branch_num_ord6_mean
    D = mean(dis(19:24));
  elseif met == 52 % branch_vol_ord6_max
    D = max(dis(1:6));
  elseif met == 53 % branch_are_ord6_max
    D = max(dis(7:12));
  elseif met == 54 % branch_len_ord6_max
    D = max(dis(13:18));
  elseif met == 55 % branch_vol_ord6_max
    D = max(dis(19:24));
  end

elseif met < 68 % Cylinder diameter distribution metrics:
  dis = max(Data.CylDiaDis(T,:),[],1)-min(Data.CylDiaDis(T,:),[],1);
  M = mean(Data.CylDiaDis(T,:),1);
  I = M > 0;
  dis(I) = dis(I)./M(I);
  if met == 56 % cyl_vol_dia10_mean
    D = mean(dis(1:10));
  elseif met == 57 % cyl_are_dia10_mean
    D = mean(dis(21:30));
  elseif met == 58 % cyl_len_dia10_mean
    D = mean(dis(41:50));
  elseif met == 59 % cyl_vol_dia10_max
    D = max(dis(1:10));
  elseif met == 60 % cyl_are_dia10_max
    D = max(dis(21:30));
  elseif met == 61 % cyl_len_dia10_max
    D = max(dis(41:50));
  elseif met == 62 % cyl_vol_dia20_mean
    D = mean(dis(1:20));
  elseif met == 63 % cyl_are_dia20_mean
    D = mean(dis(21:40));
  elseif met == 64 % cyl_len_dia20_mean
    D = mean(dis(41:60));
  elseif met == 65 % cyl_vol_dia20_max
    D = max(dis(1:20));
  elseif met == 66 % cyl_are_dia20_max
    D = max(dis(21:40));
  elseif met == 67 % cyl_len_dia20_max
    D = max(dis(41:60));
  end

elseif met < 74 % Cylinder zenith distribution metrics:
  dis = max(Data.CylZenDis(T,:),[],1)-min(Data.CylZenDis(T,:),[],1);
  M = mean(Data.CylZenDis(T,:),1);
  I = M > 0;
  dis(I) = dis(I)./M(I);
  if met == 68 % cyl_vol_zen_mean
    D = mean(dis(1:18));
  elseif met == 69 % cyl_are_zen_mean
    D = mean(dis(19:36));
  elseif met == 70 % cyl_len_zen_mean
    D = mean(dis(37:54));
  elseif met == 71 % cyl_vol_zen_max
    D = max(dis(1:18));
  elseif met == 72 % cyl_are_zen_max
    D = max(dis(19:36));
  elseif met == 73 % cyl_len_zen_max
    D = max(dis(37:54));
  end

elseif met < 92 % Surface coverage metrics:
  D = 1-mean(Data.CylSurfCov(T,:),1);
  if met == 74 % all_mean_surf
    D = D(1);
  elseif met == 75 % trunk_mean_surf
    D = D(2);
  elseif met == 76 % branch_mean_surf
    D = D(3);
  elseif met == 77 % 1branch_mean_surf
    D = D(4);
  elseif met == 78 % 2branch_mean_surf
    D = D(5);
  elseif met == 79 % trunk+branch_mean_surf
    D = D(2)+D(3);
  elseif met == 80 % trunk+1branch_mean_surf
    D = D(2)+D(4);
  elseif met == 81 % trunk+1branch+2branch_mean_surf
    D = D(2)+D(4)+D(5);
  elseif met == 82 % 1branch+2branch_mean_surf
    D = D(4)+D(5);
  elseif met == 83 % all_min_surf
    D = D(6);
  elseif met == 84 % trunk_min_surf
    D = D(7);
  elseif met == 85 % branch_min_surf
    D = D(8);
  elseif met == 86 % 1branch_min_surf
    D = D(9);
  elseif met == 87 % 2branch_min_surf
    D = D(10);
  elseif met == 88 % trunk+branch_min_surf
    D = D(6)+D(7);
  elseif met == 89 % trunk+1branch_min_surf
    D = D(6)+D(8);
  elseif met == 90 % trunk+1branch+2branch_min_surf
    D = D(6)+D(9)+D(10);
  elseif met == 91 % 1branch+2branch_min_surf
    D = D(9)+D(10);
  end
end
% End of function
end


function Sensi = sensitivity_analysis(TreeDataAll,TreeId,Inputs,OptIn,NInputs)

% Computes the sensitivity of tree attributes (e.g. total volume) to the
% changes of input parameter, the PatchDiam parameters, values. The
% sensitivity is normalized, i.e. the relative change of attribute value
% (= max change in attribute value divided by the value with the optimal
% inputs) is divided by the relative change of input parameter value. The
% sensitivity is also expressed as percentage, i.e. multiplied by 100. The
% sensitivity is computed relative PatchDiam1, PatchDiam2Min, and
% PatchDiam2Max. The sensitivity is computed only from the attributes with
% the input parameter values the closest to the optimal value. This way we
% get the local sensitivity in the neighborhood of the optimal input.
%
% Output:
% Sensi       3D-array (#trees,#attributes,#inputs)

TreeIds = unique(TreeId(:,1)); % Unique tree IDs
nt = length(TreeIds); % number of trees
A = [2 3; 1 3; 1 2]; % Keep other two inputs constant and let one varie
Sensi = zeros(nt,size(TreeDataAll,3),3); % initialization of the output
for t = 1:nt % trees
  if NInputs(t) > 1
    D = squeeze(TreeDataAll(t,1:NInputs(t),:))'; % Select the attributes for the tree
    In = squeeze(Inputs(t,1:NInputs(t),:)); % Select the inputs for the tree
    n = size(In,1); % number of different input-combinations
    I = all(In == OptIn(t,1:3),2); % Which data are with the optimal inputs
    ind = (1:1:n)';
    I = ind(I);
    for i = 1:3 % inputs
      if length(unique(In(:,i))) > 1
        dI = abs(max(In(:,i),[],2)-OptIn(t,i));
        dImin = min(dI(dI > 0)); % the minimum nonzero absolute change in inputs
        dI = dImin/OptIn(t,i); % relative change in the attributes
        K1 = abs(max(In(:,i),[],2)-min(OptIn(t,i),[],2)) < dImin+0.0001;
        K = K1 & abs(max(In(:,i),[],2)-min(OptIn(t,i),[],2)) > 0.0001;
        K = ind(K); % the inputs the closest to the optimal input
        J = all(In(K,A(i,:)) == OptIn(t,A(i,:)),2);
        J = K(J); % input i the closest to the optimal and the other two equal the optimal
        dD = max(abs(D(:,J)-D(:,I)),[],2);
        dD = dD./D(:,I); % relative change in the input
        d = dD/dI*100; % relative sensitivity as a percentage
        Sensi(t,:,i) = round(100*d)/100;
      end
    end
  end
end
% End of function
end
