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

function [TreeData,OptModels,OptInputs,OptQSM] = select_optimum(QSMs,Metric,savename)

% ---------------------------------------------------------------------
% SELECT_OPTIMUM.M       Selects optimum models based on point-cylinder model 
%                           distances or standard deviations of attributes
%
% Version 1.1.0
% Latest update     08 Oct 2019
%
% Copyright (C) 2013-2019 Pasi Raumonen
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
%   The difference between "all" and "trunk and branch" 
% is that "all" is the average of all cylinder distances which usually 
% emphasizes branch cylinder as there usually much more those, whereas 
% "trunk and branch" gives equal weight for trunk and branch cylinders. 
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
%           'tot_vol' = standard deviation of total volume
%           'trunk_vol' = standard deviation of trunk volume
%           'branch_vol' = standard deviation of branch volume
%           'trunk+branch_vol' = standard deviation of trunk plus branch volume
%           'trunk_len' = standard deviation of trunk length
%           'branch_len' = standard deviation of branch length
%           'branch_num' = standard deviation of number of branches
% savename      Optional input, name string specifying the name of the saved file
%               containing the outputs
% Outputs:
% TreeData      Similar structure array as the "treedata" in QSMs but now each
%               attribute contains the mean and std computed from the
%               models with the optimal inputs. 
% OptModels     Indexes of the models with the optimal inputs (column 1) and 
%                   the index of the optimal single model (column 2) in 
%                   "QSMs" for each tree 
% OptInputs     The optimized input parameters for each tree
% OptQSMs       The single best QSM for each tree, OptQSMs = QSMs(OptModel);
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 1.1.0, 08 Oct 2019:
% 1) Added the posibility to select the optimisation criteria or cost
%    function ("metric") that is minimised from 34 different options.
%    Previously only one option was used. The used metric is also included
%    in "OptInputs" output as one of the fields
% 2) Added OptQSM as one of the outputs

%% Collect data
m = max(size(QSMs)); % number of models
IndAll = (1:1:m)';
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

treedata = zeros(n,m); % Collect all tree attributes from all models
inputs = zeros(m,5); % collect the inputs from all models 
                    % ([PatchDiam1 PatchDiam2Min PatchDiam2Max lcyl FilRad])
dist = zeros(m,10); % collect the distances (metrics) from all models
TreeId = zeros(m,2); % collectd the tree and model indexes from all models
Keep = true(m,1); % Non-empty models
for i = 1:m
    if ~isempty(QSMs(i).cylinder)
        p = QSMs(i).rundata.inputs;
        inputs(i,:) = [p.PatchDiam1 p.PatchDiam2Min p.PatchDiam2Max p.lcyl p.FilRad];
        D = QSMs(i).pmdistance;
        % Collect distances: mean of all cylinders (cylinder-to-point distances),
        % mean of trunk, branch, 1st- and 2nd-order branch cylinders.
        % And the maximum of the previous.
        dist(i,:) =  [D.mean  D.TrunkMean  D.BranchMean  D.Branch1Mean  D.Branch2Mean ...
            D.max  D.TrunkMax  D.BranchMax  D.Branch1Max  D.Branch2Max];
        
        for j = 1:n
            treedata(j,i) = QSMs(i).treedata.(names{j});
        end
        TreeId(i,:) = [p.tree p.model];
    else
        Keep(i) = false;
    end
end
treedata = treedata(:,Keep);
inputs = inputs(Keep,:);
dist = dist(Keep,:);
TreeId = TreeId(Keep,:);
IndAll = IndAll(Keep);
TreeIds = unique(TreeId(:,1));
nt = length(TreeIds); % number of trees

% Determine the input parameter values
InputParComb = unique(inputs,'rows'); % Input parameter combinations
IV = cell(5,1);
N = zeros(5,1);
for i = 1:5
    I = unique(InputParComb(:,i));
    IV{i} = I;
    N(i) = length(I);
end

% Select the metric based on the input
if nargin > 1
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
    elseif strcmp(Metric,'tot_vol')
        met = 28;
    elseif strcmp(Metric,'trunk_vol')
        met = 29;
    elseif strcmp(Metric,'branch_vol')
        met = 30;
    elseif strcmp(Metric,'trunk+branch_vol')
        met = 31;
    elseif strcmp(Metric,'trunk_len')
        met = 32;
    elseif strcmp(Metric,'branch_len')
        met = 33;
    elseif strcmp(Metric,'branch_num')
        met = 34;
    else
        met = 1;
        Metric = 'all_mean_dis';
    end
else
    met = 1;
    Metric = 'all_mean_dis';
end

% The metric for selecting the optimal single model from the models with
% the optimal inputs is the mean point-model-distance.
best = 1;

%% Determine metric-value for each input 
% (average over number of models with the same inputs)
input = cell(nt,N(1)*N(2)*N(3)*N(4)*N(5));
distM = zeros(nt,N(1)*N(2)*N(3)*N(4)*N(5)); % average distances or volume stds
for t = 1:nt
    I = TreeId(:,1) == TreeIds(t);
    b = 0;
    for d = 1:N(1) % PatchDiam1
        J = abs(inputs(:,1)-IV{1}(d)) < 0.0001;
        for a = 1:N(2) % PatchDiam2Min
            K = abs(inputs(:,2)-IV{2}(a)) < 0.0001;
            for i = 1:N(3) % PatchDiam2Max
                L = abs(inputs(:,3)-IV{3}(i)) < 0.0001;
                for l = 1:N(4) % lcyl
                    M = abs(inputs(:,4)-IV{4}(l)) < 0.0001;
                    for f = 1:N(5) % FilRad
                        O = abs(inputs(:,5)-IV{5}(f)) < 0.0001;
                        
                        T = I&J&K&L&M&O;
                        b = b+1;
                        input{t,b} = [d a i l f];
                        % Compute the metric from all models with the same inputs
                        D = mean(dist(T,:),1); % Metric value is the mean from 
                                    % all the models with the same inputs
                        D(6:10) = 0.5*D(6:10); % Half the maximum values
                        if met == 1
                            D = D(1);
                        elseif met == 2
                            D = D(2);
                        elseif met == 3
                            D = D(3);
                        elseif met == 4
                            D = D(4);
                        elseif met == 5
                            D = D(5);
                        elseif met == 6
                            D = D(2)+D(3);
                        elseif met == 7
                            D = D(2)+D(4);
                        elseif met == 8
                            D = D(2)+D(4)+D(5);
                        elseif met == 9
                            D = D(4)+D(5);
                            
                        elseif met == 10
                            D = D(6);
                        elseif met == 11
                            D = D(7);
                        elseif met == 12
                            D = D(8);
                        elseif met == 13
                            D = D(9);
                        elseif met == 14
                            D = D(10);
                        elseif met == 15
                            D = D(7)+D(8);
                        elseif met == 16
                            D = D(7)+D(9);
                        elseif met == 17
                            D = D(7)+D(9)+D(10);
                        elseif met == 18
                            D = D(9)+D(10);
                            
                        elseif met == 19
                            D = D(1)+D(6);
                        elseif met == 20
                            D = D(2)+D(7);
                        elseif met == 21
                            D = D(3)+D(8);
                        elseif met == 22
                            D = D(4)+D(9);
                        elseif met == 23
                            D = D(5)+D(10);
                        elseif met == 24
                            D = D(2)+D(3)+D(7)+D(8);
                        elseif met == 25
                            D = D(2)+D(4)+D(7)+D(9);
                        elseif met == 26
                            D = D(2)+D(4)+D(5)+D(7)+D(9)+D(10);
                        elseif met == 27
                            D = D(4)+D(5)+D(9)+D(10);
                        
                        elseif met == 28
                            D = std(treedata(1,T));
                        elseif met == 29
                            D = std(treedata(2,T));
                        elseif met == 30
                            D = std(treedata(3,T));
                        elseif met == 31
                            D = std(treedata(2,T))+std(treedata(3,T));
                        elseif met == 32
                            D = std(treedata(5,T));
                        elseif met == 33
                            D = std(treedata(6,T));
                        elseif met == 34
                            D = std(treedata(7,T));
                        end
                        distM(t,b) = D;
                    end
                end
            end
        end
    end
end

%% Determine the optimal inputs and models
ninputs = prod(N);
OptIn = zeros(nt,5*min(ninputs,3)); % Optimal input values
OptDist = zeros(nt,min(ninputs,3)); % Smallest metric values
for i = 1:nt
    [d,J] = sort(distM(i,:));
    O = input{i,J(1)};
    OptIn(i,1:5) = [IV{1}(O(1)) IV{2}(O(2)) IV{3}(O(3)) IV{4}(O(4)) IV{5}(O(5))];
    OptDist(i,1) = d(1);
    if ninputs > 1
        O = input{i,J(2)};
        OptIn(i,6:10) = [IV{1}(O(1)) IV{2}(O(2)) IV{3}(O(3)) IV{4}(O(4)) IV{5}(O(5))];
        OptDist(i,2) = d(2);
        if ninputs > 2
            O = input{i,J(3)};
            OptIn(i,11:15) = [IV{1}(O(1)) IV{2}(O(2)) IV{3}(O(3)) IV{4}(O(4)) IV{5}(O(5))];
            OptDist(i,3) = d(3);
        end
    end
end

% Select the optimal models for each tree. That is, in the case of multiple 
% models with same inputs, select the one model with the optimal inputs 
% that has the minimum metric value among those models with the optimal inputs
OptModel = zeros(nt,min(ninputs,3)); % The indexes of the optimal single models in QSMs
OptModels = cell(nt,2); % The indexes of models in QSMs with the optimal inputs (col 1) 
                        % and the indexes of the optimal single models (col 2)
DataM = zeros(n,nt); % Mean of tree data for each tree computed from the optimal models
DataS = zeros(n,nt); % Standard deviation of tree data for each tree
DataM2 = DataM;     DataM3 = DataM;
DataS2 = DataS;     DataS3 = DataS;
for t = 1:nt
    I = TreeId(:,1) == TreeIds(t);
    J = abs(inputs(:,1)-OptIn(t,1)) < 0.0001;
    K = abs(inputs(:,2)-OptIn(t,2)) < 0.0001;
    L = abs(inputs(:,3)-OptIn(t,3)) < 0.0001;
    M = abs(inputs(:,4)-OptIn(t,4)) < 0.0001;
    O = abs(inputs(:,5)-OptIn(t,5)) < 0.0001;
    T = I&J&K&L&M&O;
    ind = IndAll(T);
    [~,T] = min(dist(ind,best));
    OptModel(t,1) = ind(T);
    OptModels{t,1} = ind;
    OptModels{t,2} = ind(T);
    DataM(:,t) = mean(treedata(:,ind),2);
    DataS(:,t) = std(treedata(:,ind),[],2);
    if ninputs > 1
        J = abs(inputs(:,1)-OptIn(t,6)) < 0.0001;
        K = abs(inputs(:,2)-OptIn(t,7)) < 0.0001;
        L = abs(inputs(:,3)-OptIn(t,8)) < 0.0001;
        M = abs(inputs(:,4)-OptIn(t,9)) < 0.0001;
        O = abs(inputs(:,5)-OptIn(t,10)) < 0.0001;
        T = I&J&K&L&M&O;
        ind = IndAll(T);
        [~,T] = min(dist(ind,best));
        OptModel(t,2) = ind(T);
        DataM2(:,t) = mean(treedata(:,ind),2);
        DataS2(:,t) = std(treedata(:,ind),[],2);
        if ninputs > 2
            J = abs(inputs(:,1)-OptIn(t,11)) < 0.0001;
            K = abs(inputs(:,2)-OptIn(t,12)) < 0.0001;
            L = abs(inputs(:,3)-OptIn(t,13)) < 0.0001;
            M = abs(inputs(:,4)-OptIn(t,14)) < 0.0001;
            O = abs(inputs(:,5)-OptIn(t,15)) < 0.0001;
            T = I&J&K&L&M&O;
            ind = IndAll(T);
            [~,T] = min(dist(ind,best));
            OptModel(t,3) = ind(T);
            DataM3(:,t) = mean(treedata(:,ind),2);
            DataS3(:,t) = std(treedata(:,ind),[],2);
        end
    end
end
if ninputs > 1
    OptModel2 = IndAll(OptModel(:,2));
    if ninputs > 2
        OptModel3 = IndAll(OptModel(:,3));
    end
end
OptModel = IndAll(OptModel(:,1));
for i = 1:nt
    OptInputs(i) = QSMs(OptModel(i,1)).rundata.inputs;
    if ninputs > 1
        OI2(i) = QSMs(OptModel2(i)).rundata.inputs;
        if ninputs > 2
            OI3(i) = QSMs(OptModel3(i)).rundata.inputs;
        end
    end
end
for i = 1:nt
    OptInputs(i).metric = Metric;
end
OptQSM = QSMs(OptModel);
DataCV = DataS./DataM*100; % Coefficient of variation
if ninputs > 1
    DataCV2 = DataS2./DataM2*100; % Coefficient of variation
    if ninputs > 2
        DataCV3 = DataS3./DataM3*100; % Coefficient of variation
    end
end

%% Display some data about optimal models 
% Decrease the number on non-zero decimals
for j = 1:nt
    DataM(:,j) = change_precision(DataM(:,j));
    DataS(:,j) = change_precision(DataS(:,j));
    DataCV(:,j) = change_precision(DataCV(:,j));
    if ninputs > 1
        DataM2(:,j) = change_precision(DataM2(:,j));
        DataS2(:,j) = change_precision(DataS2(:,j));
        DataCV2(:,j) = change_precision(DataCV2(:,j));
        if ninputs > 2
            DataM3(:,j) = change_precision(DataM3(:,j));
            DataS3(:,j) = change_precision(DataS3(:,j));
            DataCV3(:,j) = change_precision(DataCV3(:,j));
        end
    end
end

Names = names(1:n);
L = max(cellfun('length',Names))+1;
for i = 1:n
    name = Names{i};
    name(L) = ' ';
    Names{i} = name;
end

% Display optimal inputs, model and attributes for each tree
for t = 1:nt
    disp('-------------------------------')
    disp(['  Tree: ',num2str(OptInputs(t).tree),', ',OptInputs(t).name])
    if ninputs == 1
        disp(['    Metric: ',Metric])
        disp(['    Metric value:  ',num2str(1000*OptDist(t,1))])
        disp(['    Optimal inputs:  PatchDiam1 = ',num2str(OptInputs(t).PatchDiam1)])
        disp(['                  PatchDiam2Min = ',num2str(OptInputs(t).PatchDiam2Min)])
        disp(['                  PatchDiam2Max = ',num2str(OptInputs(t).PatchDiam2Max)])
        disp(['                           lcyl = ',num2str(OptInputs(t).lcyl)])
        disp(['                         FilRad = ',num2str(OptInputs(t).FilRad)])
        disp(['    Optimal model: ',num2str(OptModel(t))])
        sec = num2str(round(QSMs(OptModel(t)).rundata.time(end)));
        disp(['    Reconstruction time for the optimal model: ',sec,' seconds'])
        disp('    Attributes (mean, std, CV(%)):')
        for i = 1:n
            str = (['      ',Names{i},': ',num2str([DataM(i,t) DataS(i,t) DataCV(i,t)])]);
            disp(str)
        end
    elseif ninputs == 2
        disp('    The best two cases:')
        disp(['    Metric: ',Metric])
        disp(['    Metric values:  ',num2str(OptDist(t,:))])
        disp(['            inputs:  PatchDiam1 = ',...
            num2str([OptInputs(t).PatchDiam1 OI2(t).PatchDiam1])])
        disp(['                  PatchDiam2Min = ',...
            num2str([OptInputs(t).PatchDiam2Min OI2(t).PatchDiam2Min])])
        disp(['                  PatchDiam2Max = ',...
            num2str([OptInputs(t).PatchDiam2Max OI2(t).PatchDiam2Max])])
        disp(['                           lcyl = ',...
            num2str([OptInputs(t).lcyl OI2(t).lcyl])])
        disp(['                         FilRad = ',...
            num2str([OptInputs(t).FilRad OI2(t).FilRad])])
        disp(['    Optimal model: ',num2str(OptModel(t))])
        sec = num2str(round(QSMs(OptModel(t)).rundata.time(end)));
        disp(['    Reconstruction time for the optimal model: ',sec,' seconds'])
        disp('    Attributes (mean, std, CV(%), second best mean):')
        for i = 1:n
            str = (['      ',Names{i},':  ',num2str([DataM(i,t) ...
                DataS(i,t) DataCV(i,t) DataM2(i,t)])]);
            disp(str)
        end
    elseif ninputs > 2
        disp('    The best three cases:')
        disp(['    Metric: ',Metric])
        disp(['    Metric values:  ',num2str(OptDist(t,:))])
        disp(['            inputs:  PatchDiam1 = ',...
            num2str([OptInputs(t).PatchDiam1 OI2(t).PatchDiam1 OI3(t).PatchDiam1])])
        disp(['                  PatchDiam2Min = ',...
            num2str([OptInputs(t).PatchDiam2Min OI2(t).PatchDiam2Min OI3(t).PatchDiam2Min])])
        disp(['                  PatchDiam2Max = ',...
            num2str([OptInputs(t).PatchDiam2Max OI2(t).PatchDiam2Max OI3(t).PatchDiam2Max])])
        disp(['                           lcyl = ',...
            num2str([OptInputs(t).lcyl OI2(t).lcyl OI3(t).lcyl])])
        disp(['                         FilRad = ',...
            num2str([OptInputs(t).FilRad OI2(t).FilRad OI3(t).FilRad])])
        disp(['    Optimal model: ',num2str(OptModel(t))])
        sec = num2str(round(QSMs(OptModel(t)).rundata.time(end)));
        disp(['    Reconstruction time for the optimal model: ',sec,' seconds'])
        str = ['    Attributes (mean, std, CV(%),',...
            ' second best mean, third best mean, sensitivity):'];
        disp(str)
        for i = 1:n
            sensi = max(abs([DataM(i,t)-DataM2(i,t) DataM(i,t)-DataM3(i,t)])/DataM(i,t));
            sensi2 = 100*sensi;
            sensi = 100*sensi/DataCV(i,t);
            sensi2 = change_precision(sensi2);
            sensi = change_precision(sensi);
            str = (['      ',Names{i},':  ',num2str([DataM(i,t) DataS(i,t) DataCV(i,t)...
                DataM2(i,t) DataM3(i,t) sensi sensi2])]);
            disp(str)
        end
    end
    disp('------')
end

%% Generate TreeData sructure for optimal models
clear TreeData
TreeData = OptQSM.treedata;
for t = 1:nt
    for i = 1:n
        TreeData(t).(names{i}) = [DataM(i,t) DataS(i,t)];
    end
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


