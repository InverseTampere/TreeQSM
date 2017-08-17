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

function [TreeData,OptQSMs,OptQSM] = estimate_precision(QSMs,NewQSMs,OptModels,savename)

% Uses models with the same inputs to estimate the precision (standard
% deviation) of the results. Has two sets of models as its inputs: 
% 1) QSMs can contain models with many different input parameters for each tree 
% and OptModels contain the indexes of the models that are used here ("optimal 
% models"); 2) NewQSMs contains only models with the optimal inputs.
%
% Inputs:
% QSMs          Contain all the QSMs, possibly from multiple trees
% NewQSMs       Contains the additional models with optimal inputs, for all trees
% OptModels     Indexes of the optimal models for each tree in QSMs
% savename      Optional input, name string specifying the name of the saved
%                   file containing the outputs
% Outputs:
% TreeData      Similar structure array as the "treedata" in QSMs but now each
%                   single-number attribute contains the mean and std computed
%                   from the models with the optimal inputs. 
% OptQSMs       Contains all the models with the optimal inputs, for all trees
% OptQSM        The best model (minimum point-model distance) among the models 
%                   with the optimal inputs, for all trees

%% Reconstruct the outputs 
OptQSMs = QSMs(vertcat(OptModels{:,1})); % Optimal models from the optimization process
OptQSMs = [OptQSMs NewQSMs]; % Combine all the optimal QSMs

m = max(size(OptQSMs)); % number of models
IndAll = (1:1:m)';
% Find the first non-empty model
i = 1;
while isempty(OptQSMs(i).cylinder)
    i = i+1;
end
% Determine how many single-number attributes there are in treedata
names = fieldnames(OptQSMs(i).treedata);
n = 1;
while numel(OptQSMs(i).treedata.(names{n})) == 1
    n = n+1;
end
n = n-1;

treedata = zeros(n,m); % Collect all single-number tree attributes from all models
TreeId = zeros(m,1); % Collect tree and model indexes from all models
Dist = zeros(m,1); % Collect the distances
Keep = true(m,1); % Non-empty models
for i = 1:m
    if ~isempty(OptQSMs(i).cylinder)
        for j = 1:n
            treedata(j,i) = OptQSMs(i).treedata.(names{j});
        end
        TreeId(i) = OptQSMs(i).rundata.inputs.tree;
        Dist(i) = OptQSMs(i).pmdistance.mean;
    else
        Keep(i) = false;
    end
end
treedata = treedata(:,Keep);
TreeId = TreeId(Keep,:);
Dist = Dist(Keep);
IndAll = IndAll(Keep);
TreeIds = unique(TreeId);
nt = length(TreeIds); % number of trees

% Compute the means and standard deviations
OptModel = zeros(nt,1);
DataM = zeros(n,nt);
DataS = zeros(n,nt);
for t = 1:nt
    I = TreeId == TreeIds(t);
    ind = IndAll(I);
    dist = vertcat(Dist(ind));
    [~,J] = min(dist);
    OptModel(t) = ind(J);
    DataM(:,t) = mean(treedata(:,ind),2);
    DataS(:,t) = std(treedata(:,ind),[],2);
end
OptQSM = OptQSMs(OptModel);
DataCV = DataS./DataM*100;

%% Display some data about optimal models 
% Decrease the number of non-zero decimals
for j = 1:nt
    DataM(:,j) = change_precision(DataM(:,j));
    DataS(:,j) = change_precision(DataS(:,j));
    DataCV(:,j) = change_precision(DataCV(:,j));
end

% Display optimal inputs, model and attributes for each tree
for t = 1:nt
    disp(['  Tree: ',num2str(t)])
    disp('    Attributes (mean, std, CV(%)):')
    for i = 1:n
        str = (['     ',names{i},': ',num2str([DataM(i,t) DataS(i,t) DataCV(i,t)])]);
        disp(str)
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
if nargin == 4
    str = ['results/OptimalQSMs_',savename];
    save(str,'TreeData','OptQSMs','OptQSM')
end
