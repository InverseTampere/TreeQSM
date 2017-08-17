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

function [TreeData,OptModels,OptInputs] = select_optimum(QSMs,savename)

% Selects optimum models based on the minimum mean cylinder distance
%
% Works for single or multiple tree cases where the input QSMs contains
% multiple models for the same tree with different inputs and multiple runs
% with the same inputs
%
% Inputs:
% QSMs          Contain all the QSMs, possibly from multiple trees
% savename      Optional input, name string specifying the name of the saved file
%               containing the outputs
% Outputs:
% TreeData      Similar structure array as the "treedata" in QSMs but now each
%               attribute contains the mean and std computed from the
%               models with the optimal inputs. 
% OptModels     Indexes of the models with the optimal inputs (column 1) and the 
%                   index of the optimal single model (column 2) in "QSMs" for each tree 
% OptInputs     The optimized input parameters for each tree


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
dist = zeros(m,5); % collect the distances (metrics) from all models
TreeId = zeros(m,2); % collectd the tree and model indexes from all models
Keep = true(m,1); % Non-empty models
for i = 1:m
    if ~isempty(QSMs(i).cylinder)
        p = QSMs(i).rundata.inputs;
        inputs(i,:) = [p.PatchDiam1 p.PatchDiam2Min p.PatchDiam2Max p.lcyl p.FilRad];
        D = QSMs(i).pmdistance;
        % Collect many distances: mean of all cylinders (cylinder-to-point distances),
        % mean of trunk cylinders, mean of all branch cylinders, mean of
        % 1st- and 2nd-order branch cylinders.
        % !!! Now uses only the first one, mean of all cylinder distances,
        % for the optimization
        dist(i,:) =  [D.mean D.TrunkMean D.BranchMean D.Branch1Mean D.Branch2Mean];
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

%% Determine distances for each input 
% (average over number of models with the same inputs)
input = cell(nt,N(1)*N(2)*N(3)*N(4)*N(5));
distM = zeros(nt,N(1)*N(2)*N(3)*N(4)*N(5));
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
                        % Compute the distance as the mean of all cylinder
                        % distances
                        D = mean(dist(T,1));
                        distM(t,b) = D;
                        % Compute the distance as the mean of all cylinder
                        % distances plus proportional standard deviations of
                        % total, stem and branch volumes
                        %D1 = mean(dist(T,1));
                        %D2 = std(treedata(1,T),[],2)/mean(treedata(1,T),2)/30;
                        %D3 = std(treedata(2,T),[],2)/mean(treedata(2,T),2)/30;
                        %D4 = std(treedata(3,T),[],2)/mean(treedata(3,T),2)/30;
                        %disp([D1+D2+D3+D4 D1 D2 D3 D4])
                        %distM(t,ind) = D1+D2+D3+D4;
                    end
                end
            end
        end
    end
end

%% Determine the optimal inputs and models
OptIn = zeros(nt,5); % Optimal input values
for i = 1:nt
    [~,J] = min(distM(i,:));
    O = input{i,J};
    OptIn(i,:) = [IV{1}(O(1)) IV{2}(O(2)) IV{3}(O(3)) IV{4}(O(4)) IV{5}(O(5))];
end

% Select the optimal models for each tree. That is, in the case of multiple 
% models with same inputs, select the one model with the optimal inputs 
% that has the minimum distance among those models with the optimal inputs
OptModel = zeros(nt,1); % The indexes of the optimal single models in QSMs
OptModels = cell(nt,2); % The indexes of models in QSMs with the optimal inputs (col 1) 
%                           and the indexes of the optimal single models (col 2)
DataM = zeros(n,nt); % Mean of tree data for each tree computed from the optimal models
DataS = zeros(n,nt); % Standard deviation of tree data for each tree
for t = 1:nt
    I = TreeId(:,1) == t;
    J = abs(inputs(:,1)-OptIn(t,1)) < 0.0001;
    K = abs(inputs(:,2)-OptIn(t,2)) < 0.0001;
    L = abs(inputs(:,3)-OptIn(t,3)) < 0.0001;
    M = abs(inputs(:,4)-OptIn(t,4)) < 0.0001;
    O = abs(inputs(:,5)-OptIn(t,5)) < 0.0001;
    T = I&J&K&L&M&O;
    ind = IndAll(T);
    [~,T] = min(dist(ind,1));
    OptModel(t) = ind(T);
    OptModels{t,1} = ind;
    OptModels{t,2} = ind(T);
    DataM(:,t) = mean(treedata(:,ind),2);
    DataS(:,t) = std(treedata(:,ind),[],2);
end
OptModel = IndAll(OptModel);
clear OptInputs
for i = 1:nt
    OptInputs(i) = QSMs(OptModel(i)).rundata.inputs;
end
OptQSM = QSMs(OptModel);
DataCV = DataS./DataM*100; % Coefficient of variation

%% Display some data about optimal models 
% Decrease the number on non-zero decimals
for j = 1:nt
    DataM(:,j) = change_precision(DataM(:,j));
    DataS(:,j) = change_precision(DataS(:,j));
    DataCV(:,j) = change_precision(DataCV(:,j));
end

% Display optimal inputs, model and attributes for each tree
for t = 1:nt
    disp(['  Tree: ',num2str(t)])
    disp(['    Optimal inputs:  PatchDiam1 = ',num2str(OptInputs(t).PatchDiam1)])
    disp(['                  PatchDiam2Min = ',num2str(OptInputs(t).PatchDiam2Min)])
    disp(['                  PatchDiam2Max = ',num2str(OptInputs(t).PatchDiam2Max)])
    disp(['                           lcyl = ',num2str(OptInputs(t).lcyl)])
    disp(['                         FilRad = ',num2str(OptInputs(t).FilRad)])
    disp(['    Optimal model: ',num2str(OptModel(t))])
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
if nargin == 2
    str = ['results/OptimalQSMs_',savename];
    save(str,'TreeData','OptModels','OptInputs','OptQSM')
end
