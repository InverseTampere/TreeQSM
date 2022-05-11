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

function QSMs = make_models_parallel(dataname,savename,Nmodels,inputs)

% ---------------------------------------------------------------------
% MAKE_MODELS.M       Makes QSMs of given point clouds.
%
% Version 1.1.2
% Latest update     9 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Makes QSMs of given point clouds specified by the "dataname" and by the
% other inputs. The results are saved into file named "savename".
% Notice, the code does not save indivual QSM runs into their own .mat or
% .txt files but saves all models into one big .mat file. Same as
% MAKE_MODELS but uses parfor command (requires Parallel Computing Toolbox)
% which allows the utilization of multiple processors/cores to compute in
% parallel number of QSMs with the same inputs.
%
% Inputs:
% dataname    String specifying the .mat-file containing the point
%               clouds that are used for the QSM reconstruction.
% savename    String, the name of the file where the QSMs are saved
% Nmodels     (Optional) Number of models generated for each input
%               (cloud and input parameters). Default value is 5.
% inputs      (Optional) The input parameters structure. Can be defined
%               below as part of this code. Can also be given as a
%               structure array where each tree gets its own, possibly
%               uniquely, defined parameters (e.g. optimal parameters)
%               but each tree has to have same number of parameter values.
%
% Output:
% QSMs        Structure array containing all the QSMs generated
% ---------------------------------------------------------------------

% Changes from version 1.1.1 to 1.1.2, 18 Aug 2020:
% 1) Removed the inputs "lcyl" and "FilRad" from the inputs and the
%    calculations of number of input parameters

% Changes from version 1.1.0 to 1.1.1, 13 Jan 2020:
% 1) Changed "m = m+n;" to "m = m+n(j);" at the end of the function.

% Changes from version 1.0.0 to 1.1.0, 03 Oct 2019:
% 1) Added try-catch structure where "treeqsm" is called, so that if there
%    is an error during the reconstruction process of one tree, then the
%    larger process of making multiple QSMs from multiple tree is not
%    stopped.
% 2) Changed the way the data is loaded. Previously all the data was
%    loaded into workspace, now only one point cloud is in the workspace.
% 3) Corrected a bug where incomplete QSM was saved as complete QSM
% 4) Changed where the input-structure for each tree reconstructed
% 5) Changed the coding to separate more the results of the different
%    parallel processes (less warnings and errors)

if nargin < 2
  disp('Not enough inputs, no models generated!')
  QSMs =  struct([]);
  return
end

if nargin == 2
  Nmodels = 5; % Number of models per inputs, usually about 5 models is enough
end

%% Define the parameter values
if nargin == 3 || nargin == 2
  % The following parameters can be varied and should be optimised
  % (each can have multiple values):
  % Patch size of the first uniform-size cover:
  inputs.PatchDiam1 = [0.08 0.15];
  % Minimum patch size of the cover sets in the second cover:
  inputs.PatchDiam2Min = [0.015 0.025];
  % Maximum cover set size in the stem's base in the second cover:
  inputs.PatchDiam2Max = [0.06 0.08];

  % The following parameters can be varied and but usually can be kept as
  % shown (i.e. as little bigger than PatchDiam parameters):
  % Ball radius used for the first uniform-size cover generation:
  inputs.BallRad1 = inputs.PatchDiam1+0.02;
  % Maximum ball radius used for the second cover generation:
  inputs.BallRad2 = inputs.PatchDiam2Max+0.01;

  % The following parameters can be usually kept fixed as shown:
  inputs.nmin1 = 3; % Minimum number of points in BallRad1-balls, good value is 3
  inputs.nmin2 = 1; % Minimum number of points in BallRad2-balls, good value is 1
  inputs.OnlyTree = 1; % If "1", then point cloud contains points only from the tree
  inputs.Tria = 0; % If "1", then triangulation produces
  inputs.Dist = 1; % If "1", then computes the point-model distances

  % Different cylinder radius correction options for modifying too large and
  % too small cylinders:
  % Traditional TreeQSM choices:
  % Minimum cylinder radius, used particularly in the taper corrections:
  inputs.MinCylRad = 0.0025;
  % Child branch cylinders radii are always smaller than the parent
  % branche's cylinder radii:
  inputs.ParentCor = 1;
  % Use partially linear (stem) and parabola (branches) taper corrections:
  inputs.TaperCor = 1;
  % Growth volume correction approach introduced by Jan Hackenberg,
  % allometry: GrowthVol = a*Radius^b+c
  % Use growth volume correction:
  inputs.GrowthVolCor = 0;
  % fac-parameter of the growth vol. approach, defines upper and lower
  % boundary:
  inputs.GrowthVolFac = 2.5;

  inputs.name = 'test';
  inputs.tree = 0;
  inputs.plot = 0;
  inputs.savetxt = 0;
  inputs.savemat = 0;
  inputs.disp = 0;
end

% Compute the number of input parameter combinations
in = inputs(1);
ninputs = prod([length(in.PatchDiam1) length(in.PatchDiam2Min)...
  length(in.PatchDiam2Max)]);


%% Load data
matobj = matfile([dataname,'.mat']);
names = fieldnames(matobj);
i = 1;
n = max(size(names));
while i <= n && ~strcmp(names{i,:},'Properties')
  i = i+1;
end
I = (1:1:n);
I = setdiff(I,i);
names = names(I,1);
names = sort(names);
nt = max(size(names)); % number of trees/point clouds

%% make the models
QSMs = struct('cylinder',{},'branch',{},'treedata',{},'rundata',{},...
  'pmdistance',{},'triangulation',{});

% Generate Inputs struct that contains the input parameters for each tree
if max(size(inputs)) == 1
  for i = 1:nt
    Inputs(i) = inputs;
    Inputs(i).name = names{i};
    Inputs(i).tree = i;
    Inputs(i).plot = 0;
    Inputs(i).savetxt = 0;
    Inputs(i).savemat = 0;
    Inputs(i).disp = 0;
  end
else
  Inputs = inputs;
end

m = 1;
for t = 1:nt % trees
  disp(['Modelling tree ',num2str(t),'/',num2str(nt),' (',Inputs(t).name,'):'])
  P = matobj.(Inputs(t).name);
  qsms = cell(Nmodels,1); % save here the accepted models
  qsm = cell(Nmodels,1); % cell-structure to keep different models separate
  n = ones(Nmodels,1);
  n0 = zeros(Nmodels,1);
  k = ones(Nmodels,1);
  parfor j = 1:Nmodels % generate N models per input
    inputs = Inputs(t);
    inputs.model = j;
    while k(j) <= 5 % try up to five times to generate non-empty models
      try
        qsm{j} = treeqsm(P,inputs);
      catch
        qsm{j} = struct('cylinder',{},'branch',{},'treedata',{},...
          'rundata',{},'pmdistance',{},'triangulation',{});
        qsm{j}(ninputs).treedata = 0;
      end
      n(j) = max(size(qsm{j}));
      Empty = false(n(j),1);
      for b = 1:n(j)
        if isempty(qsm{j}(b).branch)
          Empty(b) = true;
        end
      end
      if n(j) < ninputs || any(Empty)
        n(j) = nnz(~Empty);
        k(j) = k(j)+1;
        if n(j) > n0(j)
          qsms{j} = qsm{j}(~Empty);
          n0(j) = n(j);
        end
      else
        % Successful models generated
        qsms{j} = qsm{j};
        k(j) = 10;
      end
    end
    if k(j) == 6
      disp('Incomplete run!!')
    end
  end
  % Save the models
  for j = 1:Nmodels
    QSM = qsms{j};
    a = max(size(QSM));
    QSMs(m:m+a-1) = QSM;
    m = m+n(j);
  end
  str = ['results/',savename];
  save(str,'QSMs')
end
