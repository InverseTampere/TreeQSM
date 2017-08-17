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

function QSMs = make_models(dataname,savename,Nmodels,inputs)

% Makes QSMs of given point clouds specified by the "dataname" and by the
% other inputs. The results are saved into file named "savename".
% Notice, the code does not save indivual QSM runs into their own .mat or
% .txt files but saves all models into one big .mat file.

% Inputs:
% dataname      String specifying the .mat-file containing the point
%                   clouds that are used for the QSM reconstruction.
% savename      String specifying the name of the file where the QSMs are saved
% Nmodels       (Optional) Number of models generated for each input 
%                   (cloud and input parameters). Default value is 5.
% inputs        (Optional) The input parameters structure. Can be defined
%                   below as part of this code. Can also be given as a
%                   structure array where each tree gets its own, possibly 
%                   uniquely, defined parameters (e.g. optimal parameters) 
%                   but each tree has to have same number of parameter values.
%
% Output:
% QSMs          Structure array containing all the QSMs generated

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
    clear inputs
    % The following parameters can be varied and should be optimised (can have multiple values):
    inputs.PatchDiam1 = [0.12 0.16]; % Patch size of the first uniform-size cover
    inputs.PatchDiam2Min = [0.02 0.03]; % Minimum patch size of the cover sets in the second cover
    inputs.PatchDiam2Max = [0.1]; % Maximum cover set size in the stem's base in the second cover
    inputs.lcyl = [3]; % Relative (length/radius) length of the cylinders,
    inputs.FilRad = [3]; % Relative radius for outlier point filtering,
    
    % The following parameters can be varied and but usually can be kept as
    % shown (i.e. as little bigger than PatchDiam parameters):
    inputs.BallRad1 = inputs.PatchDiam1+0.02; % Ball radius used for the first uniform-size cover generation
    inputs.BallRad2 = inputs.PatchDiam2Max+0.01; % Maximum ball radius used for the second cover generation
    
    % The following parameters can be usually kept fixed as shown:
    inputs.nmin1 = 3; % Minimum number of points in BallRad1-balls, generally good value is 3
    inputs.nmin2 = 1; % Minimum number of points in BallRad2-balls, generally good value is 1
    inputs.OnlyTree = 1; % If "1", then point cloud contains points only from the tree
    inputs.Tria = 0; % If "1", then triangulation produces
    inputs.Dist = 1; % If "1", then computes the point-model distances
    
    % Different cylinder radius correction options for modifying too large and
    % too small cylinders:
    % Traditional TreeQSM choices:
    inputs.MinCylRad = 0.0025; % Minimum cylinder radius, used particularly in the taper corrections
    inputs.ParentCor = 1; % Child branch cylinders radii are always smaller than the parent branche's cylinder radii
    inputs.TaperCor = 1; % Use partially linear (stem) and parabola (branches) taper corrections
    % Growth volume correction approach introduced by Jan Hackenberg, allometry: GrowthVol = a*Radius^b+c
    inputs.GrowthVolCor = 0; % Use growth volume correction
    inputs.GrowthVolFac = 2.5; % fac-parameter of the growth vol. approach, defines upper and lower boundary
end


%% Load data
load([dataname,'.mat']) % 'dataname.mat' contains the point clouds that are modelled
S = whos; % names of the point cloud
n = max(size(S));
I = false(n,1);
for i = 1:n
    if S(i).size(1) > 100 && S(i).size(2) == 3
        I(i) = true;
    end
end
S = S(I);

%% make the models
nt = size(S,1); % number of trees/point clouds
QSMs = struct('cylinder',{},'branch',{},'treedata',{},'rundata',{},'pmdistance',{},'triangulation',{});

% Generate Inputs struct that contains the input parameters for each tree
if max(size(inputs)) == 1
    clear Inputs
    for i = 1:nt
        Inputs(i) = inputs;
    end
else
    Inputs = inputs;
end

m = 1;
for t = 1:nt % trees
    P = eval(getfield(S(t),'name'));
    str = getfield(S(t),'name');
    j = 1;
    inputs = Inputs(t);
    inputs.name = str;
    inputs.tree = t;
    inputs.plot = 0;
    inputs.savetxt = 0;
    inputs.savemat = 0;
    inputs.disp = 1;
    while j <= Nmodels % generate N models per input
        k = 1;
        n0 = 1;
        inputs.model = j;
        while k <= 5 % try up to five times to generate non-empty models
            QSM = treeqsm(P,inputs);
            n = max(size(QSM));
            Empty = false(n,1);
            for b = 1:n
                if isempty(QSM(b).branch)
                    Empty(b) = true;
                end
            end
            if any(Empty)
                n = nnz(~Empty);
                k = k+1;
                if n > n0
                    qsm = QSM;
                    n0 = n;
                end
            else
                % Succesfull models generated
                QSMs(m:m+n-1) = QSM;
                m = m+n;
                k = 10;
            end
        end
        if k == 6
            disp('Incomplete run!!')
            QSMs(m:m+n0-1) = qsm;
            m = m+n0;
        end
        j = j+1;
        if j > 0 || floor(j/5) == ceil(j/5)
            stri = ['results/',savename];
            save(stri,'QSMs')
        end
    end
    stri = ['results/',savename];
    save(stri,'QSMs')
end
