function plot_distribution(QSM,fig,rela,cumu,dis,dis2,dis3,dis4)

% ---------------------------------------------------------------------
% PLOT_DISTRIBUTION     Plots the specified distribution(s) in the 
%                           "treedata" field of the QSM structure array.
%
% Version 1.1.0
% Latest update     3 May 2022
%
% Copyright (C) 2020-2022 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Inputs:
% QSM       The output of treeqsm function, may contain multiple models if
%               only one distribution. If multiple distributions are plotted, 
%               then only one model.
% fig       Figure number
% rela      If rela = 1, then plots relative values (%), otherwise plots 
%               absolute values
% cumu      If cumu = 1, then plot cumulative distribution
% dis       Distribution to be plotted, string name, e.g. 'VolCylDia'.
%               The name string is the one used in the "treedata"
% dis2      Optional, Second distribution to be plotted. Notice with more
%               than one distribution, only one model.
% dis3      Optional, Third distribution to be plotted
% dis4      Optional, Fourth distribution to be plotted
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 1.1.0, 3 May 2022:
% 1) Added new input "cum" for plottig the distributions as cumulative.
% 2) Added return if distributions are empty or all zero

% Generate strings for title, xlabel and ylabel:
if strcmp(dis(1:3),'Vol')
  str = 'volume';
  ylab = 'Volume (L)';
elseif strcmp(dis(1:3),'Are')
  str = 'area';
  ylab = 'Area (m^2)';
elseif strcmp(dis(1:3),'Len')
  str = 'length';
  ylab = 'Length (m)';
elseif strcmp(dis(1:3),'Num')
  str = 'number';
  ylab = 'Number';
end

if strcmp(dis(end-2:end),'Dia')
  str2 = 'diameter';
  xlab = 'diameter (cm)';
elseif strcmp(dis(end-2:end),'Hei')
  str2 = 'height';
  xlab = 'height (m)';
elseif strcmp(dis(end-2:end),'Ord')
  str2 = 'order';
  xlab = 'order';
elseif strcmp(dis(end-2:end),'Ang')
  str2 = 'angle';
  xlab = 'angle (deg)';
elseif strcmp(dis(end-2:end),'Azi')
  str2 = 'azimuth direction';
  xlab = 'azimuth direction (deg)';
elseif strcmp(dis(end-2:end),'Zen')
  str2 = 'zenith direction';
  xlab = 'zenith direction (deg)';
end

% Collect the distributions
if nargin == 5
  % Multiple QSMs, one and the same distribution
  m = max(size(QSM));
  D = QSM(1).treedata.(dis);
  n = size(D,2);
  for i = 2:m
    d = QSM(i).treedata.(dis);
    k = size(d,2);
    if k > n
      n = k;
      D(m,n) = 0;
      D(i,1:n) = d;
    elseif k < n
      D(i,1:k) = d;
    else
      D(i,:) = d;
    end
  end
  D = D(:,1:n);
else
  % One QSM, multiple distributions of the same type
  % (e.g. diameter distributions: 'NumCylDia', 'VolCylDia' and 'LenCylDia')
  m = nargin-4;
  D = QSM.treedata.(dis);
  n = size(D,2);
  if n == 0 || all(D == 0)
    return
  end
  for i = 2:m
    if i == 2
      D(m,n) = 0;
      D(i,:) = QSM.treedata.(dis2);
    elseif i == 3
      D(i,:) = QSM.treedata.(dis3);
    else
      D(i,:) = QSM.treedata.(dis4);
    end
  end
end

if rela
  % use relative value
  for i = 1:m
    D(i,:) = D(i,:)/sum(D(i,:))*100;
  end
  ylab = 'Relative value (%)';
end

if cumu
  % use cumulative distribution
  D = cumsum(D,2);
end

% Generate the bar plot
figure(fig)
if strcmp(dis(end-3:end),'hAzi') || strcmp(dis(end-3:end),'1Azi') || strcmp(dis(end-2:end),'Azi')
  bar(-170:10:180,D')
elseif strcmp(dis(end-2:end),'Zen') || strcmp(dis(end-2:end),'Ang')
  bar(10:10:10*n,D')
else
  bar(1:1:n,D')
end

% Generate the title of the plot
if strcmp(dis(end-2:end),'Ord') && ~strcmp(dis(1:3),'Num')
  tit = ['Branch ',str,' per branching order'];
elseif strcmp(dis(end-2:end),'Ord')
  tit = 'Number of branches per branching order';
elseif strcmp(dis(1:3),'Num')
  tit = ['Number of branches per ',str2,' class'];
elseif strcmp(dis(end-3),'h') || strcmp(dis(end-3),'1')
  tit = ['Branch ',str,' per ',str2,' class'];
else
  tit = ['Tree segment ',str,' per ',str2,' class'];
end
n = nargin;
if n > 5
  if ~strcmp(dis(1:3),dis2(1:3))
    if strcmp(dis(4),'C')
      tit = 'Tree segment distribution';
    else
      tit = 'Branch distribution';
    end
  elseif n > 6
    if ~strcmp(dis(1:3),dis3(1:3))
      if strcmp(dis(4),'C')
        tit = 'Tree segment distribution';
      else
        tit = 'Branch distribution';
      end
    elseif n > 7
      if ~strcmp(dis(1:3),dis4(1:3))
        if strcmp(dis(4),'C')
          tit = 'Tree segment distribution';
        else
          tit = 'Branch distribution';
        end
      end
    end
  end
end
title(tit)

% Generate the x-axis label
if strcmp(dis(end-5:end-3),'Cyl')
  xlab = ['Cylinder ',xlab];
else
  xlab = ['Branch ',xlab];
end
xlabel(xlab)

% Generate the y-axis label
ylabel(ylab);

% Tight axes and grid lines
axis tight
grid on

m = max(size(QSM));
% Add legends, if needed
if m > 1
  L = cell(m,1);
  for i = 1:m
    L{i} = ['model',num2str(i)];
  end
  legend(L,'location','best')
elseif nargin > 5
  m = nargin-4;
  L = cell(m,1);
  for i = 1:m
    if i == 1
      L{i} = dis(1:end-3);
    elseif i == 2
      L{i} = dis2(1:end-3);
    elseif i == 3
      L{i} = dis3(1:end-3);
    else
      L{i} = dis4(1:end-3);
    end
  end
  legend(L,'location','best')
end