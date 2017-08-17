function [Tmin,Tsec] = sec2min(T)

% Transforms the given number of seconds into minutes and residual seconds

Tmin = floor(T/60);
Tsec = round((T-Tmin*60)*10)/10;