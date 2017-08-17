function display_time(T1,T2,string,display)

% Display the two times given. "T1" is the time named with the "string" and
% "T2" is named "Total".

[tmin,tsec] = sec2min(T1);
[Tmin,Tsec] = sec2min(T2);
if tmin < 60 && Tmin < 60
    if tmin < 1 && Tmin < 1
        str = [string,' ',num2str(tsec),' sec.   Total: ',num2str(Tsec),' sec'];
    elseif tmin < 1
        str = [string,' ',num2str(tsec),' sec.   Total: ',num2str(Tmin),...
            ' min ',num2str(Tsec),' sec'];
    else
        str = [string,' ',num2str(tmin),' min ',num2str(tsec),...
            ' sec.   Total: ',num2str(Tmin),' min ',num2str(Tsec),' sec'];
    end
elseif tmin < 60
    Thour = floor(Tmin/60);
    Tmin = Tmin-Thour*60;
    str = [string,' ',num2str(tmin),' min ',num2str(tsec),...
    	' sec.   Total: ',num2str(Thour),' hours ',num2str(Tmin),' min'];
else
    thour = floor(tmin/60);
    tmin = tmin-thour*60; 
    Thour = floor(Tmin/60);
    Tmin = Tmin-Thour*60; 
    str = [string,' ',num2str(thour),' hours ',num2str(tmin),...
  	' min.   Total: ',num2str(Thour),' hours ',num2str(Tmin),' min'];
end
if display
    disp(str)
end