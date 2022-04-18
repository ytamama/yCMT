function [areaval,time1,val1,time2,val2]=pkarea(sacfile,pktime)
% 
% Function to calculate the area under the peak, corresponding to the 
% first arrival of a given phase
% 
% INPUTS
% sacfile : String, containing the full path to a SAC file
% pktime : Time corresponding to the maximum amplitude of the peak of a
%          given phase
% 
% OUTPUTS
% areaval : Area under the peak of the first arriving phase. If the 
%           arrival is upgoing, areaval is positive. If the arrival is 
%           downgoing, areaval is negative. 
% 
%           UNITS : 
%           Raw : Counts * s
%           Displacement : nm * s
%           Velocity : nm/s * s
%           Acceleration : nm/s^2 * s
% 
% time1 : Time of the first zero-crossing or local minimum/maximum before
%         the peak, in UTC. 
% val1 : Value of the seismogram determined at time1, either a zero
%        or a local minimum/maximum
% time2 : Time of the first zero-crossing or local minimum/maximum after 
%         the peak, in UTC.
% val2 : Value of the seismogram determined at time2, either a zero
%        or a local minimum/maximum
% 
% Last Modified : January 27, 2022 by Yuri Tamama
% 
% Uses readsac.m in csdms-contrib/slepian_oscar
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open SAC file
[sacdata,~,~,~,tims]=readsac(sacfile);

% Time at which SAC file starts
sactime=getsactime(sacfile);

% Index of peak
timinds=1:length(tims);
peakind=timinds(abs(seconds(pktime-(sactime+seconds(tims))))==min(abs(seconds(pktime-(sactime+seconds(tims))))));
% Sign of peak
peakval=sacdata(peakind);
peaksign=peakval/abs(peakval);

% Find closest zero-crossings 
% Before peak
changesign=0;
% Iterate backwards from peak to see where sign change occurs
nowind=peakind;
while changesign==0
  nowind=nowind-1;
  nowval=sacdata(nowind);
  nowsign=nowval/abs(nowval);
  
  % Sign change occurs
  if nowsign*peaksign==-1
    changesign=1;
    % Interpolate for zero-crossing
    intvals=[nowval sacdata(nowind+1)];
    inttims=[tims(nowind) tims(nowind+1)];
    ztim1=interp1(intvals,inttims,0);
    ztime1=sactime+seconds(ztim1);
    zval1=0;
    int1z=1;  % Interpolation to zero-crossing flag
   
  % Else if zero-crossing is found
  elseif nowsign*peaksign==0
    changesign=1;
    ztim1=tims(nowind);
    ztime1=sactime+seconds(ztim1);
    zval1=0;
    int1z=0;
  end
end


% Find first local minimum or maximum (opposite to the sign of 
% the peak) before peak 

% Find local maxima
if peaksign<0
  [pks,locs]=findpeaks(sacdata);
  
% Find local minima
elseif peaksign>0
  [pks,locs]=findpeaks(-1*sacdata); 
end
locs1=locs(locs<peakind);
pks1=pks(locs<peakind);
% Find local minimum/maximum closest to peak
loc1=locs1(abs(locs1-peakind)==min(abs(locs1-peakind)));
mval1=pks1(locs1==loc1);
% Find time of local min/max
mtim1=tims(loc1);
mtime1=sactime+seconds(mtim1);


% Which is closer to the peak? 
% Local minimum/maximum? 
if mtime1>ztime1
  time1=mtime1;
  val1=mval1;
  int1z=0;
  
% Or zero-crossing?
else
  time1=ztime1;
  val1=zval1;
end
tim1=seconds(time1-sactime);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After peak
changesign=0;
% Iterate forwards from peak to see where sign change occurs
nowind=peakind;
while changesign==0
  nowind=nowind+1;
  nowval=sacdata(nowind);
  nowsign=nowval/abs(nowval);
  
  % Sign change occurs
  if nowsign*peaksign==-1
    changesign=1;
    % Interpolate for zero-crossing
    intvals=[sacdata(nowind-1) nowval];
    inttims=[tims(nowind-1) tims(nowind)];
    ztim2=interp1(intvals,inttims,0);
    ztime2=sactime+seconds(ztim2);
    zval2=0;
    int2z=1;  % Interpolation to zero-crossing flag
    
  % Else if zero-crossing is found
  elseif nowsign*peaksign==0
    changesign=1;
    ztim2=tims(nowind);
    ztime2=sactime+seconds(ztim2);
    zval2=0;
    int2z=0;  
  end
end


% Find first local minimum or maximum (opposite to the sign of 
% the peak) after peak 

% Find local maxima
if peaksign<0
  [pks,locs]=findpeaks(sacdata);
  
% Find local minima
elseif peaksign>0
  [pks,locs]=findpeaks(-1*sacdata); 
end
locs2=locs(locs>peakind);
pks2=pks(locs>peakind);
% Find local minimum/maximum closest to peak
loc2=locs2(abs(locs2-peakind)==min(abs(locs2-peakind)));
mval2=pks2(locs2==loc2);
% Find time of local min/max
mtim2=tims(loc2);
mtime2=sactime+seconds(mtim2);


% Which is closer to the peak? 
% Local minimum/maximum? 
if mtime2<ztime2
  time2=mtime2;
  val2=mval2;
  int2z=0;
  
% Or zero-crossing?
else
  time2=ztime2;
  val2=zval2;
end
tim2=seconds(time2-sactime);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate area under the curve

% Obtain all points between zero-crossings, including peak
peakinds=timinds(tims>=tim1 & tims<=tim2);
peaktims=tims(peakinds);
peakvals=sacdata(peakinds).';
% If using an interpolated zero-crossing
if int1z==1
  peaktims=[tim1 peaktims];
  peakvals=[0 peakvals];
end
%
if int2z==1
  peaktims=[peaktims tim2];
  peakvals=[peakvals 0];
end

% Calculate area under the curve : trapezoidal approximation
try
  areaval=trapz(peaktims,peakvals);
catch
  keyboard
end


