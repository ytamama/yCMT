function psdist=calcpsdist(ptime,stime,pvel,degorkm)
% 
% Program to calculate the distance between station and earthquake, 
% assuming that the medium through which the seismic waves travel is 
% a Poisson solid with Poisson's ratio of 0.25. This code follows
% Equation 6.1 of Modern Global Seismology by Lay and Wallace. 
% 
% INPUTS
% parrtime : Arrival time of the P-wave, in seconds relative to the event
%            origin OR the start of the seismogram, or as a datetime. 
% sarrtime : Arrival time of the S-wave, in a format consistent with that 
%            of the P-wave
% pvel : Velocity of P-waves, in km/s
% degorkm : Do we calculate distances using degrees or in Cartesian
%           coordinates? 
%           0 : Degrees
%           1 : Cartesian (in km)
% 
% OUTPUT
% psdist : P-to-S distance, in degrees or in km
% 
% Last Modified: January 7, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if arrival times are consistent
if ~strcmp(class(ptime),class(stime))
  fprintf('Both P and S arrival times must be in a consistent format \n')
  psdist=[];
  return
end

% Difference in arrival times
if isa(ptime,'datetime')
  difftime=seconds(stime-ptime);
  
elseif isa(ptime,'double')
  difftime=stime-ptime;
end

% Equation 6.1 of Lay & Wallace to calculate distance
psdist=(difftime*pvel)/(sqrt(3)-1);
if degorkm==0
  psdist=km2deg(psdist);
end
