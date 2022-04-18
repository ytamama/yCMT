function [areaval,time1,time2,pktime]=phasearea(sacfile,phaseortime,...
    noisewin,threshval)
% 
% Function to calculate the area under the peak, corresponding to the 
% first arrival of an inputted phase. 
% 
% INPUTS
% sacfile : String, containing the full path to a SAC file
% phaseortime : Two types of inputs: 
% 
%               1. Which type of arrival? Enter one of the following strings: 
%               'P', 'S', or '_CODA'
% 
%               2. Enter the time of arrival as determined by 
%                  SeismoTech Ltd., in UTC.
% 
% noisewin : Number of seconds before the arrival of a phase to consider
%            as noise
% threshval : We consider a peak as the arrival if it is at least the mean
%             amplitude of the noise by x times. Input this value x. 
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
% time2 : Time of the first zero-crossing or local minimum/maximum after 
%         the peak, in UTC.
% pktime : Time of the peak, corresponding to the first arriving phase.
% 
% Last Modified : January 21, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find time of peak of inputted phase
[~,pktime,~]=phasepk(sacfile,phaseortime,noisewin,threshval);

% Calculate area of peak
if isempty(pktime)
  areaval=[];
  time1=[];
  time2=[];
else
  [areaval,time1,~,time2,~]=pkarea(sacfile,pktime);
end


