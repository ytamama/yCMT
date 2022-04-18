function [maxdat,mindat]=sacbounds(sacfiles,atims,arrtime)
% 
% Function to retrieve the maximum and minimum data values in the inputted
% SAC files. 
% 
% INPUT
% sacfiles : Cell array containing the full paths to the SAC files we wish
%            to process. 
% atims : Do we use only a segment of the SAC file, X seconds before to Y 
%         seconds after a given arrival? Enter [X Y], where X is the number
%         of seconds before arrival, and Y is number of seconds after 
%         arrival. Default: []
% arrtime : If using only a segment wrt arrival time, specify that arrival
%           time in UTC. Default: []
% 
% OUTPUTS
% maxdat : Maximum data value across all inputted SAC files
% mindat : Minimum data value across all inputted SAC files
% 
% Last Modified: January 8, 2022 by Yuri Tamama
% 
% Uses readsac.m from csdms-contrib/slepian_oscar
% Uses defval.m from csdms-contrib/slepian_alpha
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default variables
defval('atims',[])
defval('arrtime',[])

% Initial values for minimum and maximum data values
maxdat=0;
mindat=0;

% Iterate through all SAC files
for s=1:length(sacfiles)
  nowsac=sacfiles{s};
  % Y-axis bounds
  [sacdata,~]=readsac(nowsac);
  % If using segment wrt arrival time
  if ~isempty(atims) && ~isempty(arrtime)
    % Start time of SAC file
    [sactime,tims]=getsactime(nowsac);
    difftime=seconds(arrtime-sactime);
    starttim=difftime-atims(1);
    fintim=difftime+atims(2);
    % Indices spanning segment
    timinds=1:length(tims);
    arrinds=timinds((tims>=starttim) & (tims<=fintim));
    sacdata=sacdata(arrinds);
  end
  if s==1
    maxdat=max(sacdata);
    mindat=min(sacdata);
  else 
    if max(sacdata)>=maxdat
      maxdat=max(sacdata);
    end
    if min(sacdata)<=mindat
      mindat=min(sacdata);
    end
  end
end





