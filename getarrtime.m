function arrtime=getarrtime(sacvars,phase)
% 
% Function to retrieve the arrival time of the inputted phase within the
% inputted SAC file. This arrival time is chosen by SeismoTech Ltd. 
% 
% INPUTS
% sacvars : One of two inputs : 
% 
%           1 : A string, containing a SAC file including its full path
%           2 : A three-element cell array, with the following entries : 
%               
%               evtnum : ID number of the earthquake
%               stanum : ID number of the station recording the event
%                 comp : Component, either 'Z', 'N', 'E', 'R', or 'T'
% 
% phase : String, containing one of three phases: 
%         'P', 'S', or '_CODA'
% 
% OUTPUT
% arrtime : Arrival time of phase in UTC
% 
% Last Modified : January 17, 2022 by Yuri Tamama
% 
% Uses readsac.m in csdms-contrib/slepian_oscar
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF a SAC file
if ischar(sacvars)
    
  % Read SAC header  
  [~,hdrinfo]=readsac(sacvars);
  sactime=getsactime(sacvars);
 
  % Try to find the arrival time in the SAC file itself
  if strcmp(phase,'P')
    ptim=hdrinfo.T0;
    if ptim==-12345
      arrtime=[];
    else
      arrtime=sactime+seconds(ptim);
    end
  elseif strcmp(phase,'S')
    stim=hdrinfo.T1;
    if stim==-12345
      arrtime=[];
    else
      arrtime=sactime+seconds(stim);
    end
  elseif strcmp(phase,'_CODA')
    ctim=hdrinfo.T2;
    if ctim==-12345
      arrtime=[];
    else
      arrtime=sactime+seconds(ctim);
    end
  else
    fprintf('Invalid phase \n');
    arrtime=[];
    return
  end

  % If arrival time is NOT in SAC file, consult the arrival catalog! 
  if isempty(arrtime)
    evtnum=str2double(hdrinfo.KEVNM);
    stanum=str2double(hdrinfo.KSTNM);
    if isnan(evtnum) || isnan(stanum)
      keyboard
    end
    comp=strtrim(hdrinfo.KCMPNM);
    sacvars2={evtnum,stanum,comp};
    arrtypesin={phase};
    pickrow=getpickinfo(sacvars2,arrtypesin);
  
    % If no arrival, return nothing
    if isempty(pickrow)
      arrtime=[];
      return
    
    else
      arrtime=pickrow.ArrivalTime;
    end
  end

  
% If a cell array input
elseif iscell(sacvars)
    
  % Retrieve information about this phase from catalog
  pickrow=getpickinfo(sacvars,{phase});
    
  % If no arrival, return nothing
  if isempty(pickrow)
    arrtime=[];
    return
    
  else
    arrtime=pickrow.ArrivalTime;
  end
end

