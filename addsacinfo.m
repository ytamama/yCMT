function addsacinfo(sacfiles,evtnum,stanum)
% 
% Function to add information pertaining to event location, station 
% location, and picked arrivals to inputted SAC files. 
% 
% INPUTS
% sacfiles : Cell array, containing the full paths to SAC files. Each SAC
%            file is inputted as a string. 
% evtnum : ID number, corresponding to the event recorded on the SAC files
% stanum : ID number, corresponding to the station at which these SAC
%          files were recorded. 
% 
% 
% Last Modified : January 3, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Retrieve event information
evtinfo=getevtinfo(evtnum);
evlo=evtinfo.LON_dec;
evla=evtinfo.LAT_dec;
evdp=evtinfo.DEP_KM;
mag=evtinfo.MAG;

% Retrive station information
stainfo=getstainfo(stanum);
stlo=stainfo.Longitude;
stla=stainfo.Latitude;
stel=-1*stainfo.Altitude; 

% Iterate through SAC files and add necessary information
for f=1:length(sacfiles)
  sacfile=sacfiles{f};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Add event info
  % ID number
  evnumcmd=sprintf('echo "r %s ; chnhdr KEVNM %d ; w %s ; q" | /usr/local/sac/bin/sac',...
    sacfile,evtnum,sacfile);
  [status,cmdout]=system(evnumcmd);
      
  % Coordinates
  evcmd=sprintf('echo "r %s ; chnhdr EVLO %f EVLA %f EVDP %f ; w %s ; q" | /usr/local/sac/bin/sac',...
     sacfile,evlo,evla,evdp,sacfile);
  [status,cmdout]=system(evcmd);
    
  % Mmagnitude
  magcmd=sprintf('echo "r %s ; chnhdr MAG %f ; w %s ; q" | /usr/local/sac/bin/sac',...
    sacfile,mag,sacfile);
  [status,cmdout]=system(magcmd);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Add station info
  stacmd=sprintf('echo "r %s ; chnhdr STLO %f STLA %f STEL %f ; w %s ; q" | /usr/local/sac/bin/sac',...
     sacfile,stlo,stla,stel,sacfile);
  [status,cmdout]=system(stacmd);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Obtain information on picked arrivals
  
  % Component
  [~,hdrinfo]=readsac(sacfile);
  comp=strtrim(hdrinfo.KCMPNM);
  
  % Arrivals
  sacvars={evtnum,stanum,comp};
  picktbl=getpickinfo(sacvars,{});
  [numpicks,~]=size(picktbl);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Add picks information to SAC files
  
  % Starting time of SAC file
  sactime=getsactime(sacfile);
  
  % Iterate through picks
  for p=1:numpicks
    pickrow=picktbl(p,:);
    picktype=strtrim(pickrow.ArrivalType);
    picktime=seconds(pickrow.ArrivalTime-sactime);
    pickqual=strtrim(pickrow.ArrivalQuality{1});
    if strcmp(pickqual,'?')
      pickqual=-12345;
    else
      try
        pickqual=str2double(pickqual);
      catch
        keyboard
      end
    end
    
    % P
    if strcmp(picktype,'P')
      pickcmd=sprintf('echo "r %s ; chnhdr T0 %f USER0 %d ; w %s ; q" | /usr/local/sac/bin/sac',...
        sacfile,picktime,pickqual,sactime);
        
    % S
    elseif strcmp(picktype,'S')
      pickcmd=sprintf('echo "r %s ; chnhdr T1 %f USER1 %d ; w %s ; q" | /usr/local/sac/bin/sac',...
        sacfile,picktime,pickqual,sactime);
        
    % _CODA
    elseif strcmp(picktype,'_CODA')
      pickcmd=sprintf('echo "r %s ; chnhdr T2 %f USER0 %2 ; w %s ; q" | /usr/local/sac/bin/sac',...
        sacfile,picktime,pickqual,sactime);
        
    % ? (Unknown phase)
    elseif strcmp(picktype,'?')
      pickcmd=sprintf('echo "r %s ; chnhdr T3 %f USER3 %d ; w %s ; q" | /usr/local/sac/bin/sac',...
        sacfile,picktime,pickqual,sactime);
    end
    [status,cmdout]=system(pickcmd);
  end  
end


