function [figname,seisfig,allax,snrtex]=ploteq_seis(sacfiles,varargin)
%
% Function to plot a seismogram, recording an earthquake at a particular
% station, at one or more components. 
% 
% INPUTS
% sacfiles : Cell array, containing the SAC files we wish to plot. All SAC
%            files should contain the following entries: 
% 
%            KEVNM : ID Number of the earthquake recorded by the SAC file,
%            EVLA : Latitude of earthquake origin
%            EVLO : Longitude of earthquake origin
%            EVDP : Depth of earthquake origin, in km
%            KSTNM : ID Number of the station recording the SAC file
%            STLA : Latitude of station
%            STLO : Longitude of station
%            STEL : Elevation of station, in km
%            KCMPNM : Component recording the SAC file 
%            T0 : Time of arrival of the picked P-wave, in seconds 
%                 relative to the start of the SAC file
%            USER0 : Quality of P-wave arrival, chosen by SeismoTech Ltd. 
%            T1 : Time of arrival of the picked S-wave, in seconds 
%                 relative to the start of the SAC file
%            USER1 : Quality of S-wave arrival, chosen by SeismoTech Ltd. 
%            T2 : Time of arrival of the picked coda, in seconds 
%                 relative to the start of the SAC file
%            USER2 : Quality of coda arrival, chosen by SeismoTech Ltd. 
% 
% 
% Input the following as varargin: 
% 
% 'procInfo' : Cell array, containing information relevant to processing the 
%              inputted SAC files. 
% 
%              Format: {a; b; c; d; e} where
%  
%              a : 0 or 1, indicating if the data were already deconvolved 
%                  / filtered / etc. using the following information. 
% 
%                  0 : Not yet, please process
%                  1 : Already done
% 
%              b : Corner frequencies in Hz, through which we deconvolve our 
%                  data. Define a four-element vector with frequencies f1, 
%                  f2, f3, f4 in Hz, where f2 is at least twice f1, and f4 
%                  at least twice f3. Deconvolve the seismograms if we enter 
%                  this vector. If we don't wish to do so, enter an empty 
%                  vector. 
% 
%              c : Full path to SACPZ file, with which we deconvolve. Leave 
%                  empty if not deconvolving. 
% 
%              d : Do we deconvolve our SAC file into displacement, 
%                  velocity, or acceleration?
%                  0 : Displacement 
%                  1 : Velocity 
%                  2 : Acceleration 
%                  [] : Not deconvolving
% 
%              e : Frequencies through which we filter our data. Define a 2 
%                  element vector, where the first one is our lower bound 
%                  and the second one is our upper bound. Both in Hz. 
% 
%              Default: {}
% 
% 'markArrs' : Do we wish to mark arrivals? Enter a 3x1 cell:
% 
%              The first element is a cell array, containing strings
%              indicating which phases to mark: 
% 
%              'P' : P-arrivals
%              'S' : S-arrivals
%              {'P';'S'} : P and S arrivals
% 
%              The second element is if these arrivals are contained in 
%              the SAC files as variables, or if we need to consult a
%              catalog: 
% 
%              0 : In SAC files
%              1 : Need to consult a catalog. 
% 
%              The third element specifies whether we wish to limit the 
%              time axis of our seismograms, based on the time of the 
%              first-inputted arrival in the first element. Enter a 2x1 
%              vector, specifying how many seconds before to how many 
%              seconds after this arrival at which we plot our seismograms. 
% 
%              Default: {}
% 
% 'arrivalArea' : Do we want to plot the area under the peak corresponding
%                 to each marked arrival? 
% 
%                 0 : No [Default]
%                 1 : Yes
% 
% 'SNRmethod' : Calculate and plot the signal-to-noise ratio of each 
%               arrival? If so, then for how many seconds before or after? 
%               In what way do we calculate our SNR? Enter a 3-element 
%               vector, where the first element is the definition of SNR 
%               used: 
% 
%               1 : Amplitude: max(abs(signal)) / max(abs(noise)) 
%               2 : Variance: variance(signal) / variance(noise)
%               3 : Root mean square: RMS(signal) / RMS(noise)
% 
%           [4 x] : Peak-to-noise ratio: 
%                   Absolute amplitude of the first or second peak following 
%                   the arrival time, that is at least x times the mean peak 
%                   amplitude of the preceding noise. 
% 
%               For each definition of SNR, the second element (or third, 
%               if the first element is [4 x]) is: 
% 
%               1, 2, 3 : # of seconds before and after, which we consider 
%                         as signal and noise
%               4 : # of seconds before arrival, to consider as the noise
% 
%               Default: {}
%
% 'saveDirectory' : Full path to the directory at which we wish to save 
%                   our map. Enter an empty string if we don't wish to 
%                   save it. 
% 
%                   Default: ''
% 
% 'saveExtension' : What type of file do we save our seismograms as?
%                  'png' : .png 
%                  'pdf' : .pdf 
%                     '' : Default
% 
% 'seeFigure' : Figure visibility on or off? 
%               0 : Off [Default]
%               1 : On
% 
% 'makeFigure' : Make a new figure just for this, or plot on an existing
%                figure? 
%                0 : Plot on an existing figure 
%                1 : Make a new figure [Default]
% 
% 'cartCoords' : Do we assign Cartesian coordinates to the event and 
%                station corresponding to the inputted SAC files? If so, 
%                enter the following : 
% 
%                [[evx, evy, evz]; [stx, sty, stz]]
% 
%                evx : X-coordinate of event
%                evy : Y-coordinate of event 
%                evz : Z-coordinate of event (depth is positive)
%                stx : X-coordinate of station
%                sty : Y-coordinate of station
%                stz : Z-coordinate of station (depth is positive)
% 
% OUTPUT
% figname : Name of our seismogram and the full path to it. 
% seisfig : Current figure handle, containing our seismogram.
% allax : All axes objects
% snrtex : Text handle, to text indicating values of SNR of arrivals, 
%          if applicable. 
% 
% Uses readsac.m, in csdms-contrib/slepian_oscar
% Uses defval.m, in csdms-contrib/slepian_alpha
% Uses figdisp.m in csdms-contrib/slepian_alpha
% 
% Last Modified: January 8, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse through inputs
p=inputParser;
p.addRequired('sacfiles',@(x) iscellstr(x));
p.addParameter('procInfo',{},@(x) isempty(x) || (iscell(x) && length(x)==5));
p.addParameter('markArrs',{},@(x) isempty(x) || (iscell(x) && length(x)==3));
p.addParameter('arrivalArea',0,@(x) x==0 || x==1);
p.addParameter('SNRmethod',[],@(x) isempty(x) || isnumeric(x));
p.addParameter('saveDirectory','',@(x) isempty(x) || (ischar(x) && exist(x)==7));
p.addParameter('saveExtension','',@(x) isempty(x) || strcmp(x,'png') || strcmp(x,'pdf'));
p.addParameter('seeFigure',0,@(x) x==0 || x==1);
p.addParameter('makeFigure',1,@(x) x==0 || x==1);
p.addParameter('cartCoords',[],@(x) isempty(x) || isnumeric(x));
p.parse(sacfiles,varargin{:});
% 
parameters=p.Results;
procinfo=parameters.procInfo;
markarrs=parameters.markArrs;
arrarea=parameters.arrivalArea;
getsnr=parameters.SNRmethod;
savedir=parameters.saveDirectory;
saveext=parameters.saveExtension;
seefig=parameters.seeFigure;
makefig=parameters.makeFigure;
cartcoords=parameters.cartCoords;

% Default output if no arrivals or SNR values plotted.
snrtex='';

% Iterate through phases, if applicable
picktbl=[];
if ~isempty(markarrs)
  % Which phases to mark?
  markphase=markarrs{1};
  % Use SAC files or catalog?
  usecat=markarrs{2};
  % Narrow down time-axis window based on arrival?
  wintimes=markarrs{3};
  if length(markphase)>1
    winphase='P';
  else
    winphase=markphase;
  end
  
  %%%%%%%%%%%%%
  % Retrieve time of arrival around which to center SAC files
  
  % Arrivals in SAC files
  if usecat==0
    for s=1:length(sacfiles)
      sacfile=sacfiles{s};
      nowtbl=getpickinfo(sacfile,{winphase});
      if ~isempty(nowtbl)
        arrtime=nowtbl.ArrivalTime;
        break
      end
    end
      
  % Arrivals in catalog
  elseif usecat==1
    sacfile=sacfiles{1};
    [~,hdrinfo]=readsac(sacfile);
    nowevt=str2double(hdrinfo.KEVNM);
    nowsta=str2double(hdrinfo.KSTNM);
    sacvars={nowevt,nowsta,[]};
    nowtbl=getpickinfo(sacvars,{winphase});
    if ~isempty(nowtbl)
      arrtime=nowtbl.ArrivalTime;
    end
  end
  
  
  %%%%%%%%%%%%%
  % Retrieve information about arrivals to mark 
  
  % Arrivals in SAC files
  if usecat==0
    picktbl=[];
    for s=1:length(sacfiles)
      sacfile=sacfiles{s};
      if length(markphase)==1
        nowtbl=getpickinfo(sacfile,{markphase});
      else
        nowtbl=getpickinfo(sacfile,markphase);
      end
      picktbl=vertcat(picktbl,nowtbl);
    end
      
  % Arrivals in catalog
  elseif usecat==1
    if length(markphase)==1
      picktbl=getpickinfo(sacvars,{markphase});
    else
      picktbl=getpickinfo(sacvars,markphase);
    end
  end
end


% Extract information on calculating SNR
if ~isempty(getsnr)
  snrtype=getsnr(1);
  if snrtype==4
    threshval=getsnr(2);
    snrwin=getsnr(3);
    snrdef=[snrtype threshval];
  else
    snrwin=getsnr(2);
    snrdef=snrtype;
  end
end


% Retrieve event information
tempfile=sacfiles{1};
evtinfo=getevtinfo(tempfile);
if isempty(evtinfo)
  keyboard
end
evtnum=evtinfo.EVNT;
evtmag=evtinfo.MAG;
evtdep=evtinfo.DEP_KM;
% Origin time of event
origtime=getevttime(evtnum);

% Retrieve station information
stainfo=getstainfo(tempfile);
if isempty(stainfo)
  keyboard
end
stanum=stainfo.StationID;

% Calculate distance and azimuth (Cartesian or Great Circle)
% Great Circle
if isempty(cartcoords)
  evtlon=evtinfo.LON_dec;
  evtlat=evtinfo.LAT_dec;
  evcoord=[evtlon, evtlat];
  stalon=stainfo.Longitude;
  stalat=stainfo.Latitude;
  stacoord=[stalon,stalat];
  calcmode=2;

% Cartesian
else
  % Retrieve Cartesian coordinates
  evcoords=cartcoords(1,:);
  evcoord=evcoords(1:2);
  stacoords=cartcoords(2,:);
  stacoord=stacoords(1:2);
  calcmode=1;
end
[az,dist]=azdist(evcoord,stacoord,calcmode);


% If applicable, filter and/or deconvolve SAC files
if ~isempty(procinfo)
  % Processing information
  tfreqs=procinfo{2};
  pzfile=procinfo{3};
  valtype=procinfo{4};
  ffreqs=procinfo{5};
  if procinfo{1}==0
    dcinfo={tfreqs;pzfile;valtype};
    newsacs=procsacs(sacfiles,dcinfo,ffreqs);
    sacfiles=newsacs;
  end
end

% Get bounds for y-axis across all 3 components
if exist('arrtime','var')
  [maxval,minval]=sacbounds(sacfiles,wintimes,arrtime);
else
  [maxval,minval]=sacbounds(sacfiles,[],[]);
end
% Axis scaling!
if maxval>=1000 || minval<=-1000
  scaleval=1000;
  scalestr='k';
  maxval=maxval/1000;
  minval=minval/1000;
else
  scaleval=1;
  scalestr='';
end

% Assign bounds for time (x) axis
if exist('arrtime','var')
  % ...based on time of a plotted, picked arrival
  mintime=arrtime-seconds(wintimes(1));
  maxtime=arrtime+seconds(wintimes(2));
  mintick=seconds(mintime-origtime);
  maxtick=seconds(maxtime-origtime);
else
  % or common start and end times across all SAC files
  file1=sacfiles{1};
  [ustarttime,utims]=getsactime(file1);
  ufintime=ustarttime+seconds(max(utims));
  mintick=seconds(ustarttime-origtime);
  maxtick=seconds(ufintime-origtime);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize figure
if makefig==1
  if seefig==0
    seisfig=figure('visible','off');
  else
    seisfig=figure();
  end
  seisfig.Units='normalized';
  seisfig.Position(3)=0.55;
  seisfig.Position(4)=0.75;
else
  seisfig=gcf;
end


% Locations and sizes of all seismogram axes
numcomps=length(sacfiles);
switch numcomps
  case 1
    % nowax.Position(2)
    startpos=[0.12];
    % nowax.Position(4)
    heightpos=[0.75];
  case 2
    startpos=[0.5 0.15];
    heightpos=[0.35 0.35];  
  case 3
    startpos=[0.6 0.4 0.2];
    heightpos=[0.2 0.2 0.2];
  case 4
    startpos=[.675 0.5 0.325 0.15];
    heightpos=[0.175 0.175 0.175 0.175];    
  case 5
    startpos=[0.71 0.57 0.43 0.29 0.15];
    heightpos=[0.14 0.14 0.14 0.14 0.14]; 
end


% Load seismic data from each file and plot!
allax=[];
comps={};
for c=1:numcomps
  % Read SAC file
  sacfile=sacfiles{c};
  [sacdata,hdrinfo,~,~,tims]=readsac(sacfile);
  % Scale of data
  sacdata=sacdata./scaleval;
  % Component
  comp=hdrinfo.KCMPNM;
  comps=vertcat(comps,comp);
  
  % Align time (x) axis so that the event origin corresponds to 0
  sacstart=getsactime(sacfile);
  sacdiff=seconds(sacstart-origtime);
  tims=tims+sacdiff;  
  
  % Initialize axes for seismogram
  nowax=axes;
  allax=vertcat(allax,nowax);
  nowax.Units='normalized';
  hold on
  nowax.Position(2)=startpos(c);
  nowax.Position(4)=heightpos(c);
  % X-axis limits
  nowax.XLim=[mintick maxtick];
  % No x-axis label unless bottommost plot
  if c<numcomps
    nowax.XTickLabel={};
  else
    starttimestr=datestr(origtime);
    nowax.XLabel.String=sprintf('Time (s) since %s UTC (Event Origin)',...
      starttimestr); 
  end
  % Y-axis limits and label
  nowax.YLim=[round(1.2*minval,2,'significant')... 
    round(1.2*maxval,2,'significant')];
  nowax.YTick=[round(minval,2,'significant'), 0,...
    round(maxval,2,'significant')];


  % Check data type
  datatype=hdrinfo.IDEP;
  % Unknown
  if datatype==-12345
    nowax.YLabel.String=sprintf('%s (%scts)',comp,scalestr);
  % Displacement
  elseif datatype==6
    nowax.YLabel.String=sprintf('%s (%snm)',comp,scalestr);
    valstr='D';
  % Velocity
  elseif datatype==7
    nowax.YLabel.String=sprintf('%s (%snm/s)',comp,scalestr);
    valstr='V';
  % Acceleration
  elseif datatype==8
    nowax.YLabel.String=sprintf('%s (%snm/s^2)',comp,scalestr);
    valstr='A';
  end
  
  % Plot zero horizontal line
  zeroline=plot(tims,zeros(length(tims),1));
  zeroline.Color=[0.5 0.5 0.5];
  zeroline.LineStyle='--';
  zeroline.HandleVisibility='off';
  hold on
    
  % Plot seismic data
  nowplot=plot(tims,sacdata);
  nowplot.Color=[0 0 0];
  nowplot.HandleVisibility='off';
  hold on

  
  % Plot P and/or S wave picks, if desired
  if ~isempty(picktbl)
    [numpicks,~]=size(picktbl);

    % Legend label
    lgdstrs={};
    % Phases of arrivals
    arrphases={};
    % SNR values
    snrvals=[];
    % Time windows of all arrivals for SNR calculation
    snrwins={};
    % Plot colors of SNR plots
    plotclrs=[];
    % Legend for area below peaks, if necessary?
    areastrs={};
    % defout outputs
    definfo=[];
    
    % Iterate through each pick
    for p=1:numpicks
        
      pickrow=picktbl(p,:);
      % Arrival quality, type, time, and component
      phase=pickrow.ArrivalType{1};
      arrphases=vertcat(arrphases,phase);
      nowarrtime=pickrow.ArrivalTime;
      nowarrtim=seconds(nowarrtime-origtime);
      pickcomp=pickrow.Component{1};
      
      % Duller plotting color if the arrival was NOT picked on the 
      % current plotted component
      if ~strcmp(pickcomp,comp) || (strcmp(comp,'R') || strcmp(comp,'T'))
        if strcmpi(phase,'P')
          plotclr=[0.45 0 0.3];
          plotclrs=vertcat(plotclrs,plotclr);
        elseif strcmpi(phase,'S')
          plotclr=[0.45 0.3 0];
          plotclrs=vertcat(plotclrs,plotclr);
        end
      % Brighter color if the arrival is picked on the plotted component
      else
        if strcmpi(phase,'P')
          plotclr=[1 0 0.75];
          plotclrs=vertcat(plotclrs,plotclr);    
        elseif strcmpi(phase,'S')
          plotclr=[1 0.5 0];
          plotclrs=vertcat(plotclrs,plotclr); 
        end
      end
      
      % Plot arrival time as a vertical line
      arrline=plot([nowarrtim nowarrtim],[min(nowax.YLim) max(nowax.YLim)]);
      arrline.LineStyle='--';
      arrline.HandleVisibility='on';
      arrline.Color=plotclr;
      hold on
      
      % Plot area below arrival if applicable
      if arrarea==1
        [areaval,time1,time2]=phasearea(sacfile,phase,snrwin,threshval);
        % Points corresponding to both ends of the peak
        tim1=seconds(time1-origtime);
        tim2=seconds(time2-origtime);
        pktims=tims(tims>=tim1 & tims<=tim2);
        pkvals=sacdata(tims>=tim1 & tims<=tim2).';
        % If earlier time is interpolated
        if min(tims)~=tim1
          pktims=[tim1, pktims];
          pkvals=[0, pkvals];
        end
        
        % If later time is interpolated to 0
        if max(tims)~=tim2
          pktims=[pktims, tim2];
          pkvals=[pkvals, 0];
        end
        
        % Plot area
        pkplot=fill(pktims,pkvals,plotclr);
        pkplot.EdgeAlpha=0;
        pkplot.FaceAlpha=0.2;
        pkplot.FaceColor=plotclr;
        pkplot.HandleVisibility='off';
        
        % Add area value 
        if datatype==-12345
          areastr=sprintf('%s Area: %.2f ct^2',phase,areaval);     
        elseif datatype==6
          areastr=sprintf('%s Area: %.2f nm^2',phase,areaval);
        elseif datatype==7
          areastr=sprintf('%s Area: %.2f (nm/s)^2',phase,areaval);
        elseif datatype==8
          areastr=sprintf('%s Area: %.2f (nm/s^2)^2',phase,areaval);
        end
        areastrs=horzcat(areastrs,areastr);
      end
      
      % Add arrivals to legend
      lgdstr=sprintf('%s',phase);
      lgdstrs=horzcat(lgdstrs,lgdstr);
      
      % Add areas below peak to plot
      areatex=text(1,1,areastrs);
      areatex.Units='normalized';
      areatex.FontSize=8;
      areatex.Position(1)=0.02;
      areatex.Position(2)=0.95;
      areatex.HandleVisibility='off';
      hold on
      
      % Calculate SNR and get time window
      if ~isempty(getsnr)
        [snrval,snrtims,defout]=seis2snr(sacfile,nowarrtime,snrdef,snrwin);
        defout=defout./scaleval;
        snrtims=vertcat(snrtims,nowarrtime);
        snrtims={snrtims};
        snrvals=vertcat(snrvals,snrval);
        snrwins=vertcat(snrwins,snrtims);
        definfo=vertcat(definfo,defout);
      end
    end
    
    % Set axes limits again
    nowax.XLim=[mintick maxtick];
    nowax.YLim=[round(1.2*minval,2,'significant')... 
     round(1.2*maxval,2,'significant')];
    nowax.YTick=[round(minval,2,'significant'), 0,...
      round(maxval,2,'significant')];
    
    % If requested, plot a visual showing SNR of arrivals
    if ~isempty(getsnr)
      snrtex=plotsnr(snrvals,snrtype,snrwins,arrphases,nowax,origtime,...
        definfo,plotclrs);
      snrtex.HandleVisibility='off';
    end

    % Plot legend
    plotlgd=legend(nowax,lgdstrs);
    plotlgd.Location='northeast';
    

  % No picks: Just plot seismic data
  else
    nowplot=plot(tims,sacdata);
    nowplot.Color=[0 0 0];
    nowplot.HandleVisibility='off';
    hold on
  end
  
  % Set axes limits again
  nowax.XLim=[mintick maxtick];
  nowax.YLim=[round(1.2*minval,2,'significant')... 
    round(1.2*maxval,2,'significant')];
  nowax.YTick=[round(minval,2,'significant'), 0,...
    round(maxval,2,'significant')];
  

  % Plot Title
  if c==1
    % Event and station
    titlestr1=sprintf('Event ID %d: M = %.2f, Depth = %.3f km',evtnum,...
      evtmag,evtdep);
    titlestr2=sprintf('Station %d:',stanum);
    deltastr=' \Delta';
    % Get P-S distance
    if ~isempty(cartcoords)
      psdist=getpsdist(evtnum,stanum,0);
    else
      psdist=getpsdist(evtnum,stanum,1); 
    end
    
    
    
    % Distance and azimuth
    if isempty(psdist)
      if isempty(cartcoords)
        deltastr2=sprintf('=%.2f%s, Azimuth = %.2f%s',...
          round(dist,2),char(176),round(az,2),char(176));
      else
        deltastr2=sprintf('=%.2f km, Azimuth = %.2f%s',...
          round(dist,2),round(az,2),char(176));
      end
    else
      if isempty(cartcoords)
        deltastr2=sprintf('=%.2f%s, Azimuth = %.2f%s, PS = %.2f%s',...
          round(dist,2),char(176),round(az,2),char(176),...
          round(psdist,2),char(176));
      else
        deltastr2=sprintf('=%.2f km, Azimuth = %.2f%s, PS = %.2f km',...
          round(dist,2),round(az,2),char(176),round(psdist,2));
      end
    end
    titlestr2=strcat(titlestr2,deltastr,deltastr2);
    % No deconvolution or filtering
    if isempty(procinfo)
      nowax.Title.String={titlestr1;titlestr2};
      
    % Both deconvlution and filtering
    elseif ~isempty(tfreqs) && ~isempty(ffreqs)
      titlestr3=sprintf('Deconvolution: [%.2f %.2f %.2f %.2f], Bandpass %.2f-%.2f Hz',...
        tfreqs(1),tfreqs(2),tfreqs(3),tfreqs(4),ffreqs(1),ffreqs(2));
      nowax.Title.String={titlestr1;titlestr2;titlestr3};
        
    % Deconvolved but not filtered
    elseif ~isempty(tfreqs) && isempty(ffreqs)
      titlestr3=sprintf('Deconvolution: [%.2f %.2f %.2f %.2f] Hz',...
        tfreqs(1),tfreqs(2),tfreqs(3),tfreqs(4));
      nowax.Title.String={titlestr1;titlestr2;titlestr3};
    
    % Filtered but not deconvolved
    elseif ~isempty(ffreqs) && isempty(tfreqs)
      titlestr3=sprintf('Bandpass %.2f-%.2f Hz',ffreqs(1),ffreqs(2));
      nowax.Title.String={titlestr1;titlestr2;titlestr3};
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure name
if ~isempty(savedir)
  % Check save extension
  if isempty(saveext)
    figname='';
    return
  end
  
  % Components
  compstr='';
  for c=1:length(comps)
    compstr=sprintf('%s%s',compstr,comps{c});
  end
  
  % Arrivals or not
  if ~isempty(markphase)
    if length(markphase)==1
      arrstr=sprintf('.%s',markphase);
    else
      arrstr=sprintf('.PS');
    end
  else
    arrstr='';
  end

  % No deconvolution or filtering
  if isempty(procinfo)
    procstr='';
      
  % Both deconvlution and filtering
  elseif ~isempty(tfreqs) && ~isempty(ffreqs)
    procstr=sprintf('.%s.DC%.1f%.1f%.1f%.1f_FI%.1f%.1f',valstr,...
      tfreqs(1),tfreqs(2),tfreqs(3),tfreqs(4),ffreqs(1),ffreqs(2));
        
  % Deconvolved but not filtered
  elseif ~isempty(tfreqs) && isempty(ffreqs)
    procstr=sprintf('.%s.DC%.1f%.1f%.1f%.1f',valstr,...
      tfreqs(1),tfreqs(2),tfreqs(3),tfreqs(4));
    
  % Filtered but not deconvolved
  elseif ~isempty(ffreqs) && isempty(tfreqs)
    procstr=sprintf('.FI%.1f%.1f',ffreqs(1),ffreqs(2));
  end
  
  % SNR NOT calculated
  if isempty(getsnr)
    snrstr='';
  else
    % Amplitude
    if snrwin==1
      snrstr=sprintf('.SNRAMP%.2f',snrtype);
        
    % Variance
    elseif snrwin==2
      snrstr=sprintf('.SNRVAR%.2f',snrtype);  
        
    % Root mean square
    elseif snrwin==3
      snrstr=sprintf('.SNRRMS%.2f',snrtype);
    end
  end
  
  % Full figure name
  fname=sprintf('EVT%d.STA%d.SEIS.%s%s%s%s',evtnum,stanum,compstr,...
    arrstr,procstr,snrstr);
   
  % Save figure
  if strcmp(saveext,'png')
    figname=sprintf('%s.%s',fname,'png');
    print(figname,'-dpng');
  elseif strcmp(saveext,'pdf')
    setenv('EPS',savedir);
    fname=sprintf('%s.eps',fname);
    figname=figdisp(fname,[],[],2,[],'epstopdf');
  end
else
  figname='';
end


