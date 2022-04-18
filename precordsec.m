function [recfig,recax,recname]=precordsec(evtnum,procinfo,threshval,...
  noisewin,varargin)
% 
% Function to plot a record section of the P-waves of a given event
% (recorded in the Tajik Basin by SeismoTech Ltd.), recorded by 
% stations and with a signal to noise ratio at or above a given
% threshold value. Note that all seismograms are in the vertical component
% and will be deconvolved to the displacement-time domain. 
% 
% INPUTS
% evtnum : ID number corresponding to the event
% procinfo : To obtain displacement SAC files, input information needed to
%            deconvolve our raw SAC files: 
% 
%            Format: {a; b; c; d} where
% 
%            a : Corner frequencies in Hz, through which we deconvolve our 
%                data. Define a four-element vector with frequencies f1, 
%                f2, f3, f4 in Hz, where f2 is at least twice f1, and f4 
%                at least twice f3. Deconvolve the seismograms if we enter 
%                this vector. If we don't wish to do so, enter an empty 
%                vector. 
% 
%            b : Leave empty. Different stations may require different 
%                polezero files with which to deconvolve, which will be 
%                handled in the function. 
% 
%            c : Input 0 to deconvolve into displacement.
% 
%            d : Frequencies through which we filter our data. Define a 2 
%                element vector, where the first one is our lower bound 
%                and the second one is our upper bound. Both in Hz. 
% 
% threshval : To find the peak corresponding to the first-arriving P-wave,
%             we imposed a threshold such that the amplitude of this peak
%             should exceed the mean absolute amplitude of the noise by 
%             this factor. Input this number (see phasepk.m)
% 
% noisewin : To find the peak corresponding to the first-arriving P-wave,
%            we compare the amplitude of the P-wave to the mean absolute
%            amplitude of the noise. Enter the number of seconds before
%            the onset of the P-wave (as determined by SeismoTech Ltd.)
%            to use as the noise. 
% 
% Input the following as varargin: 
% 
% 'mint' : How many seconds prior to the P-wave arrival do we wish to begin 
%          our record section? Default: 0.5
% 
% 'maxt' : How many seconds after the P-wave arrival do we wish to end our 
%          record section? Default: 1
%
% 'scaleopt' : How do we normalize our seismograms? 
%              1 : Normalize each one relative to its maximum amplitude 
%                  that we plot [Default]
%              2 : Normalize all seismograms relative to the maximum 
%                  amplitude that we plot, across all seismograms
% 
% 'scaleval' : To what value should we normalize the maximum amplitude? Default: 10
% 
% 'shownum' : How many seismograms do we show? Default: 20
%             To show all seismograms, enter []. 
%             Note that plotting a subset prioritizes the seismograms 
%             with the highest P-wave SNR
%             
% 'seisorder' : How do we order the seismograms? 
%               1. By Azimuth
%               2. Arbitrarily, with even spacing between each one
% 
% 'makefig' : Do we plot this record section on a new figure, or on an 
%             existing figure? 
% 
%             0 : On an existing figure
%             1 : Make a new figure
% 
% 'savedir' : Full path to the directory where we wish to save our record
%             section. Input '' if we do not wish to save our figure
% 
% 'clrscheme' : What color scheme do we use to plot our seismograms?
% 
%               1 : Dark grey for all seismograms [Default]
%   
%               2 : Red for upwards P-waves, blue for downwards P-waves
% 
%               3 : Blue for upwards P-waves, red for downwards P-waves
% 
% 'addhorz' : Do we wish to add a faint line indicating the "zero" for 
%             each seismogram? 
%             0 : No [default]
%             1 : Yes
% 
% 'labelp' : Do we wish to mark the peak corresponding to the P-wave
%            with a dot?
%            0 : No [Default]
%            1 : Yes
% 
% 
% OUTPUTS
% recfig : Figure handle to record section
% recax : Axis handle containing the record section
% recname : Full path to saved record section. '' if figure is not saved
%    
% Last Modified : March 15, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse through inputs
p=inputParser;
p.addRequired('evtnum',@(x) isnumeric(x));
p.addRequired('procinfo',@(x) (iscell(x) && length(x)==4))
p.addRequired('threshval',@(x) isnumeric(x));
p.addRequired('noisewin',@(x) isnumeric(x));
p.addParameter('mint',0.5,@(x) isnumeric(x));
p.addParameter('maxt',1,@(x) isnumeric(x));
p.addParameter('scaleopt',1,@(x) isnumeric(x));
p.addParameter('scaleval',10,@(x) isnumeric(x));
p.addParameter('shownum',20,@(x) isnumeric(x));
p.addParameter('seisorder',1,@(x) isnumeric(x));
p.addParameter('makefig',0,@(x) isnumeric(x));
p.addParameter('savedir','',@(x) isnumeric(x));
p.addParameter('clrscheme',1,@(x) isnumeric(x));
p.addParameter('addhorz',0,@(x) isnumeric(x));
p.addParameter('labelp',0,@(x) isnumeric(x));
p.parse(evtnum,procinfo,threshval,noisewin,varargin{:});
% 
parameters=p.Results;
mint=parameters.mint;
maxt=parameters.maxt;
scaleopt=parameters.scaleopt;
scaleval=parameters.scaleval;
shownum=parameters.shownum;
seisorder=parameters.seisorder;
makefig=parameters.makefig;
savedir=parameters.savedir;
clrscheme=parameters.clrscheme;
addhorz=parameters.addhorz;
labelp=parameters.labelp;

% Initialize figure
if makefig==1
  recfig=figure();
  recfig.Units='normalized';
  recfig.Position(1)=0.1;
  recfig.Position(2)=0.05;
  recfig.Position(3)=0.65;
  recfig.Position(4)=0.75;
elseif makefig==0
  recfig=gcf;
end
recax=axes;
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Event coordinates
evtrow=getevtinfo(evtnum);
% Easting and Northing coordinates, in m
evte=evtrow.EASTING;
evtn=evtrow.NORTHING; 

% Retrieve P-wave picks in the Z-component!
sacvars={evtnum;[];'Z'};
arrtypesin={'P'};
picktbl=getpickinfo(sacvars,arrtypesin);
[numpicks,~]=size(picktbl);


% Iterate through each pick
snrvals=[];
sacfiles={};
stanums=[];
mintplot=[];
maxtplot=[];
for p=1:numpicks
  % Arrival information
  pickrow=picktbl(p,:);
  nowsta=pickrow.Station;
  nowcomp=pickrow.Component{1};
  nowarr=pickrow.ArrivalType{1};
  if ~strcmp(nowcomp,'Z')
    continue
  end
  if ~strcmp(nowarr,'P')
    continue
  end
  
  % Time of onset of P-wave, as chosen by SeismoTech Ltd. 
  starrtime=pickrow.ArrivalTime;
    
  % Get polezero file for this station for deconvolution
  pzfile=getpzfile(nowsta);
  procinfo{2}=pzfile;
  
  % Obtain SAC file
  sacfile=getsac_TOTAL(evtnum,nowsta,{'Z'},procinfo);
  sacfile=sacfile{1};
  
  % Calculate the SNR of the first-arriving P-wave, if necessary
  [snrval,~,~]=seis2snr(sacfile,starrtime,[4 threshval],noisewin);
  % Skip this P-wave arrival if its signal is too weak OR 
  % if the SNR is empty, meaning the SAC file has no meaningful data
  if isempty(snrval) || abs(snrval)<abs(threshval)
    continue; 
  end
  
  % Add to list of SAC files
  sacfiles=vertcat(sacfiles,sacfile);
  
  % Add to list of SNR files
  snrvals=vertcat(snrvals,snrval);
  
  % Add to list of station numbers
  stanums=vertcat(stanums,nowsta);
  
  % Retrieve calculated time of onset of the P-wave
  [~,time1,~,~]=phasearea(sacfile,starrtime,noisewin,threshval);
    
  % Retrieve data immediately surrounding the onset of the P-wave
  [~,~,~,~,tims]=readsac(sacfile);
  sacst=getsactime(sacfile);
  sactims=sacst+seconds(tims);
  ptims=tims(...
    sactims>=(time1-seconds(mint)) & sactims<=(time1+seconds(maxt))); 
  ptims=ptims-min(ptims);
  ptims=ptims-mint;

  % Limits for time axis
  if isempty(mintplot)
    mintplot=abs(min(ptims));
  else
    nowmint=abs(min(ptims));
    if nowmint<mintplot
      mintplot=nowmint;
    end
  end
  if isempty(maxtplot)
    maxtplot=abs(max(ptims));
  else
    nowmaxt=abs(max(ptims));
    if nowmaxt<maxtplot
      maxtplot=nowmaxt;
    end
  end
end


% If plotting only a subset of P-wave arrivals, select the ones with 
% highest SNR
if ~isempty(shownum)
  ptbl=table(abs(snrvals),sacfiles,stanums,snrvals);
  ptbl=sortrows(ptbl,1,'descend');
  [numr,~]=size(ptbl);
  if numr<shownum
    shownum=numr;
  end
  snrvals=ptbl.snrvals(1:shownum);
  stanums=ptbl.stanums(1:shownum);
  sacfiles=ptbl.sacfiles(1:shownum);
end
numfiles=length(sacfiles);

% Retrieve maximum amplitude across all plotted seismograms, if needed
maxamp=[];
if scaleopt==2
  for f=1:numfiles
    nowsta=stanums(f);
    nowsac=sacfiles{f};
    pickrow=picktbl(picktbl.Station==nowsta,:);
    [nump,~]=size(pickrow);
    if nump>1
      keyboard
    end
    starrtime=pickrow.ArrivalTime;
      
    % Retrieve calculated time of onset of the P-wave
    [~,time1,~,~]=phasearea(nowsac,starrtime,noisewin,threshval);
    
    % Retrieve data immediately surrounding the onset of the P-wave
    [sacdata,~,~,~,tims]=readsac(nowsac);
    sacst=getsactime(nowsac);
    sactims=sacst+seconds(tims);
    pdata=sacdata(...
      sactims>=(time1-seconds(mint)) & sactims<= (time1+seconds(maxt)));
    
    if isempty(maxamp)
      maxamp=max(abs(pdata));
    else
      nowamp=max(abs(pdata));
      if nowamp>maxamp
        maxamp=nowamp;
      end
    end
  end
end


% Plot the selected seismograms
for f=1:numfiles
  nowsta=stanums(f);
  nowsac=sacfiles{f};
  nowsnr=snrvals(f);
  pickrow=picktbl(picktbl.Station==nowsta,:);
  [nump,~]=size(pickrow);
  if nump>1
    keyboard
  end
  starrtime=pickrow.ArrivalTime;
    
  % Retrieve calculated time of onset of the P-wave
  [~,time1,~,~]=phasearea(nowsac,starrtime,noisewin,threshval);
    
  % Retrieve data immediately surrounding the onset of the P-wave
  [sacdata,~,~,~,tims]=readsac(nowsac);
  sacst=getsactime(nowsac);
  sactims=sacst+seconds(tims);
  ptims=tims(...
    sactims>=(time1-seconds(mint)) & sactims<= (time1+seconds(maxt))); 
  ptims=ptims-min(ptims);
  ptims=ptims-mint;
  pdata=sacdata(...
    sactims>=(time1-seconds(mint)) & sactims<= (time1+seconds(maxt)));

  % Normalize amplitude
  if scaleopt==1
    nowamp=max(abs(pdata));
    pdata=scaleval*(pdata/nowamp);
    
  elseif scaleopt==2
    pdata=scaleval*(pdata/maxamp);
  end
  
  if seisorder==1
    % Calculate azimuth
    starow=getstainfo(nowsta);
    stae=starow.Easting;
    stan=starow.Northing;
    [azval,~]=azdist([evte evtn],[stae stan],1);
    % Add azimuth to vertical coordinate
    pdata=pdata+azval;
      
  elseif seisorder==2
    % Separate seismograms by 10 units
    pdata=pdata+10*(f-1);
  end
  
  % Add a 'zero line' for each seismogram
  if addhorz==1
    if seisorder==1
      horzline=plot(ptims,azval*ones(length(ptims),1));
    elseif seisorder==2
      horzline=plot(ptims,10*(f-1)*ones(length(ptims),1));   
    end
    horzline.LineWidth=0.5;
    horzline.Color=[0.8 0.8 0.8];
  end
  
  % Plot seismogram
  try
    seisplot=plot(ptims,pdata);
  catch
    keyboard
  end
  
  % Color scheme
  if clrscheme==1
    seisplot.Color=[0.2 0.2 0.2];
    
  elseif clrscheme==2
    if nowsnr<0
      seisplot.Color=[0 0.3 1];
    else
      seisplot.Color=[1 0 0];
    end
    
  elseif clrscheme==3
    if nowsnr<0
      seisplot.Color=[1 0 0];
    else
      seisplot.Color=[0 0.3 1];
    end
  end
  seisplot.LineWidth=1;
  
  % Plot P-wave peak!
  if labelp==1
    % Retrieve time and amplitude of P-wave
    [pkval,pktime,~]=phasepk(nowsac,starrtime,noisewin,threshval);
    pktim=seconds(pktime-sacst);
    ptims2=tims(...
      sactims>=(time1-seconds(mint)) & sactims<=(time1+seconds(maxt))); 
    pktim=pktim-min(ptims2)-mint;
    
    % Scale amplitude of P-wave
    if scaleopt==1
      pkval=scaleval*(pkval/nowamp);
    elseif scaleopt==2
      pkval=scaleval*(pkval/maxamp);
    end
    if seisorder==1
      pkval=pkval+azval;
    elseif seisorder==2
      pkval=pkval+10*(f-1);
    end
    
    % Mark P-wave peak
    pkplot=plot(pktim,pkval);
    pkplot.Marker='.';
    pkplot.MarkerSize=20;
    if clrscheme==1
      pkplot.Color=[1 0.2 0.8];
    else
      pkplot.Color=[0 0 0];
    end
  end
end


% Adjust axes limits
recax.XLim=[-1*mintplot maxtplot];
if seisorder==1
  recax.YLim=[0 370];
  recax.YTick=0:45:360;
  recax.YLabel.String=sprintf('Azimuth (%s)',char(176));
elseif seisorder==2
  recax.YLim=[-10 10*numfiles];
  recax.YTick=[];
end


% Add line indicating P-wave onset
pline=plot([0 0],[min(recax.YLim) max(recax.YLim)]);
if clrscheme==1
  pline.Color=[1 0.2 0.8];
else
  pline.Color=[0.2 0.2 0.2];
end
pline.LineStyle='--';
pline.LineWidth=1.5;

% Axes Labels and Title
recax.XLabel.String='Time (s) wrt Onset of P-wave';
recax.FontSize=12;
recax.Title.String=sprintf('Record Section of P-wave Arrivals: Event %d',...
  evtnum);


% Save figure if requested AND we make a new figure
if exist(savedir)==7 && makefig==1
  ffreqs=procinfo{4};
  recname=sprintf(...
    'EVT%d_PRECORDSEC_%.2gto%.2gHz_SNRTHR%.2g_NOISEWIN%.2gs.eps',...
    evtnum,ffreqs(1),ffreqs(2),threshval,noisewin);
  setenv('EPS',savedir);
  recname=figdisp(recname,[],[],2,[],'epstopdf');
  recname=sprintf('%s.pdf',recname(1:end-4));
else
  recname='';
end


