function [ufig,uax,uname]=udistr(inptmat,issynth,makefig,clrscheme,...
  addlbl)
% 
% Function to plot the distribution of measured or synthetic areas under 
% the P-wave at various stations about a seismic event. 
% 
% INPUTS
% inptmat : Full path to the .mat file resulting from moment tensor 
%           inversion using fociMT, containing the input parameters 
%           of the inversion. 
%
% issynth : Are these areas synthetic (i.e. derived from a run with 
%           synthetic data, or real data)?
%           0 : Measured 
%           1 : Synthetic 
% 
% makefig : Do we make a new figure to plot our distribution, or do we 
%           plot on an existing figure? 
%           0 : Plot on existing figure
%           1 : Plot on new figure
% 
% clrscheme : What color scheme do we use to mark stations with upwards
%             P-wave polarity (corresponding to compressional quadrants 
%             in the focal mechanism) and downwards P-wave polarity
%             (corresponding to dilatational quadrants in the focal 
%             mechanism)?
% 
%             1 : Dark grey for upwards polarity, white for downwards 
%                 polarity, grey for stations whose signals are too small 
%                 for a definite polarity [Default]
% 
%             2 : Red for upwards polarity, blue for downwards polarity,
%                 lighter shade for stations whose signals are too small 
%                 for a definite polarity 
% 
%             3 : Blue for upwards polarity, red for downwards polarity,
%                 lighter shade for stations whose signals are too small 
%                 for a definite polarity, consistent with the 
%                 online application pyfocmec at 
%                 https://share.streamlit.io/cedrictwardz/pyfocmec/main/
%                 app.py?fbclid=IwAR3IJlDcd_p5gtPoo9ppLOK2dL5oT_pL65y1cM0dml8_-qFQ7pFCpBdY7Os
% 
% addlbl : Do we add labels to each station, to indicate the station 
%          number? 
%          0 : No [Default]
%          1 : Yes
% 
% OUTPUTS
% ufig : Figure handle of azimuthal plot of areas 
% uax : Axis handle
% uname : Full path to ufig, saved in the same directory as the
%         corresponding moment tensor solution. The figure is saved only 
%         if makefig = 1. Otherwise, uname is ''
% 
% References : 
% Kwiatek, G., Martínez‐Garzón, P. & Bohnhoff, M. (2016) HybridMT: A 
% MATLAB/Shell Environment Package for Seismic Moment Tensor Inversion and 
% Refinement. Seismological Research Letters, 87, 964–976. 
% doi:10.1785/0220150251
% 
% Uses defval.m and figdisp.m in csdms-contrib/slepian_alpha
% 
% Last Modified : March 10, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default values
defval('clrscheme',1);
defval('addlbl',0);

% Load .MAT files
load(inptmat);

% Event information + number of stations
evtnum=str2double(Input{1}.event_id);
evtn=0.001*Input{1}.e_northing;
evte=0.001*Input{1}.e_easting;
numsta=length(Input{1}.S_EASTING);

% Obtain farthest coordinates between event and station
staes=0.001*Input{1}.S_EASTING-evte;
stans=0.001*Input{1}.S_NORTHING-evtn;
maxcoord=1.1*max(abs(vertcat(staes,stans)));

% Obtain a lower threshold for area under the P-wave
allomegas=Input{1}.OMEGA;
omegathresh=0.01*max(abs(allomegas));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot distribution of areas

% Initialize figure
if makefig==1
  ufig=figure();
  ufig.Units='normalized';
  ufig.Position(1)=0.15;
  ufig.Position(2)=0.2;
  ufig.Position(3)=0.45;
  ufig.Position(4)=0.65;
  
elseif makefig==0
  ufig=gcf;
end

% Initialize axes
uax=axes;
pbaspect([1 1 1])
hold on
% 
uax.XLabel.String='Relative Easting (km)';
uax.YLabel.String='Relative Northing (km)';
if issynth==0
  uax.Title.String=sprintf('Event %d, Areas Under the P-wave',evtnum);
elseif issynth==1
  uax.Title.String=sprintf('Event %d, Synthetic Areas Under the P-wave',...
    evtnum); 
end

% Plot dotted lines for X = 0 and Y = 0
x0line=plot([0 0],[-1*maxcoord maxcoord]);
x0line.LineWidth=1;
x0line.Color=[0.5 0.5 0.5];
x0line.LineStyle='--';
x0line.HandleVisibility='off';
hold on
% 
y0line=plot([-1*maxcoord maxcoord],[0 0]);
y0line.LineWidth=1;
y0line.Color=[0.5 0.5 0.5];
y0line.LineStyle='--';
y0line.HandleVisibility='off';

% Plot 'N' to signify where North is
nmarker=text(0.03*maxcoord,0.93*maxcoord,'N');
nmarker.Color=[0.3 0.3 0.3];
nmarker.FontSize=12;
nmarker.HandleVisibility='off';

% Plot marker for event
evtplot=scatter(0,0);
evtplot.MarkerEdgeColor=[0 0 0];
evtplot.MarkerFaceColor=[0 0 0];
evtplot.Marker='+';
evtplot.LineWidth=4;
evtplot.SizeData=220;
evtplot.HandleVisibility='off';

% Iterate through each measured moment and plot! 
for s=1:numsta
  % Station coordinates and moment
  omegaval=Input{1}.OMEGA(s);
  stae=0.001*(Input{1}.S_EASTING(s))-evte;
  stan=0.001*(Input{1}.S_NORTHING(s))-evtn;

  % Color and shape of marker, based on polarity and the threshold value 
  % of the area under the P-wave
  if abs(omegaval)>omegathresh
    % Upwards polarity
    if omegaval>0
      if clrscheme==1
        mclr=[0.2 0.2 0.2];
      elseif clrscheme==2
        mclr=[1 0 0];
      elseif clrscheme==3
        mclr=[0 0.3 1];  
      end
      mmarker='^';
      msize=100;
      
    % Downwards polarity
    else
      if clrscheme==1
        mclr=[1 1 1];
      elseif clrscheme==2
        mclr=[0 0.3 1];
      elseif clrscheme==3
        mclr=[1 0 0];  
      end
      mmarker='v';
      msize=100;
    end
    
  % Signal too small for definite polarity
  else
    if clrscheme==1
      mclr=[0.75 0.75 0.75];
      mmarker='o';
      msize=36;
    elseif clrscheme==2
      if omegaval>0
        mclr=[1 0.55 0.55];
        mmarker='^';
      else
        mclr=[0.55 0.55 1];  
        mmarker='v';
      end
      msize=50;
    elseif clrscheme==3
      if omegaval>0
        mclr=[0.55 0.55 1]; 
        mmarker='^';
      else
        mclr=[1 0.55 0.55];  
        mmarker='v';
      end
      msize=50;
    end
  end
  
  % Plot station and measured moment
  omegaplot=scatter(stae,stan);
  hold on
  omegaplot.MarkerEdgeColor=[0 0 0];
  omegaplot.LineWidth=1;
  omegaplot.MarkerFaceColor=mclr;
  omegaplot.Marker=mmarker;
  omegaplot.SizeData=msize;
  
  % Plot station number?
  if addlbl==1
    nowlbl=Input{1}.Station{s};
    stalbl=text(stae,stan,nowlbl);
    stalbl.HorizontalAlignment='center';
    stalbl.VerticalAlignment='top'; 
  end
end


% Axes limits
uax.XLim=[-1*maxcoord maxcoord];
uax.YLim=[-1*maxcoord maxcoord];

% Font Size
uax.FontSize=12;

% Save Figure
if makefig==1
  % Directory where we save our figure
  [savedir,~,~]=fileparts(inptmat);
    
  % Name of figure
  if issynth==0
    uname=sprintf(...
      'OMEGADISTR_%s_EVT%d_DISP.eps',solstr,evtnum);
  elseif issynth==1
    uname=sprintf(...
      'OMEGADISTR_SYNTH%s_EVT%d_DISP.eps',solstr,...
      evtnum);
  end
  
  % Save figure
  setenv('EPS',savedir);
  uname=figdisp(uname,[],[],2,[],'epstopdf');
else
  uname='';
end


