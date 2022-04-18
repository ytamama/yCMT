function [mapfig,mapax,mapname]=focimt_map(evtnums,solntype,varargin)
% 
% Function to (1) invert for the moment tensor of a seismic event, 
% recorded by SeismoTech Ltd. and originating in the Tajik Basin, and 
% (2) plot the resulting beach ball diagram on a map of the Tajik Basin, 
% featuring fault lines and topography. 
% 
% INPUTS
% evtnums : ID number(s) of seismic event(s) plotted. 
% solntype : Which type of moment tensor do we plot? 
%            1 : Full
%            2 : Deviatoric
%            3 : Double-couple
% 
% Input the following as varargin: 
% 
% 'solmats' : Full path to the .mat file containing the solution to our 
%             inversion using fociMT, or a cell array containing the 
%             paths to each file as a string. 
% 
% 'inputfile' : Input file to use in our inversion for fociMT, or a cell 
%               array containing the paths to each file as a string
% 
% 'procinfo' : A 4 x 1 cell array telling us how to process our seismic 
%              data
%           
%              Format: {a; b; c; d} where
% 
%              a : Corner frequencies in Hz, through which we deconvolve 
%                  data. Define a four-element vector with frequencies f1, 
%                  f2, f3, f4 in Hz, where f2 is at least twice f1, and f4 
%                  at least twice f3. Deconvolve the seismograms if we enter 
%                  this vector. If we don't wish to do so, enter an empty 
%                  vector. 
% 
%              b : Leave empty. Different stations may require different 
%                  polezero files with which to deconvolve, which will be 
%                  handled in the function. 
% 
%              c : Input 0 to deconvolve into displacement.
% 
%              d : Frequencies through which we filter our data. Define a 2 
%                  element vector, where the first one is our lower bound 
%                  and the second one is our upper bound. Both in Hz. 
% 
% 'snrthresh' : To find the peak corresponding to the P-wave arrival,
%               we imposed a threshold such that the amplitude of this peak
%               should exceed the mean absolute amplitude of the noise by 
%               this factor. Input this number (see phasepk.m)
% 
% 'noisewin' : To find the peak corresponding to the P-wave,
%              we compare the amplitude of the P-wave to the mean absolute
%              amplitude of the noise. Enter the number of seconds before
%              the onset of the P-wave (as determined by SeismoTech Ltd.)
%              to use as the noise. 
% 
% 'rho' : What is the density of the Tajik Basin, in kg/m^3? By default, 
%         we use the density of the model of Chapman et al. 2017, 
%         doi: 10.1111/bre.12381, 2500 kg/m3.
% 
% 'mapzoom' : Do we wish to zoom our map into the vicinity of the events
%             whose moment tensors we plot? 
%             0 : No [Default]
%             # : Yes, add a 'buffer' of this many degrees longitude 
%                 and latitude on either side of the farthest events. 
% 
% 'clrscheme' : What color scheme do we use to mark stations with upwards
%               P-wave polarity (corresponding to compressional quadrants 
%               in the focal mechanism) and downwards P-wave polarity
%               (corresponding to dilatational quadrants in the focal 
%               mechanism)?
% 
%               1 : Dark grey for upwards polarity, white for downwards 
%                   polarity, grey for stations whose signals are too small 
%                   for a definite polarity [Default]
%   
%               2 : Red for upwards polarity, blue for downwards polarity,
%                   white for stations whose signals are too small for a 
%                   definite polarity 
% 
% 'mechdiam' : Diameter of the beach ball diagram, in degrees (longitude or
%              latitude). Default: 0.05
% 
% 'evtlbl' : Do we wish to attach a label to each beach ball, showing 
%            the event's ID number? 
%            0 : No
%            1 : Yes
% 
% 'makefig' : Do we wish to create a new figure for our map, or plot on an 
%             existing figure? 
%             0 : Plot on existing figure 
%             1 : Create a new figure handle [Default]
% 
% 'savedir' : Where do we save our resulting map? Input '' if not saving 
%             our map (default). Note that if makefig = 0, we will not 
%             save our resulting map in this function. 
% 
% 'axpos' : 4-element vector for the position of the axes in the current
%           figure, assuming the units of the axes are normalized. 
% 
% 'needcbar' : Do we need the colorbar indicating topography or no? 
%              0 : No
%              1 : Yes
% 
% 'axlbl' : Do we need axis labels and ticks, and a plot title? 
%           0 : No
%           1 : Yes
% 
% 
% OUTPUTS
% mapfig : Figure handle of the resulting map
% mapax : Axis handle of the resulting map
% mapname : Full path to the saved figure. Return '' if the figure is 
%           not saved. 
% 
% References : 
% Kwiatek, G., Martínez‐Garzón, P. & Bohnhoff, M. (2016) HybridMT: A 
% MATLAB/Shell Environment Package for Seismic Moment Tensor Inversion and 
% Refinement. Seismological Research Letters, 87, 964–976. 
% doi:10.1785/0220150251
% 
% Uses focalmech.m by James Conder (2022) to plot the beachball diagrams. 
% James Conder (2022). focalmech(fm, centerX, centerY, diam, varargin) 
% (https://www.mathworks.com/matlabcentral/fileexchange/61227-focalmech-fm-centerx-centery-diam-varargin), 
% MATLAB Central File Exchange. Retrieved February 14, 2022.
% 
% Conder, J.A. and C.A. Arciniegas, Conjugate Faulting in the Wabash
%     Valley Fault Zone exhibited by the Nov 20, 2012 M3.6 earthquake, 
%     a Mt Carmel Late Aftershock, Seismological Research Letters, 88,
%     1203-1209, doi:10.1785/0220170021, 2017
% 
% Topography data from GEBCO Compilation Group. (2021) GEBCO 2021 Grid. 
% doi:10.5285/c6612cbe-50b3-0cff-e053-6c86abc09f8f
% 
% Fault lines from Gagala et al. 2020 and Kufner et al. 2018
% 
% We set the default density of the Tajik Basin is 2500 kg/m3, based on a 
% model by Chapman et al. 2017, doi: 10.1111/bre.12381
% 
% Uses defval.m and figdisp.m in csdms-contrib/slepian_alpha
% 
% Last Modified: March 29, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse through inputs
p=inputParser;
p.addRequired('evtnums',@(x) isnumeric(x));
p.addRequired('solntype',@(x) isnumeric(x));
p.addParameter('solmats',{},@(x) isempty(x) || iscell(x) || exist(x)==2);
p.addParameter('inputfile',{},@(x) isempty(x) || iscell(x) || exist(x)==2);
p.addParameter('procinfo',{},@(x) isempty(x) ||...
  (iscell(x) && length(x)==4))
p.addParameter('snrthresh',2,@(x) isnumeric(x));
p.addParameter('noisewin',3,@(x) isnumeric(x));
p.addParameter('rho',2500,@(x) isnumeric(x));
p.addParameter('mapzoom',0,@(x) isnumeric(x) && x>=0);
p.addParameter('clrscheme',1,@(x) isnumeric(x));
p.addParameter('mechdiam',0.05,@(x) isnumeric(x) && x>=0);
p.addParameter('evtlbl',0,@(x) x==0 || x==1)
p.addParameter('makefig',1,@(x) x==0 || x==1);
p.addParameter('savedir','',@(x) isempty(x) || exist(x)==7);
p.addParameter('axpos',[],@(x) isempty(x) || (isnumeric(x) && length(x)==4));
p.addParameter('needcbar',1,@(x) x==0 || x==1);
p.addParameter('axlbl',1,@(x) x==0 || x==1);
p.parse(evtnums,solntype,varargin{:});
% 
parameters=p.Results;
solmats=parameters.solmats;
inputfile=parameters.inputfile;
procinfo=parameters.procinfo;
snrthresh=parameters.snrthresh;
noisewin=parameters.noisewin;
rho=parameters.rho;
mapzoom=parameters.mapzoom;
clrscheme=parameters.clrscheme;
mechdiam=parameters.mechdiam;
evtlbl=parameters.evtlbl;
makefig=parameters.makefig;
savedir=parameters.savedir;
axpos=parameters.axpos;
needcbar=parameters.needcbar;
axlbl=parameters.axlbl;
minmaxfreq=procinfo{4};

% Adjust inputs if necessary
if makefig==0
  savedir='';
end

% Load velocity model of the Tajik Basin for the inversion
velmdl=loadvelmdl();

% Plot the map of the Tajik Basin 
maprange=2;
addfaults=2;
addbrdr=0;
[mapfig,mapax,~]=plottopo_tajik(maprange,addfaults,addbrdr,makefig,'');
mapax.Units='normalized';
if ~isempty(axpos)
  mapax.Position(1)=axpos(1);
  mapax.Position(2)=axpos(2);
  mapax.Position(3)=axpos(3);
  mapax.Position(4)=axpos(4);
end
hold on

% Zoom into map?
if mapzoom>0
  evtnum=evtnums(1);
  evtinfo=getevtinfo(evtnum);
  minlon=evtinfo.LON_dec;
  maxlon=evtinfo.LON_dec;
  minlat=evtinfo.LAT_dec;
  maxlat=evtinfo.LAT_dec;
  if length(evtnums)>1
    for v=2:length(evtnums)
      evtnum=evtnums(v);
      evtinfo=getevtinfo(evtnum);
      nowlon=evtinfo.LON_dec;
      nowlat=evtinfo.LAT_dec;
      % 
      if nowlon<minlon
        minlon=nowlon;
      end
      if nowlon>maxlon
        maxlon=nowlon;
      end
      if nowlat<minlat
        minlat=nowlat;
      end
      if nowlat>maxlat
        maxlat=nowlat;
      end
    end
  end
  minlon=minlon-mapzoom;
  maxlon=maxlon+mapzoom;
  minlat=minlat-mapzoom;
  maxlat=maxlat+mapzoom;
  % Adjust Axes
  mapax.XLim=[minlon maxlon];
  mapax.YLim=[minlat maxlat];  
  % Adjust Axes labels
  if abs(maxlon-minlon)<0.02
    minxtick=min(unique(floor(mapax.XTick)));
    maxxtick=max(unique(ceil(mapax.XTick)));
    mapax.XTick=minxtick:0.0025:maxxtick;
  elseif abs(maxlon-minlon)<0.05
    minxtick=min(unique(floor(mapax.XTick)));
    maxxtick=max(unique(ceil(mapax.XTick)));
    mapax.XTick=minxtick:0.01:maxxtick;
  elseif abs(maxlon-minlon)<0.1
    minxtick=min(unique(floor(mapax.XTick)));
    maxxtick=max(unique(ceil(mapax.XTick)));
    mapax.XTick=minxtick:0.02:maxxtick;
  elseif abs(maxlon-minlon)<0.5
    minxtick=min(unique(floor(mapax.XTick)));
    maxxtick=max(unique(ceil(mapax.XTick)));
    mapax.XTick=minxtick:0.05:maxxtick;
  elseif abs(maxlon-minlon)<2
    minxtick=min(unique(floor(mapax.XTick)));
    maxxtick=max(unique(ceil(mapax.XTick)));
    mapax.XTick=minxtick:0.25:maxxtick;
  else
    minxtick=min(unique(floor(mapax.XTick)));
    maxxtick=max(unique(ceil(mapax.XTick)));
    mapax.XTick=minxtick:maxxtick;
  end
  % 
  if abs(maxlat-minlat)<0.02
    minytick=min(unique(floor(mapax.YTick)));
    maxytick=max(unique(ceil(mapax.YTick)));
    mapax.YTick=minytick:0.0025:maxytick;
  elseif abs(maxlat-minlat)<0.05
    minytick=min(unique(floor(mapax.YTick)));
    maxytick=max(unique(ceil(mapax.YTick)));
    mapax.YTick=minytick:0.01:maxytick;
  elseif abs(maxlat-minlat)<0.1
    minytick=min(unique(floor(mapax.YTick)));
    maxytick=max(unique(ceil(mapax.YTick)));
    mapax.YTick=minytick:0.02:maxytick;
  elseif abs(maxlat-minlat)<0.5
    minytick=min(unique(floor(mapax.YTick)));
    maxytick=max(unique(ceil(mapax.YTick)));
    mapax.YTick=minytick:0.05:maxytick;
  elseif abs(maxlat-minlat)<2
    minytick=min(unique(floor(mapax.YTick)));
    maxytick=max(unique(ceil(mapax.YTick)));
    mapax.YTick=minytick:0.25:maxytick;
  else
    minytick=min(unique(floor(mapax.YTick)));
    maxytick=max(unique(ceil(mapax.YTick)));
    mapax.YTick=minytick:maxytick;
  end
end

% Aspect ratio of map
mapaspect=mapax.PlotBoxAspectRatio;
mapaspect=mapaspect(2);
hold on
mapaspect=mapax.PlotBoxAspectRatio(2);
mapaspect=mapaspect*(max(mapax.XLim)-min(mapax.XLim))/(max(mapax.YLim)-min(mapax.YLim));


% Plot the moment tensor of each event onto the map
for v=1:length(evtnums)
  evtnum=evtnums(v);
  
  % Invert for the moment tensor if we do not already have a solution
  if isempty(solmats) 
    % Add path for inversion code
    codedir=addfocimtpath(0);
    if isempty(inputfile) 
      % Generate input file for inversion
      focimtdir=getfocimtdir(evtnum,'minmaxfreq',minmaxfreq);
      [inputfile,~]=focimt_1dvel_inpt(evtnum,rho,procinfo,snrthresh,...
        noisewin,focimtdir,[]);
    end
    
    % Run inversion
    [solmat,~,~]=focimt_1dvel(inputfile,velmdl,[],[]);
  else
    solmat=solmats{v};
    codedir='';
  end
  
  % Load solution
  load(solmat)
  % Full
  if solntype==1
    soln=Solution{1}.full;
      
  % Deviatoric
  elseif solntype==2
    soln=Solution{1}.deviatoric; 
      
  % Double-couple
  elseif solntype==3
    soln=Solution{1}.dc;
  end
  Mxx=soln.MXX(1);
  Mxy=soln.MXX(2);
  Mxz=soln.MXX(3);
  Myy=soln.MXX(4);
  Myz=soln.MXX(5);
  Mzz=soln.MXX(6);
  fm=[Mxx Myy Mzz Mxy Mxz Myz];
  
  % Load coordinates of event
  evtinfo=getevtinfo(evtnum);
  evtlon=evtinfo.LON_dec;
  evtlat=evtinfo.LAT_dec;
  
  % Plot focal mechanism on map
  % Don't label
  if evtlbl==0
    if solntype<3
      focalmech(fm,evtlon,evtlat,mechdiam,mapaspect,'XYZ') 
    elseif solntype==3
      focalmech(fm,evtlon,evtlat,mechdiam,mapaspect,'XYZ','DC')
    end
    
  % Label with event ID
  elseif evtlbl==1
    if solntype<3
      focalmech(fm,evtlon,evtlat,mechdiam,mapaspect,'XYZ','text',...
        sprintf('%d',evtnum),'FontSize',12)
    elseif solntype==3
      focalmech(fm,evtlon,evtlat,mechdiam,mapaspect,'XYZ','DC','text',...
        sprintf('%d',evtnum),'FontSize',12)
    end
  end
  hold on
end

% Recolor beach balls according to color-coding
numobj=length(mapax.Children);
for o=1:numobj
  nowobj=mapax.Children(o);
  if isa(nowobj,'matlab.graphics.primitive.Patch')
    % Compression
    if sum(nowobj.FaceColor)==0
      if clrscheme==1
        nowobj.FaceColor=[0.1 0.1 0.1]; 
      elseif clrscheme==2
        nowobj.FaceColor=[1 0.35 0.35];
      end
        
    % Dilatation   
    elseif sum(nowobj.FaceColor)==3
      if clrscheme==1
        nowobj.FaceColor=[1 1 1];
      elseif clrscheme==2
        nowobj.FaceColor=[0.25 0.35 1];
      end
    end
    % Transparency
    nowobj.FaceAlpha=0.85;
  end
end

% Adjust map title if necessary
if axlbl==1
  if length(evtnums)==1
    mapax.Title.String=sprintf('Event %d',evtnums);
  end
  
elseif axlbl==0
  mapax.Title.String='';
  mapax.XTick=[];
  mapax.XLabel.String='';
  mapax.YTick=[];
  mapax.YLabel.String='';
end


% Remove color bar if necessary
if needcbar==0
  for c=1:length(mapfig.Children)
    nowobj=mapfig.Children(c);
    if isa(nowobj,'matlab.graphics.illustration.ColorBar')
      delete(nowobj)
      break  
    end
  end
end

% Save Map?
if ~isempty(savedir) && exist(savedir)==7
  if length(evtnums)==1
    evtstr=sprintf('EVT%d_',evtnums);
  else
    evtstr=sprintf('EVTS%dto%d_',min(evtnums),max(evtnums));
  end
  % 
  if solntype==1
    solnstr='FULL';
  elseif solntype==2
    solnstr='DEV';
  elseif solntype==3
    solnstr='DC';
  end
  % 
  filterfreqs=procinfo{4};
  % 
  if mapzoom==0
    mapname=sprintf(...
      'FOCMECHMAP_%s%s_SNRTHR%.2f_NOISE%.2fs_%.2fto%.2fHz.eps',...
      evtstr,solnstr,snrthresh,noisewin,min(filterfreqs),...
      max(filterfreqs));
  else
    mapname=sprintf(...
      'FOCMECHMAP_%s%s_SNRTHR%.2f_NOISE%.2fs_%.2fto%.2fHz_MAPZOOM.eps',...
      evtstr,solnstr,snrthresh,noisewin,min(filterfreqs),...
      max(filterfreqs));
  end
  setenv('EPS',savedir);
  mapname=figdisp(mapname,[],[],2,[],'epstopdf');
  mapname=sprintf('%s.pdf',mapname(1:end-4));
else
  mapname='';
end


% Remove directories containing inversion code from path
if ~isempty(codedir)
  rmpath(genpath(codedir))
end



