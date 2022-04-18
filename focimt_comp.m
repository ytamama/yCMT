function [compfigs,compnames,originpt,solmat,paramat,inptmat]=...
  focimt_comp(evtnum,mttype,solntype,varargin)
% 
% Function to either invert for the moment tensor of, or generate 
% synthetic, characteristic moment tensors from the distribution of 
% stations about an earthquake, recorded in the Tajik Basin by 
% SeismoTech Ltd. Note that we conduct all inversions using fociMT, and 
% the beach ball diagrams are plotted using focalmech.m. 
% 
% Note that all inversions are conducted with seismograms in the 
% displacement-time domain, deconvolved and filtered to 3 to 6 Hz. 
% Specifically, we use first-arriving P-waves with a signal-to-noise 
% ratio of at least 2, where the noise is defined as the seismic data
% 3s before the P-wave arrival determined by SeismoTech Ltd. 
% 
% INPUTS
% evtnum : ID Number of seismic event
% 
% mttype : Which type of inversion do we wish to conduct? 
%          0 : Invert for the moment tensor of our real data
%          1 : Strike-slip
%          2 : Reverse 
%          3 : Normal 
%          4 : Strike-slip (rotated 45 degrees)
%          5 : Vertical dip slip fault (Rotated)
%          6 : Vertical dip slip fault
%          7 : Reverse (Rotated)
% 
% solntype : Which moment tensor solution do we wish to feature? 
%            1 : Full
%            2 : Deviatoric
%            3 : Double-couple
%           [] : All three (only for mttype = 0)
% 
% Input the following as varargin: 
% 
% 'inputfile' : Input file to use in our inversion for fociMT, based on 
%               real data. If we don't have this file, we would create it 
%               during this routine. 
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
% 
% 'threshval' : To find the peak corresponding to the P-wave arrival,
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
% 'savedir' : Full path to the directory to save our inversion,
%             inputted as a string. 
% 
% 'evte' : The easting of the event in m, which we wish to use in this 
%          synthetic inversion instead of the actual easting of the event.
% 
% 'evtn' : The northing of the event in m, which we wish to use in this 
%          synthetic inversion instead of the actual northing of the event.
% 
% 'evtz' : The depth of the event in m, which we wish to use in this 
%          synthetic inversion instead of the actual depth of the event
%          
% 'staes' : The easting coordinates of each station, in m, which we wish
%           to use in this synthetic inversion instead of the actual 
%           coordinates. If we wish to use the same coordinate for all 
%           events, input just one number 
% 
% 'stans' : The northing coordinates of each station, in m, which we wish
%           to use in this synthetic inversion instead of the actual 
%           coordinates. If we wish to use the same coordinate for all 
%           events, input just one number
% 
% 'stazs' : The depth coordinates of each station, in m, which we wish
%           to use in this synthetic inversion instead of the actual 
%           coordinates. If we wish to use the same coordinate for all 
%           events, input just one number. Note that these coordinates 
%           are positive the higher the altitude. 
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
%               2 : Red for upwards polarity, blue for downwards polarity
% 
%               3 : Blue for upwards polarity, red for downwards polarity
% 
% 
% References
% Characteristic moment tensors are from Dahlen & Tromp (1998), with the 
% coordinates converted into the Cartesian System using Aki & Richards 
% (2002). 
% 
% Last Modified : March 16, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For stations at 0 m depth, having an event depth of +/- 0 to 0.5 
% doesn't work for some reason. The magnitude of the depth has to be 
% at least 1 m

% Parse through inputs
p=inputParser;
p.addRequired('evtnum',@(x) isnumeric(x));
p.addRequired('mttype',@(x) isnumeric(x));
p.addRequired('solntype',@(x) isempty(x) || isnumeric(x));
p.addParameter('inputfile',{},@(x) isempty(x) || exist(x)==2);
p.addParameter('procinfo',{},@(x) isempty(x) ||... 
  (iscell(x) && length(x)==4))
p.addParameter('threshval',2,@(x) isnumeric(x));
p.addParameter('noisewin',3,@(x) isnumeric(x));
p.addParameter('rho',2500,@(x) isnumeric(x));
p.addParameter('savedir',{},@(x) isempty(x) || exist(x)==7);
p.addParameter('evte',{},@(x) isempty(x) || isnumeric(x));
p.addParameter('evtn',{},@(x) isempty(x) || isnumeric(x));
p.addParameter('evtz',{},@(x) isempty(x) || isnumeric(x));
p.addParameter('staes',{},@(x) isempty(x) || isnumeric(x));
p.addParameter('stans',{},@(x) isempty(x) || isnumeric(x));
p.addParameter('stazs',{},@(x) isempty(x) || isnumeric(x));
p.addParameter('clrscheme',1,@(x) isnumeric(x));
p.parse(evtnum,mttype,solntype,varargin{:});
% 
parameters=p.Results;
inputfile=parameters.inputfile;
procinfo=parameters.procinfo;
threshval=parameters.threshval;
noisewin=parameters.noisewin;
rho=parameters.rho;
savedir=parameters.savedir;
evte=parameters.evte;
evtn=parameters.evtn;
evtz=parameters.evtz;
staes=parameters.staes;
stans=parameters.stans;
stazs=parameters.stazs;
clrscheme=parameters.clrscheme;


% Add paths containing code
codedir=addfocimtpath(mttype);

% Create parent directory to store results of inversion
if isempty(procinfo)
  minmaxfreq=[];
else
  minmaxfreq=procinfo{4};
end
focimtdir=getfocimtdir(evtnum,'minmaxfreq',minmaxfreq);

% Load velocity model of the Tajik Basin for the inversion
velmdl=loadvelmdl();

% If making NO changes to our data
if mttype==0 && (isempty(evte) && isempty(evtn) && isempty(evtz)...
   && isempty(staes) && isempty(stans) && isempty(stazs))
  issynth=0;
else
  issynth=1;
end

% Type of moment tensor to invert for / plot
if isempty(solntype) && issynth==1
  % Default for synthetic inversion: full moment tensor
  solntype=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct input file
if isempty(inputfile) || exist(inputfile)==0
  inputfile=focimt_1dvel_inpt(evtnum,rho,procinfo,threshval,noisewin,...
    focimtdir,[]);
end
originpt=inputfile;

% Construct new input file, if we wish to alter where the event and 
% stations are located
if ~(isempty(evte) && isempty(evtn) && isempty(evtz) && isempty(staes)...
   && isempty(stans) && isempty(stazs))

  oldinput=inputfile;
  [inputdir,~,~]=fileparts(oldinput);
  newinputfile=synth_focimt_inpt(oldinput,inputdir,'minmaxfreq',...
    minmaxfreq,'evte',evte,'evtn',evtn,'evtz',evtz,'staes',staes,...
    'stans',stans,'stazs',stazs);
  inputfile=newinputfile;
end

% Run fociMT.m, calculating the theoretical areas under the P-waves 
% based on the characteristic moment tensor
try
  focimt(inputfile,'VelocityModel',velmdl,'ProjectDir',focimtdir);
catch
  return
end
% Moment tensor solution
solmat=fullfile(focimtdir,'Solution.mat');
inptmat=fullfile(focimtdir,'Input.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run fociMT again if we wish to generate synthetic characteristic 
% moment tensors

% Prepare the directory to save our results, if necessary
if isempty(savedir)
  % Real Data (non-synthetic moment tensor)
  if mttype==0
    savedir=focimtdir;
      
  % Strike-slip
  elseif mttype==1
    savedir=fullfile(focimtdir,'STRIKESLIP');
  
  % Reverse 
  elseif mttype==2
    savedir=fullfile(focimtdir,'REVERSE');
  
  % Normal
  elseif mttype==3
    savedir=fullfile(focimtdir,'NORMAL');
    
  % Strike-slip (rotated 45 degrees)
  elseif mttype==4
    savedir=fullfile(focimtdir,'SSROT');
      
  % Vertical dip slip (rotated)
  elseif mttype==5
    savedir=fullfile(focimtdir,'VERTSLIPROT');
  
  % Vertical dip slip
  elseif mttype==6
    savedir=fullfile(focimtdir,'VERTSLIP');
      
  % Reverse (rotated)
  elseif mttype==7
    savedir=fullfile(focimtdir,'REVHORZ');
  end
  if exist(savedir)==0
    mkdir(savedir);
  end
end

% Run the synthetic inversion if necessary
if mttype>0
  paramat=fullfile(focimtdir,'Params.mat');
  [synthsol,synthin,~,synthinfile]= synth_focimt(inptmat,paramat,...
    savedir,'minmaxfreq',minmaxfreq,'calcorth',solntype,'solmat',...
    solmat,'evte',evte,'evtn',evtn,'evtz',evtz,'staes',staes,...
    'stans',stans,'stazs',stazs);
  inputfile=synthinfile;
end


% Final moment tensor solution
solmat=fullfile(savedir,'Solution.mat');
inptmat=fullfile(savedir,'Input.mat');
paramat=fullfile(savedir,'Params.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot a "dashboard figure" of radiation pattern, distribution of 
% stations, moment tensor solution, and beachball diagram

compfigs=[];
compnames={};
solntypes=[1 2 3];

for s=1:length(solntypes)
  nowtype=solntypes(s);
  if isempty(solntype) || nowtype==solntype

    % Initialize figure
    compfig=figure();
    compfig.Units='normalized';
    compfig.Position(1)=0.14;
    compfig.Position(2)=0.05;
    if issynth==1
      compfig.Position(3)=0.8;
    elseif issynth==0
      compfig.Position(3)=0.86;
    end
    compfig.Position(4)=0.85;
    figclr=get(compfig,'Color');

    % First, plot the radiation pattern in the top left 
    [~,projax,~]=stnproj(inputfile,nowtype,mttype,'velmdl',velmdl,...
      'makefig',0);
    hold on
    projax.Units='normalized';
    projax.Position(2)=0.27;
    if issynth==1
      projax.Position(1)=0.03;
      projax.Position(3)=0.31;
    elseif issynth==0
      projax.Position(1)=0.02;
      projax.Position(3)=0.28; 
    end
    projax.Title.String='';
    projax.XLabel.String='Takeoff Angle Radiation Pattern';
    projax.FontSize=11;
    projax.XColor=[1 1 1];
    projax.YColor=[1 1 1];
    projax.XLabel.Color=[0 0 0];

    % Next to it, plot the geographical distribution of stations
    [~,uax,~]=udistr(inptmat,issynth,0,clrscheme,0);
    hold on
    uax.Units='normalized';
    uax.Position(2)=0.27;
    if issynth==1
      uax.Position(1)=0.39;
      uax.Position(3)=0.31;
    elseif issynth==0
      uax.Position(1)=0.326;
      uax.Position(3)=0.28;  
    end
    uax.XColor=[0.15 0.15 0.15];
    uax.YColor=[0.15 0.15 0.15];
    uax.XTick=[];
    uax.YTick=[];
    uax.Title.String='';
    uax.XLabel.String='Geographical Radiation Pattern';
    uax.YLabel.String='';
    uax.XLabel.Color=[0 0 0];
    uax.FontSize=11;
    uax.Box='on';

    % Plot the moment tensor solution on the bottom left
    [~,mtax,~]=printmt(solmat,inptmat,nowtype,'makefig',0);
    hold on
    mtax.Units='normalized';
    if issynth==1
      mtax.Position(1)=0.03;
      mtax.Position(2)=0.07;
      mtax.Position(3)=0.31;
      mtax.Position(4)=0.28;
    elseif issynth==0
      mtax.Position(1)=0.02;
      mtax.Position(2)=0.08;
      mtax.Position(3)=0.28;
      mtax.Position(4)=0.26;
    end
    mtax.XAxis.Visible='on';
    mtax.XTick=[];
    mtax.YAxis.Visible='on';
    mtax.YTick=[];
    mtax.Box='on';
    mtax.XLabel.String='Moment Tensor Solution';
    mtax.FontSize=11;


    % Plot depth distribution of stations : East-West binning
    [~,depax_ew,~]=depdistr_focimt(inptmat,issynth,'numbins',10,...
      'makefig',0,'bindir',0,'rotplot',0);
    hold on
    depax_ew.Units='normalized';
    if issynth==1
      depax_ew.Position(1)=0.39;
      depax_ew.Position(2)=0.07;
      depax_ew.Position(3)=0.31;
      depax_ew.Position(4)=0.28;
    elseif issynth==0
      depax_ew.Position(1)=0.326;
      depax_ew.Position(2)=0.08;
      depax_ew.Position(3)=0.28;
      depax_ew.Position(4)=0.26;
    end
    depax_ew.Title.String='';
    depax_ew.FontSize=11;
    depax_ew.YAxisLocation='right';

    % Plot depth distribution of stations : North-South binning
    [~,depax_ns,~]=depdistr_focimt(inptmat,issynth,'numbins',10,...
      'makefig',0,'bindir',1,'rotplot',1);
    hold on
    depax_ns.Units='normalized';
    if issynth==1
      depax_ns.Position(1)=0.75;
      depax_ns.Position(2)=0.41;
      depax_ns.Position(3)=0.18;
      depax_ns.Position(4)=0.525;
    elseif issynth==0
      depax_ns.Position(1)=0.645;
      depax_ns.Position(2)=0.415;
      depax_ns.Position(3)=0.15;
      depax_ns.Position(4)=0.527;
    end
    depax_ns.Title.String='';
    depax_ns.YLabel.Rotation=90;
    depax_ns.YLabel.Position(1)=1.12;
    depax_ns.XAxisLocation='bottom';
    depax_ns.XTickLabelRotation=90;
    depax_ns.FontSize=11;
    depax_ns.XDir='reverse';
    depax_ns.YDir='normal';
    depax_ns.YTickLabelRotation=90;
    depax_ns.YAxisLocation='right';


    % Plot the actual beachball
    load(solmat);
    if nowtype==1
      nowsoln=Solution{1}.full;
    elseif nowtype==2
      nowsoln=Solution{1}.deviatoric;
    elseif nowtype==3
      nowsoln=Solution{1}.dc;   
    end
    mtvec=nowsoln.MXX;
    % Double couple
    if nowtype==3
      bbax=mt2focplot(mtvec,clrscheme,'coord','XYZ','isdc',1);  
    else
      bbax=mt2focplot(mtvec,clrscheme,'coord','XYZ','isdc',0);
    end   
    bbax.Title.String='';
    bbax.Units='normalized';
    bbax.XLabel.String='Beach ball';
    if issynth==1
      bbax.Position(1)=0.454;
      bbax.Position(2)=0.07;
      bbax.Position(4)=0.28;
    elseif issynth==0
      bbax.Position(1)=0.332;
      bbax.Position(2)=0.09;
      bbax.Position(4)=0.24;
    end
    bbax.XAxis.Color=figclr;
    bbax.YAxis.Color=figclr;
    bbax.XLabel.Color=[0 0 0];
    bbax.FontSize=11;
    
    
    % Plot the record-section of P-waves, by azimuth
    if issynth==0
      [~,recax,~]=precordsec(evtnum,procinfo,threshval,noisewin,...
        'mint',0.2,'maxt',0.4,'scaleopt',1,'scaleval',4.5,...
        'shownum',20,'seisorder',2,'makefig',0,'clrscheme',3,...
        'addhorz',1,'labelp',1);
      recax.Units='normalized';
      recax.YAxisLocation='right';
      recax.Box='on';
      recax.Position(1)=0.84;
      recax.Position(2)=0.415;
      recax.Position(3)=0.138;
      recax.Position(4)=0.527;
      recax.Title.String='';
      recax.XLabel.String='Time (s) wrt P-wave Onset';
      recax.FontSize=11;
      recax.XLabel.FontSize=11;
    
      % Plot beach ball on a map (zoomed out)
      mapax=axes;
      axpos=[0.84 0.08 0.138 0.26];
      [mapfig,mapax,mapname]=focimt_map(evtnum,nowtype,'solmat',...
        {solmat},'procinfo',procinfo,'mapzoom',0.1,'clrscheme',1,'mechdiam',...
        0.04,'makefig',0,'axpos',axpos);
      mapax.FontSize=11;
      mapax.Title.String='';
    
      % Delete the colorbar
      for o=1:length(compfig.Children)
        nowobj=compfig.Children(1);
        if strcmp(class(nowobj),'matlab.graphics.illustration.ColorBar')
          delete(nowobj)
        end
      end
    end

    % Figure Title! 
    if mttype>0
      faulttypes={'Strike Slip';'Reverse';'Normal';'Strike Slip (Rotated)';...
        'Vertical Dip Slip (Rotated)';'Vertical Dip Slip';'Reverse (Rotated)'};
      if nowtype==1
        uax.Title.String=sprintf('Full Moment Tensor Solution of a %s Fault',...
          faulttypes{mttype});
      elseif nowtype==2
        uax.Title.String=sprintf('Deviatoric Moment Tensor Solution of a %s Fault',...
          faulttypes{mttype});
      elseif nowtype==3
        uax.Title.String=sprintf(...
          'Constrained Double Couple Moment Tensor Solution of a %s Fault',...
          faulttypes{mttype});
      end

    elseif mttype==0
      if nowtype==1
        if isempty(minmaxfreq)
          uax.Title.String=sprintf('Full Moment Tensor Solution of Event %d',...
            evtnum);
        else
          uax.Title.String=sprintf(...
            'Full Moment Tensor Solution of Event %d (%.2g - %.2g Hz)',...
            evtnum,min(minmaxfreq),max(minmaxfreq));
        end
      elseif nowtype==2
        if isempty(minmaxfreq)
          uax.Title.String=sprintf(...
            'Deviatoric Moment Tensor Solution of Event %d',...
            evtnum);
        else
          uax.Title.String=sprintf(...
            'Deviatoric Moment Tensor Solution of Event %d (%.2g - %.2g Hz)',...
            evtnum,min(minmaxfreq),max(minmaxfreq));
        end
      elseif nowtype==3
        if isempty(minmaxfreq)
          uax.Title.String=sprintf(...
            'Constrained Double Couple Moment Tensor Solution of Event %d',...
            evtnum);
        else
          uax.Title.String=sprintf(...
            'Constrained Double Couple Moment Tensor Solution of Event %d (%.2g - %.2g Hz)',...
            evtnum,min(minmaxfreq),max(minmaxfreq));
        end
      end
    end
    uax.Title.Units='normalized';
    if issynth==1
      uax.Title.Position(1)=0.4;
      uax.Title.Position(2)=1.04;
    elseif issynth==0
      uax.Title.Position(1)=0.7;
      uax.Title.Position(2)=1.04;
    end
    uax.Title.FontSize=14;


    % Save figure as PDF and PNG
    if nowtype==1
      solnstr='FULL';
    elseif nowtype==2
      solnstr='DEV';
    elseif nowtype==3
      solnstr='DC';
    end
    % 
    if mttype==0
      mtstr='';
    elseif mttype==1
      mtstr='_SS';
    elseif mttype==2
      mtstr='_REV';
    elseif mttype==3
      mtstr='_NORMAL';
    elseif mttype==4
      mtstr='_SSROT';
    elseif mttype==5
      mtstr='_VERTSLIPROT';
    elseif mttype==6
      mtstr='_VERTSLIP';
    elseif mttype==7
      mtstr='_REVROT';
    end
    if issynth==0
      compname=sprintf('FOCIMTDASHBOARD_EVT%d%s_%s.png',evtnum,mtstr,solnstr);
    elseif issynth==1
      compname=sprintf('FOCIMTDASHBOARD_EVT%d%s_%s_SYNTH.png',evtnum,mtstr,...
        solnstr);  
    end
    pause(2);
    compname=fullfile(savedir,compname);
    saveas(compfig,compname)
    
    % 
    compfigs=vertcat(compfigs,compfig);
    compnames=vertcat(compnames,compname);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove paths
rmpath(genpath(codedir))



