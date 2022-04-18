function [focfig,focax,figname]=plot_focimt(evtnum,resdir,solntype,...
  jkpar,issynth,varargin)
% 
% Function to plot the focal mechanism(s) generated by fociMT, for the 
% seismic events recorded by SeismoTech Ltd. 
% 
% INPUTS
% evtnum : ID number of the event, whose moment tensor we inverted for 
%          using fociMT.
% 
% resdir : Directory containing the results of the inversions conducted 
%          by fociMT. We will save our plots in this directory.
% 
% solntype : Which moment tensor solution do we plot? 
%            1 : Full moment tensor inversion
%            2 : Deviatoric moment tensor inversion
%            3 : Double couple moment tensor inversion
%            [] : Results from all three inversions. 
% 
% jkpar : With regards to jacknifing... 
%         [] : NOT part of a jacknife test
%         0 : Jacknife test, as encoded in fociMT
%         ID # of station : Part of a manual jacknife test, where we 
%                           exclude one station at a time. Input the ID 
%                           number of said station. 
%
% issynth : Is this focal mechanism generated by a synthetic inversion? 
%           0 : No
%           1 : Yes
% 
% 
% Input the following as varargin: 
% 
% 'clrscheme' : For plotting using focalmech.m, what color scheme do we use 
%               to mark stations with upwards P-wave polarity (corresponding 
%               to compressional quadrants in the focal mechanism) and 
%               downwards P-wave polarity (corresponding to dilatational 
%               quadrants in the focal mechanism)?
% 
%               1 : Dark grey for upwards polarity, white for downwards 
%                   polarity, grey for stations whose signals are too small 
%                   for a definite polarity [Default]
% 
%               2 : Red for upwards polarity, blue for downwards polarity
% 
%               3 : Blue for upwards polarity, red for downwards polarity
% 
%               [] : Default, use when plotting the figures made by 
%                    fociMT
% 
% 'makefig' : Do we make a new figure to plot our distribution, or do we 
%             plot on an existing figure? 
%             0 : Plot on existing figure
%             1 : Plot on new figure [Default]
% 
% 'savefig' : Do we save our figure? 
%             0 : No [Default]
%             1 : Yes, save in resdir
%   
%
% OUTPUTS
% focfig : Figure handle to beach ball diagram
% focax : Axis handle to beach ball diagram
% figname : Full path to the figure, saved in resdir. '' if not saved. 
% 
% References: 
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
% Uses figdisp.m and defval.m in csdms-contrib/slepian_alpha
% 
% Last Modified : March 27, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse through inputs
p=inputParser;
p.addRequired('evtnum',@(x) isnumeric(x));
p.addRequired('resdir',@(x) exist(x)==7);
p.addRequired('solntype',@(x) isnumeric(x));
p.addRequired('jkpar',@(x) isempty(x) || isnumeric(x));
p.addRequired('issynth',@(x) isnumeric(x));
p.addParameter('clrscheme',{},@(x) isempty(x) || isnumeric(x));
p.addParameter('makefig',1,@(x) isnumeric(x));
p.addParameter('savefig',0,@(x) isnumeric(x));
p.parse(evtnum,resdir,solntype,jkpar,issynth,varargin{:});
% 
parameters=p.Results;
clrscheme=parameters.clrscheme;
makefig=parameters.makefig;
savefig=parameters.savefig;


% Initialize figure
if makefig==1
  focfig=figure();
  focfig.Units='normalized';
  if isempty(solntype)
    focfig.Position(1)=0.1;
    focfig.Position(2)=0.1;
    focfig.Position(3)=0.75;
    focfig.Position(4)=0.5;
  else
    focfig.Position(2)=0.35;
    focfig.Position(3)=0.35;
    focfig.Position(4)=0.5;
  end
elseif makefig==0
  focfig=gcf;
end


% Plot using focalmech.m
% FULL
if isempty(solntype) || solntype==1
  % Load full moment tensor solution
  solmat=fullfile(resdir,'Solution.mat');
  load(solmat);
  fullsoln=Solution{1}.full;
  fullmt=fullsoln.MXX;
  m11=fullmt(1);
  m12=fullmt(2);
  m13=fullmt(3);
  m22=fullmt(4);
  m23=fullmt(5);
  m33=fullmt(6);
  mw=fullsoln.MW;
  diam=4;

  % Plot moment tensor using focalmech.m
  fm=[m11, m22, m33, m12, m13, m23];
  % Skip if it doesn't work
  try
    focax=mt2focplot(fm,clrscheme);  
    hold on
  catch
    keyboard
  end
  
  % Axis bounds
  maxlim=1.02*0.5*diam;
  focax.XLim=[-1*maxlim maxlim];
  focax.YLim=[-1*maxlim maxlim];
  
  % Remove axes ticks
  focax.XTick=[];
  focax.YTick=[];
  
  % Label Plot
  % One of three
  if isempty(solntype)
    focax.XLabel.String='Full';
    focax.FontSize=12;
    focax.Position(1)=0.04;
    focax.Position(3)=0.28;
    focax.XLim=[-1*maxlim maxlim];
    focax.YLim=[-1*maxlim maxlim];
    pbaspect([1 1 1]);
      
  % Only plot
  elseif solntype==1
    if issynth==0
      titlestr1='Full Moment Tensor';
    else
      titlestr1='Synthetic Full Moment Tensor';
    end
    titlestr2=sprintf('Solution for Event %d',evtnum);
    if isempty(jkpar)
      focax.Title.String={titlestr1;titlestr2};
    elseif jkpar==0
      focax.Title.String={titlestr1;titlestr2;'Jacknife Test'};
    else
      focax.Title.String={titlestr1;titlestr2;...
        sprintf('Jacknife Test (sans Station %d)',jkpar)};
    end
  end
  
  % Plot a cross where our event would be
  evtmarker=scatter(0,0);
  evtmarker.Marker='+';
  evtmarker.MarkerEdgeColor=[0 0 0];
  evtmarker.SizeData=300;
  evtmarker.LineWidth=4;
end

% DEVIATORIC
if isempty(solntype) || solntype==2
  % Load deviatoric moment tensor solution
  solmat=fullfile(resdir,'Solution.mat');
  load(solmat);    
  devsoln=Solution{1}.deviatoric;
  devmt=devsoln.MXX;
  m11=devmt(1);
  m12=devmt(2);
  m13=devmt(3);
  m22=devmt(4);
  m23=devmt(5);
  m33=devmt(6);
  mw=devsoln.MW;
  diam=4;
    
  % Initialize axes
  focax=axes;
  pbaspect([1 1 1]);
  
  % Plot moment tensor using focalmech.m
  fm=[m11, m22, m33, m12, m13, m23];
  % Skip if it doesn't work
  try
    focax=mt2focplot(fm,clrscheme);  
    hold on
  catch
    keyboard
  end
  
  % Axis bounds
  maxlim=1.02*0.5*diam;
  focax.XLim=[-1*maxlim maxlim];
  focax.YLim=[-1*maxlim maxlim];
  
  % Remove axes ticks
  focax.XTick=[];
  focax.YTick=[];
  
  % Label Plot
  % One of three plots
  if isempty(solntype)
    focax.Position(1)=0.358;
    focax.Position(3)=0.28;
    focax.XLabel.String='Deviatoric';
    focax.FontSize=12;
    focax.XLim=[-1*maxlim maxlim];
    focax.YLim=[-1*maxlim maxlim];
    pbaspect([1 1 1]);
    if issynth==0
      titlestr1='Deviatoric Moment Tensor';
    else
      titlestr1='Synthetic Deviatoric Moment Tensor';
    end
      
  % Only plot on figure
  else
    if issynth==0
      titlestr1='Deviatoric Moment Tensor';
    else
      titlestr1='Synthetic Deviatoric Moment Tensor';
    end
  end
  titlestr2=sprintf('Solution for Event %d',evtnum);
  % 
  if isempty(jkpar)
    focax.Title.String={titlestr1;titlestr2};
  elseif jkpar==0
    focax.Title.String={titlestr1;titlestr2;'Jacknife Test'};
  else
    focax.Title.String={titlestr1;...
      sprintf('Jacknife Test (sans Station %d)',jkpar)};
  end
  if isempty(solntype)
    focax.Title.Units='normalized';
    focax.Title.Position(2)=1.03;
  end
  
  % Plot a cross where our event would be
  evtmarker=scatter(0,0);
  evtmarker.Marker='+';
  evtmarker.MarkerEdgeColor=[0 0 0];
  evtmarker.SizeData=300;
  evtmarker.LineWidth=4;
end


% DOUBLE COUPLE
if isempty(solntype) || solntype==3
  % Load double couple moment tensor solution
  solmat=fullfile(resdir,'Solution.mat');
  load(solmat);
  dcsoln=Solution{1}.dc;
  dcmt=dcsoln.MXX;
  m11=dcmt(1);
  m12=dcmt(2);
  m13=dcmt(3);
  m22=dcmt(4);
  m23=dcmt(5);
  m33=dcmt(6);
  mw=dcsoln.MW;
  diam=4;
    
  % Initialize axes
  focax=axes;
  pbaspect([1 1 1]);

  % Plot moment tensor using focalmech.m
  fm=[m11, m22, m33, m12, m13, m23];
  % Skip if it doesn't work
  try
    focax=mt2focplot(fm,clrscheme);  
    hold on
  catch
    keyboard
  end
  
  % Axis bounds
  maxlim=1.02*0.5*diam;
  focax.XLim=[-1*maxlim maxlim];
  focax.YLim=[-1*maxlim maxlim];
  
  % Remove axes ticks
  focax.XTick=[];
  focax.YTick=[];
  
  % Label Plot
  % One of three
  if isempty(solntype)
    focax.Position(1)=0.675;
    focax.Position(3)=0.28;
    focax.XLabel.String='Double Couple';
    focax.FontSize=12;
    focax.XLim=[-1*maxlim maxlim];
    focax.YLim=[-1*maxlim maxlim];
    pbaspect([1 1 1]);
      
  % Only plot
  elseif solntype==3
    if issynth==0
      titlestr1='Double Couple Moment Tensor';
    else
      titlestr1='Synthetic Double Couple Moment Tensor'; 
    end
    titlestr2=sprintf('Solution for Event %d',evtnum);
    if isempty(jkpar)
      focax.Title.String={titlestr1;titlestr2};
    elseif jkpar==0
      focax.Title.String={titlestr1;titlestr2;'Jacknife Test'};
    else
      focax.Title.String={titlestr1;titlestr2;...
        sprintf('Jacknife Test (sans Station %d)',jkpar)};
    end
  end
  
  % Plot a cross where our event would be
  evtmarker=scatter(0,0);
  evtmarker.Marker='+';
  evtmarker.MarkerEdgeColor=[0 0 0];
  evtmarker.SizeData=300;
  evtmarker.LineWidth=4;
end


% Save figure
if savefig==1
  if isempty(jkpar)
    jkstr='';
  elseif jkpar==0
    jkstr='_JACKNIFE';
  else
    jkstr=sprintf('_JACKNIFE.NOSTA%d',jkpar);
  end
  %
  if issynth==0
    synthstr='';
  else
    synthstr='_SYNTH';
  end 
  %
  if plottype==1
    plotstr='';
  elseif plottype==2
    plotstr='fociMTPLOT_';
  end
  %
  if isempty(solntype)
    figname=sprintf('EVT%d_%sfociMTSolutions%s%s.eps',evtnum,plotstr,...
      synthstr,jkstr);
  elseif solntype==1
    figname=sprintf('EVT%d_%sfociMTFULLSOLN%s%s.eps',evtnum,plotstr,...
      synthstr,jkstr);
  elseif solntype==2
    figname=sprintf('EVT%d_%sfociMTDEVSOLN%s%s.eps',evtnum,plotstr,...
      synthstr,jkstr);
  elseif solntype==3
    figname=sprintf('EVT%d_%sfociMTDCSOLN%s%s.eps',evtnum,plotstr,...
      synthstr,jkstr);
  end
  setenv('EPS',resdir);
  figname=figdisp(figname,[],[],2,[],'epstopdf');
else
  figname='';
end



