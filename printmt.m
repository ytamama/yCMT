function [mtfig,mtax,mtname]=printmt(solmat,inptmat,solntype,varargin)
% 
% Function to "print out" the moment tensor solution of a seismic event, 
% as solved by fociMT
% 
% INPUTS
% solmat : Full path to the .mat file, containing the moment tensor
%          solution by fociMT
% 
% inptmat : Full path to the .mat file, containing the input parameters
%           processed by fociMT
% 
% solntype : Which moment tensor solution do we plot? 
%            1 : Full moment tensor inversion
%            2 : Deviatoric moment tensor inversion
%            3 : Double couple moment tensor inversion
% 
% Input the following as varargin : 
% 
% 'makefig' : Do we make a new figure to plot our distribution, or do we 
%             plot on an existing figure? 
%             0 : Plot on existing figure
%             1 : Plot on new figure [Default]
% 
% 'savedir' : In which directory do we save our inversion results and 
%             our radiation pattern? 
% 
%             Default: the same directory as the one containing inptmat
% 
%             If makefig = 0, don't save this figure
% 
% 'fontsize' : Font size to use for our values. Default: 10.5
% 
% 
% OUTPUT
% mtfig : Figure handle
% mtax : Axis handle
% mtname : Name of saved figure. '' if figure was not saved. 
% 
% References
% Kwiatek, G., Martínez‐Garzón, P. & Bohnhoff, M. (2016) HybridMT: A 
% MATLAB/Shell Environment Package for Seismic Moment Tensor Inversion and 
% Refinement. Seismological Research Letters, 87, 964–976. 
% doi:10.1785/0220150251
% 
% Hanks, T.C. & Kanamori, H. (1979) A moment magnitude scale. 
% Journal of Geophysical Research: Solid Earth, 84, 2348?2350. 
% doi:10.1029/JB084iB05p02348
% 
% Uses figdisp.m in csdms-contrib/slepian_alpha
% 
% Last Modified : March 20, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse through inputs
p=inputParser;
p.addRequired('solmat',@(x) exist(x)==2);
p.addRequired('inptmat',@(x) exist(x)==2);
p.addRequired('solntype',@(x) x==1 || x==2 || x==3);
p.addParameter('makefig',1,@(x) x==0 || x==1);
p.addParameter('savedir','',@(x) isempty(x) || exist(x)==7);
p.addParameter('fontsize',10.5,@(x) isnumeric(x));
p.parse(solmat,inptmat,solntype,varargin{:});
% 
parameters=p.Results;
makefig=parameters.makefig;
savedir=parameters.savedir;
fontsize=parameters.fontsize;

% Load solution
load(solmat);
% Full
if solntype==1
  solution=Solution{1}.full;

% Deviatoric
elseif solntype==2
  solution=Solution{1}.deviatoric;

% Double couple
elseif solntype==3
  solution=Solution{1}.dc;
end

% Moment tensor
mtensor=solution.MXX;
mxx=mtensor(1);
mxy=mtensor(2);
mxz=mtensor(3);
myy=mtensor(4);
myz=mtensor(5);
mzz=mtensor(6);
% Scalar seismic moment
m0=solution.M0;
% Total seismic moment: square-rooted squared sum of eigenvalues
mt=solution.MT;
% Moment magnitude, from Hanks + Kanamori (1979)
mw=(2/3)*log10(m0)-6.07;


% Calculate errors
[rmsu,rmso]=getMTerrs(solmat,inptmat,solntype);

% Fault plane solution, in degrees
f1=solution.F1;
strike=f1(1);
dip=f1(2);
rake=f1(3);

% Isotropic, CLVD, and double couple component
isoprc=solution.ISO;
clvdprc=solution.CLVD;
dcprc=solution.DC;


% Initialize figure
if makefig==1
  mtfig=figure();
elseif makefig==0
  mtfig=gcf;
end
mtax=axes;
mtax.YAxis.Visible='off';
mtax.XAxis.Visible='off';
mtax.Units='normalized';
mtax.Position(4)=0.75;

% Moment Tensor
mxxtxt1=text(0.05,0.86,'Mxx','FontSize',fontsize,'FontName','fixedwidth');
myxtxt1=text(0.05,0.76,'Myx','FontSize',fontsize,'FontName','fixedwidth');
mzxtxt1=text(0.05,0.66,'Mzx','FontSize',fontsize,'FontName','fixedwidth');
mxytxt1=text(0.15,0.86,'Mxy','FontSize',fontsize,'FontName','fixedwidth');
myytxt1=text(0.15,0.76,'Myy','FontSize',fontsize,'FontName','fixedwidth');
mzytxt1=text(0.15,0.66,'Mzy','FontSize',fontsize,'FontName','fixedwidth');
mxztxt1=text(0.25,0.86,'Mxz','FontSize',fontsize,'FontName','fixedwidth');
myztxt1=text(0.25,0.76,'Myz','FontSize',fontsize,'FontName','fixedwidth');
mzztxt1=text(0.25,0.66,'Mzz','FontSize',fontsize,'FontName','fixedwidth');
eqltxt=text(0.33,0.76,'=','FontSize',fontsize,'FontName','fixedwidth');
% 
if mxx>=0
  mxxtxt2=text(0.38,0.86,sprintf(' %.2e',mxx),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  mxxtxt2=text(0.38,0.86,sprintf('%.2e',mxx),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
if mxy>=0
  myxtxt2=text(0.38,0.76,sprintf(' %5.2e',mxy),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  myxtxt2=text(0.38,0.76,sprintf('%5.2e',mxy),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
if mxz>=0
  mzxtxt2=text(0.38,0.66,sprintf(' %5.2e',mxz),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  mzxtxt2=text(0.38,0.66,sprintf('%5.2e',mxz),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
% 
if mxy>=0
  mxytxt2=text(0.58,0.86,sprintf(' %.2e',mxy),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  mxytxt2=text(0.58,0.86,sprintf('%.2e',mxy),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
if myy>=0
  myytxt2=text(0.58,0.76,sprintf(' %.2e',myy),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  myytxt2=text(0.58,0.76,sprintf('%.2e',myy),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
if myz>=0
  mzytxt2=text(0.58,0.66,sprintf(' %.2e',myz),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  mzytxt2=text(0.58,0.66,sprintf('%.2e',myz),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
% 
if mxz>=0
  mxztxt2=text(0.78,0.86,sprintf(' %.2e',mxz),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  mxztxt2=text(0.78,0.86,sprintf('%.2e',mxz),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
if myz>=0
  mzytxt2=text(0.78,0.76,sprintf(' %.2e',myz),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  mzytxt2=text(0.78,0.76,sprintf('%.2e',myz),'FontSize',fontsize,...
    'FontName','fixedwidth');
end
if mzz>=0
  mzztxt2=text(0.78,0.66,sprintf(' %.2e',mzz),'FontSize',fontsize,...
    'FontName','fixedwidth');
else
  mzztxt2=text(0.78,0.66,sprintf('%.2e',mzz),'FontSize',fontsize,...
    'FontName','fixedwidth');
end

% Scalar moment + Error
m0txt=text(0.05,0.52,sprintf('M0: %.2e',m0),'FontSize',fontsize,...
  'FontName','fixedwidth');
mttxt=text(0.36,0.52,sprintf('MT: %.2e',mt),'FontSize',fontsize,...
  'FontName','fixedwidth');
mwtxt=text(0.67,0.52,sprintf('MW: %.2e',mw),'FontSize',fontsize,...
  'FontName','fixedwidth');
rmsutxt=text(0.05,0.42,sprintf('RMSU: %.2f%%',rmsu),'FontSize',fontsize,...
  'FontName','fixedwidth');
rmsotxt=text(0.36,0.42,sprintf('RMSO: %.2f%%',rmso),'FontSize',fontsize,...
  'FontName','fixedwidth');


% Strike, dip, and rake
strtxt1=text(0.05,0.28,sprintf('\\phi_s'),'FontSize',fontsize+3,...
  'FontName','fixedwidth','Color',[0 0 0]);
strtxt2=text(0.1,0.28,sprintf(': %.2f%s',strike,char(176)),...
  'FontSize',fontsize,'FontName','fixedwidth','Color',[0 0 0]);
diptxt1=text(0.36,0.28,sprintf('\\delta'),'FontSize',fontsize+3.5,...
  'FontName','fixedwidth','Color',[0 0 0]);
diptxt2=text(0.41,0.28,sprintf(': %.2f%s',dip,char(176)),...
    'FontSize',fontsize,'FontName','fixedwidth','Color',[0 0 0]);
raketxt1=text(0.67,0.28,sprintf('\\lambda'),...
    'FontSize',fontsize+3.5,'FontName','fixedwidth','Color',[0 0 0]);
raketxt2=text(0.72,0.28,sprintf(': %.2f%s',rake,char(176)),...
    'FontSize',fontsize,'FontName','fixedwidth','Color',[0 0 0]);

% Fraction of isotropic, CLVD, and double couple component
isotxt=text(0.05,0.14,sprintf('ISO: %.2f%%',isoprc),'FontSize',fontsize,...
  'FontName','fixedwidth','Color',[0 0 0]);
clvdtxt=text(0.36,0.14,sprintf('CLVD: %.2f%%',clvdprc),'FontSize',...
  fontsize,'FontName','fixedwidth','Color',[0 0 0]);
dctxt=text(0.67,0.14,sprintf('DC: %.2f%%',dcprc),'FontSize',fontsize,...
  'FontName','fixedwidth','Color',[0 0 0]);


%%%%%%%%%%%%
% Save figure
if makefig==1
  evtid=Solution{1}.event_id;
  mtname=sprintf('MOMENTTENSORSOLN_EVT%s.eps',evtid);
  if isempty(savedir)
    [savedir,~,~]=fileparts(solmat);
  end
  setenv('EPS',savedir);
  mtname=figdisp(mtname,[],[],2,[],'epstopdf');
    
elseif makefig==0
  mtname='';
end




