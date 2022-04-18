function focax=mt2focplot(mtvec,clrscheme,varargin)
% 
% Function to plot a beach ball diagram from the components of a 
% moment tensor, in XYZ convention
% 
% INPUTS
% mtvec : Vector containing the 6 moment tensor components in this order:
% 
%         [Mxx Mxy Mxz Myy Myz Mzz]
% 
% 
% clrscheme : For plotting using focalmech.m, what color scheme do we use 
%             to mark stations with upwards P-wave polarity (corresponding 
%             to compressional quadrants in the focal mechanism) and 
%             downwards P-wave polarity (corresponding to dilatational 
%             quadrants in the focal mechanism)?
% 
%             1 : Dark grey for upwards polarity, white for downwards 
%                 polarity, grey for stations whose signals are too small 
%                 for a definite polarity [Default]
% 
%             2 : Red for upwards polarity, blue for downwards polarity
% 
%             3 : Blue for upwards polarity, red for downwards polarity
% 
% Input the following as varargin
% 
% coord : 'XYZ' or 'RTF'? Default is 'XYZ'
% isdc : 0 for not double couple, 1 for double couple
% 
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
% 
% Last Modified : April 4, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse through inputs
p=inputParser;
p.addRequired('mtvec',@(x) isnumeric(x));
p.addRequired('clrscheme',@(x) isnumeric(x));
p.addParameter('coord','XYZ',@(x) ischar(x));
p.addParameter('isdc',0,@(x) isnumeric(x));
p.parse(mtvec,clrscheme,varargin{:});
% 
parameters=p.Results;
coord=parameters.coord;
isdc=parameters.isdc;

% Get vector components
m11=mtvec(1);
m12=mtvec(2);
m13=mtvec(3);
m22=mtvec(4);
m23=mtvec(5);
m33=mtvec(6);

% Initialize axes
focax=axes;
pbaspect([1 1 1]);

% Plot moment tensor using focalmech.m
fm=[m11, m22, m33, m12, m13, m23];


% Plot using focalmech
if isdc==0
  if strcmp(coord,'XYZ')
    focalmech(fm,0,0,4,'XYZ');
  else
    focalmech(fm,0,0,4);
  end
else
  if strcmp(coord,'XYZ')
    focalmech(fm,0,0,4,'XYZ','DC');
  else
    focalmech(fm,0,0,4,'DC');
  end
end
hold on;

% Adjust colors of quadrants
numobj=length(focax.Children);
for o=1:numobj
  nowobj=focax.Children(o);
  if isa(nowobj,'matlab.graphics.primitive.Patch')
    % Compression
    if sum(nowobj.FaceColor)==0
      if clrscheme==1
        nowobj.FaceColor=[0.3 0.3 0.3]; 
      elseif clrscheme==2
        nowobj.FaceColor=[1 0.5 0.5];
      elseif clrscheme==3
        nowobj.FaceColor=[0.7 0.8 1];  
      end

    % Dilatation
    elseif sum(nowobj.FaceColor)==3
      if clrscheme==1
        nowobj.FaceColor=[1 1 1];
      elseif clrscheme==2
        nowobj.FaceColor=[0.7 0.8 1];
      elseif clrscheme==3
        nowobj.FaceColor=[1 0.5 0.5];  
      end
    end
  end
end

% Axis bounds
maxlim=1.02*0.5*4;
focax.XLim=[-1*maxlim maxlim];
focax.YLim=[-1*maxlim maxlim];

% Remove axes ticks
focax.XTick=[];
focax.YTick=[];
focax.Visible='off';

% 
evtmarker=scatter(0,0);
evtmarker.Marker='+';
if clrscheme>1
  evtmarker.MarkerEdgeColor=[0 0 0];
else
  evtmarker.MarkerEdgeColor=[1 0 0];  
end
evtmarker.SizeData=80;
evtmarker.LineWidth=2;



