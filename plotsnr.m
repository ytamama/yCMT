function snrtex=plotsnr(snrvals,snrdef,snrwins,arrtypes,axhdl,starttime,...
  definfo,snrclrs)
% 
% Function to plot representations of the signal to noise ratio (SNR) of
% one or more picked phase arrivals on the corresponding seismogram.
% 
% INPUTS
% snrvals : Vector of SNR values to add to a seismogram
% snrdef : How were the SNR values calculated? 
%          1 : Amplitude: max(abs(signal - mean(signal))) / max(abs(noise - mean(noise))) 
%          2 : Variance: variance(signal) / variance(noise-mean(noise))
%          3 : Root mean square: RMS(signal) / RMS(noise-mean(noise))
%          4 : Peak to noise: abs(chosen peak following arrival) /
%          mean(abs(peak height of the noise))
% 
% snrwins : Cell array, with each entry corresponding to the inputted
%           SNR values in snrvals. Each entry is a 3x1 vector, where: 
%
%           The first element is the starting time of the noise window, in 
%           UTC. 
%           The second element is (snrdef = 1,2,3) the ending time of the
%           signal window OR (snrdef = 4) time of the peak, 
%           following the picked arrival, chosen as the signal in SNR 
%           calculation. Both times are in UTC. 
%           The third element is the picked arrival time, in UTC (also 
%           the ending time of the noise window). 
% 
% arrtypes : For what arrivals did we calculate the SNR? Enter as a cell
%            array, with one entry per inputted arrival. 
%            Example: {'P';'S'};
% 
% axhdl : Present axis handle 
% 
% starttime : Time corresponding to the 'zero' in the seismogram.
%
% definfo : Outputs definfo from each SNR calculation in seis2snr.m, based 
%           on the inputted SNR calculation type snrdef: 
% 
%           1 : [max(abs(noise)), max(abs(signal))]
%           2 : [variance(noise), variance(signal)]
%           3 : [RMS(noise), RMS(signal)]
%           4 : [mean(abs(peak height of the noise)), 
%               height of signal peak selected for SNR calculation]
% 
% snrclrs : A matrix containing 1x3 vectors, representing the color needed
%           to plot each SNR value. 
% 
% 
% axhdl : Current axis handle
% starttime : Time corresponding to the very start of the plot, in UTC. 
% 
% OUTPUT
% snrtex : Text handle
% 
% Last Modified: December 29, 2021 by Yuri Tamama
% 
% Uses boxtex.m in csdms-contrib/slepian_alpha
% Uses defval.m in csdms-contrib/slepian_alpha
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bounds of plot axes
minx=min(axhdl.XLim);
maxx=max(axhdl.XLim);
miny=min(axhdl.YLim);
maxy=max(axhdl.YLim);


snrstrs={};
numvals=length(snrvals);
maxlen=0;
% Iterate through all inputted SNR values
for s=1:numvals
  % Construct string containing SNR value
  nowsnr=snrvals(s);
  nowarr=arrtypes{s};
  nowstr=sprintf('%s SNR: %.2f',nowarr,nowsnr);
  snrstrs=horzcat(snrstrs,nowstr);
  % Check length
  nowlen=length(nowstr);
  if nowlen>maxlen
    maxlen=nowlen;
  end
  
  % Plot color
  plotclr=snrclrs(1,:);
  
  % Plot rectangular box of noise
  snrwin=snrwins{s};
  noisex1=snrwin(1);
  noisex1=seconds(noisex1-starttime);
  noisex2=snrwin(3);
  noisex2=seconds(noisex2-starttime);
  noisey1=-1*definfo(1);
  noisey2=definfo(1);
  %
  noisex=[noisex1 noisex1 noisex2 noisex2];
  noisey=[noisey1 noisey2 noisey2 noisey1];
  noisebox=fill(noisex,noisey,plotclr);
  noisebox.EdgeAlpha=0;
  noisebox.FaceAlpha=0.05;
  noisebox.HandleVisibility='off';
  hold on
  
  
  % SNRDEF < 4 : Plot box around signal
  if snrdef<4
    % Signal
    sigx1=noisex2;
    sigx2=snrwin(2);
    sigx2=seconds(sigx2-starttime);
    sigy1=-1*definfo(2);
    sigy2=definfo(2);
    %
    sigx=[sigx1 sigx1 sigx2 sigx2];
    sigy=[sigy1 sigy2 sigy2 sigy1];
    sigbox=area(sigx,sigy,sigy1);
    sigbox.ShowBaseLine='off';
    sigbox.EdgeAlpha=0;
    sigbox.FaceAlpha=0.05;
    sigbox.FaceColor=plotclr;
    sigbox.HandleVisibility='off';

    
  % SNRDEF < 4: Dotted lines for signal
  elseif snrdef==4
    % Vertical line at signal peak
    peaktim=snrwin(2);
    peaktim=seconds(peaktim-starttime);
    peaklinev=plot([peaktim peaktim],[miny maxy]);
    peaklinev.LineWidth=1;
    peaklinev.Color=plotclr;
    peaklinev.HandleVisibility='off';
    
    % Horizontal line at signal peak height across whole seismogram
    peaklineh=plot([minx maxx],[definfo(2) definfo(2)]);
    peaklineh.LineWidth=1;
    peaklineh.Color=plotclr;
    peaklineh.HandleVisibility='off';
  end
  hold on
end

% Fix Y and X axis limits
axhdl.XLim=[minx maxx];
axhdl.YLim=[miny maxy];

% List SNR in box on seismogram
snrtex=text(1,1,snrstrs);
snrtex.Units='normalized';
snrtex.FontSize=8;
snrtex.Position(1)=0.02;
snrtex.Position(2)=0.15;
hold on

% Fix Y and X axis limits
axhdl.XLim=[minx maxx];
axhdl.YLim=[miny maxy];

