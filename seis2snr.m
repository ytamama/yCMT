function [snrval,snrtims,defout]=seis2snr(sacfile,phaseortime,snrdef,winval)
% 
% Function to calculate the signal-to-noise ratio (SNR) of a particular
% arrival within a seismogram. 
% 
% INPUTS
% sacfile : SAC file, containing the arrival
% phaseortime : Two types of inputs: 
% 
%               1. Which type of arrival? Enter one of the following strings: 
%               'P', 'S', or '_CODA'
% 
%               2. Enter the onset time of arrival, as determined by 
%                  SeismoTech Ltd., in UTC.
% 
% snrdef : Which definition of SNR do we use? 
%          1 : Amplitude: max(abs(signal)) / max(abs(noise)) 
%          2 : Variance: variance(signal) / variance(noise)
%          3 : Root mean square: RMS(signal) / RMS(noise)
% 
%      [4 x] : Peak-to-noise ratio: 
%              Absolute amplitude of the first or second peak following 
%              the arrival time, that is at least x times the mean absolute
%              amplitude of the preceding noise. 
% 
% winval : For the following values of snrdef, input the following:
%          1-3: Number of seconds before and after the arrival to consider
%               as the noise and signal (e.g. 0.5)
%          4 : Number of seconds before the arrival to consider as the 
%              noise (e.g. 1s)
% 
% OUTPUT
% snrval : Value of the SNR
% snrtims : Two-element vector. The following are the elements: 
% 
%           First : Start time of the noise window, in UTC. 
%           Second : snrdef = 1,2,3: The ending time of the signal 
%                                    window, in UTC
%                    snrdef = [4 x]: The time of the chosen signal peak, 
%                                    in UTC
% 
% defout : Variable output, depending on the inputted value of snrdef: 
%          1 : [max(abs(noise)), max(abs(signal))]
%          2 : [variance(noise), variance(signal)]
%          3 : [RMS(noise), RMS(signal)]
%          4 : [mean(abs(peak height of the noise)), 
%              height of signal peak selected for SNR calculation]
%                   
% 
% Last Modified: January 5, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Retrieve snrdef values
if length(snrdef)>1
  threshval=snrdef(2);
  snrdef=snrdef(1);
end

% Initialize snrtims
snrtims=[];

% Get starting time of SAC file
[sactime,tims]=getsactime(sacfile);
% Get arrival time of phase
if ischar(phaseortime)
  arrtime=getarrtime(sacfile,phaseortime);
else
  arrtime=phaseortime;
end

% If no arrival time, return nothing
if isempty(arrtime)
  snrval=[];
  snrtims=[];
  defout=[];
  return
end
    
% If ending time of SAC file is before the arrival time, skip!
difftime=seconds(arrtime-sactime);
fintime=sactime;
fintime.Second=fintime.Second+max(tims);
if arrtime>=fintime
  snrval=[];
  snrtims=[];
  defout=[];
  return;
end

% Load SAC data
sacdata=readsac(sacfile);

% Noise window
timinds=1:length(tims);
% Indices
noiseinds=timinds(tims>=(difftime-winval) & tims<difftime);  
% Noise
noisedat=sacdata(noiseinds);  
% Time at the start of noise
noisest=arrtime;
noisest.Second=noisest.Second-winval;
snrtims=vertcat(snrtims,noisest);

if snrdef<4
  % Signal window
  sigtims=timinds(tims>=difftime & tims<=(difftime+winval));
  sigdat=sacdata(sigtims);
  % Time at the end of signal window
  sigend=arrtime;
  sigend.Second=sigend.Second+winval;
  snrtims=vertcat(snrtims,sigend);
  
  % Calculate SNR value
  % Amplitude
  if snrdef==1
    snrval=max(abs(sigdat))/max(abs(noisedat));
    defout=max(abs(noisedat));
      
  % Variance
  elseif snrdef==2
    snrval=var(sigdat)/var(noisedat);  
    defout=var(noisedat);
      
  % RMS
  elseif snrdef==3
    snrval=rms(sigdat)/rms(noisedat);  
    defout=rms(noisedat);
  end 

elseif snrdef==4
  % Determine the amplitudes correpsonding to the peak and noise, as well
  % as the time of the peak
  [pkval,pktime,noisepk]=phasepk(sacfile,phaseortime,winval,threshval);
  if isempty(pkval) || isempty(pktime) || isempty(noisepk)
    snrval=[];
    snrtims=[];
    defout=[];
    return
  end
    
  % Calculate SNR of peak : the sign of the SNR reflects whether the peak
  % is upgoing or downgoing
  snrval=pkval/noisepk;
  snrtims=vertcat(snrtims,pktime);
  defout=[noisepk pkval];
  
end


