function [pkval,pktime,noisepk]=phasepk(sacfile,phaseortime,noisewin,...
  threshval)
% 
% Function to retrieve the time and amplitude of the peak, corresponding
% to the first arrival of an inputted phase
% 
% We follow the following protocol to find this peak, given a pre-decided
% time of arrival for a given phase: 
% 
% 1. Locate the noise preceding the time of arrival. Calculate the mean
%    absolute amplitude of all the peaks in the noise. 
% 
% 2. Find the first peak following the time of arrival in the seismogram.
%    Calculate the ratio of the peak amplitude to the mean absolute 
%    amplitude of the noise, or the SNR.
% 
% 3. Decide if the SNR is at/above or below a particular threshold. 
%    If the former, then we have our peak and SNR. 
%    If the latter, then check the next peak and calculate its SNR. 
%    If Peak #2 also has a SNR below the threshold, then output the
%    SNR for Peak #1. If not, then output the SNR for Peak #2.
% 
% 
% INPUTS
% sacfile : String, containing the full path to a SAC file
% phaseortime : Two types of inputs: 
% 
%               1. Which type of arrival? Enter one of the following strings: 
%               'P', 'S', or '_CODA'
% 
%               2. Enter the time of arrival, in UTC.
% 
% noisewin : Number of seconds before the arrival of a phase to consider
%            as noise
% threshval : The threshold we use for deciding whether a peak is an
%             arrival. Input a number. 
% 
% OUTPUTS
% pkval : Amplitude of the peak, corresponding to the first arrival of 
%         the inputted phase. 
%         
% pktime : Time corresponding to the maximum amplitude of the peak, in UTC
% noisepk : Mean absolute amplitude of the noise
% 
% Last Modified : January 10, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get starting time of SAC file
[sactime,tims]=getsactime(sacfile);
% Get chosen arrival time of phase
if ischar(phaseortime)
  arrtime=getarrtime(sacfile,phaseortime);
else
  arrtime=phaseortime;
end

% If no arrival time, return nothing
if isempty(arrtime)
  pkval=[];
  pktime=[];
  noisepk=[];
  return
end
    
% If ending time of SAC file is before the arrival time, skip!
difftime=seconds(arrtime-sactime);
fintime=sactime;
fintime.Second=fintime.Second+max(tims);
if arrtime>=fintime
  pkval=[];
  pktime=[];
  noisepk=[];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify peaks in noise data

% Noise
sacdata=readsac(sacfile);
timinds=1:length(tims);
noiseinds=timinds(tims>=(difftime-noisewin) & tims<difftime);  
noisedat=sacdata(noiseinds);

% Peaks
try
  [npks1,nlocs1]=findpeaks(noisedat);
catch
  keyboard
end
%
try
  [npks2,nlocs2]=findpeaks(-1*noisedat);
catch
  keyboard 
end
noisemat=horzcat(vertcat(abs(npks1),-1*npks2),vertcat(nlocs1,nlocs2));
% Exit if no peaks
if isempty(noisemat)
  pkval=[];
  pktime=[];
  noisepk=[];
  return
end

% Calculate the mean absolute peak amplitude of noise
noisepks=noisemat(:,1);
noisepk=mean(abs(noisepks));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arrival

% Data following arrival time
arrinds=timinds(tims>=difftime);  
arrtims=tims(arrinds);
arrdat=sacdata(arrinds);

% Peaks
try
  [apks1,alocs1]=findpeaks(arrdat);
catch
  keyboard
end
%
try
  [apks2,alocs2]=findpeaks(-1*arrdat);
catch
  keyboard 
end
arrmat=horzcat(vertcat(abs(apks1),-1*apks2),vertcat(alocs1,alocs2));
% Exit if no peaks
if isempty(arrmat)
  pkval=[];
  pktime=[];
  noisepk=[];
  return
end
arrmat=sortrows(arrmat,2);

% Find the first peak following the picked arrival time
pk1=arrmat(1,1);
loc1=arrmat(1,2);
% Check if the signal-to-noise ratio exceeds the threshold
if abs(pk1)/noisepk>=threshval
  pkval=pk1;
  pktime=sactime+seconds(arrtims(loc1));
  return
end

% Otherwise, check the second peak
pk2=arrmat(2,1);
loc2=arrmat(2,2);
% Check if the signal-to-noise ratio exceeds the threshold
if abs(pk2)/noisepk>=threshval
  pkval=pk2;
  pktime=sactime+seconds(arrtims(loc2));
  return
  
% Otherwise, return the larger peak
else
  if abs(pk2)>abs(pk1)
    pkval=pk2;
    pktime=sactime+seconds(arrtims(loc2));
    return
  else
    pkval=pk1;
    pktime=sactime+seconds(arrtims(loc1));
    return
  end
end


