function [inputfile,nostas]=focimt_1dvel_inpt(evtnum,rho,procinfo,...
  threshval,noisewin,savedir,jkpar)
% 
% Function to construct the 1D velocity model input file to use in 
% FociMT to reconstruct the focal mechanism for an event recorded by 
% SeismoTech Ltd.
% 
% INPUTS
% evtnum : ID number corresponding to the event
% rho : Approx. density of the region, in kg/m3. In a thermodynamic model
%       of the Tajik Basin by Chapman et al. 2017, doi: 10.1111/bre.12381
%       the density is set to be 2500 kg/m3.
% 
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
% savedir : Directory to save our resulting input file, inputted as a 
%           string. 
% 
% jkpar : Do we wish to run our inversion with a jacknife test, where we
%         run multiple inversions where we remove just one station at a 
%         time? 
% 
%         [] : No [default]
%         1 : Yes
% 
% OUTPUT
% inputfile : if jkpar = [] or nonexistent, inputfile is a string 
%             containing the full path to the input file used in fociMT. 
% 
%             if jkpar = 1, inputfile is a cell array of strings, 
%             containing the full path to an input file using all but one 
%             station. We have n files total, where n is the number of 
%             stations. 
% 
% nostas : [] if jkpar = [] or nonexistent. 
% 
%          if jkpar = 1, then nostas is a vector the same dimensions as 
%          inputfile, where each value is the station that was excluded in 
%          composing each entry of inputfile. 
%          
% References: 
% Kwiatek, G., Martínez‐Garzón, P. & Bohnhoff, M. (2016) HybridMT: A 
% MATLAB/Shell Environment Package for Seismic Moment Tensor Inversion and 
% Refinement. Seismological Research Letters, 87, 964–976. 
% doi:10.1785/0220150251
% 
% Last Modified : March 6, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default for jkpar
if ~exist('jkpar','var')
  jkpar=[];
end

% Retrieve event information
evtrow=getevtinfo(evtnum);
% Easting and Northing coordinates, in m
evte_val=evtrow.EASTING;
evtn_val=evtrow.NORTHING; 
evte=0;
evtn=0;
% Z should increase in value with increasing altitude
evtz=-1000*evtrow.z_km;
% Cell array, containing event header information
evthdr={evtnum; evtn; evte; evtz; rho};

% Retrieve P-wave picks in the Z-component!
sacvars={evtnum;[];'Z'};
arrtypesin={'P'};
picktbl=getpickinfo(sacvars,arrtypesin);
[numpicks,~]=size(picktbl);

% Iterate through P-wave picks and assemble input file information
stanums=[];
omegas=[];
staes=[];
stans=[];
stazs=[];
snrvals=[];
% Calculate areas under the P-waves
for p=1:numpicks
  % Arrival information
  pickrow=picktbl(p,:);
  nowsta=pickrow.Station;
  nowcomp=pickrow.Component{1};
  if ~strcmp(nowcomp,'Z')
    continue
  end
  nowarr=pickrow.ArrivalType{1};
  if ~strcmp(nowarr,'P')
    continue
  end
  % Time of onset of P-wave, as chosen by SeismoTech Ltd. 
  starrtime=pickrow.ArrivalTime;
  
  % Station coordinates, in m
  starow=getstainfo(nowsta);
  stae=round(starow.Easting-evte_val,4);
  stan=round(starow.Northing-evtn_val,4);
  % Z-coordinate should increase with higher altitude
  staz=-1000*starow.Altitude;
  
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
  
  % Calculate 'omega' (Area under first arriving P-wave, with sign)
  [omega,~,~]=phasearea(sacfile,starrtime,noisewin,threshval);
  % Convert omega from nm * s to m * s
  omega=omega*(1e-9);
  
  % Populate arrays! 
  stanums=vertcat(stanums,nowsta);
  staes=vertcat(staes,stae);
  stans=vertcat(stans,stan);
  stazs=vertcat(stazs,staz);
  omegas=vertcat(omegas,omega);
  snrvals=vertcat(snrvals,snrval);
end


% IF we have at least 132 stations as input, fociMT somehow refuses 
% to work... The upper threshold seems to be 128 stations. 
% SO, remove the stations with the lowest P-wave SNR until we have 128
% stations: 
if length(stanums)>=132
  snrvals=abs(snrvals);
  inptmat=horzcat(snrvals,stanums,staes,stans,stazs,omegas);
  inptmat=sortrows(inptmat,1,'descend');
  inptmat2=inptmat(1:128,:);
  % 
  stanums=inptmat2(:,2);
  staes=inptmat2(:,3);
  stans=inptmat2(:,4);
  stazs=inptmat2(:,5);
  omegas=inptmat2(:,6);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct 1D velocity model input file

% Frequencies through which we filtered the seismograms
ffreqs=procinfo{4};

% No jacknife test
try
  inputfile=bldinpt(evthdr,stanums,omegas,staes,stans,stazs,...
    0,savedir,ffreqs);
catch
  keyboard
end

% No jacknife test
if isempty(jkpar)
  nostas=[];
  return;

% Jacknife test
elseif jkpar==1

  inputfile={inputfile};
  nostas=[];
  
  % Iterate through all stations
  for s=1:length(stanums)
    % Station number we "kick out"
    missingsta=stanums(s);
    nostas=vertcat(nostas,missingsta);
    
    % Inputs to input file, minus that corresponding to the station we 
    % "kick out"
    stanums_new=stanums(stanums~=missingsta);
    omegas_new=omegas(stanums~=missingsta);
    staes_new=staes(stanums~=missingsta);
    stans_new=stans(stanums~=missingsta);
    stazs_new=stazs(stanums~=missingsta);
    
    % Put each test in its own folder 
    stafldr=sprintf('NOSTA_%d',missingsta);
    savedir_new=fullfile(savedir,stafldr);
    if exist(savedir_new)==0
      mkdir(savedir_new);
    end
    
    % Build input file
    try
      inptfile_new=bldinpt(evthdr,stanums_new,omegas_new,staes_new,...
        stans_new,stazs_new,0,savedir_new,ffreqs);     
    catch
      keyboard 
    end
    
    % Add to output array
    inputfile=vertcat(inputfile,inptfile_new);
  end  
end



