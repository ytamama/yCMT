function inputfile=bldinpt(evthdr,stanums,omegas,staes,stans,stazs,...
  issynth,savedir,minmaxfreq)
% 
% Function to write the 1D velocity input file to use for FociMT to 
% derive the focal mechanism of an event
% 
% INPUTS
% evthdr : Information to add to the header of the file, entered as a 
%          cell with the following entries : 
% 
%          evtnum : ID number corresponding to the event
%          evtn : Northing of event location, in m
%          evte : Easting of event location, in m
%          evtz : z-coordinate (altitude) of event location, in m. This 
%                 value should increase with increasing altitude. 
%          rho : Density of area, in kg/m3
% 
% stanums : Cell array, containing ID numbers of stations
% omegas : Vector containing areas under first arriving P-wave in the 
%          Z-component, for each station in stanums, in m*s
% staes : Easting of each station in stanums, in m
% stans : Northing of each station in stanums, in m
% stazs : z-coordinate (altitude) of each station in stanums, in m. This 
%         value should increase with increasing altitude.
% 
% issynth : Are we making an input file for a synthetic inversion? If so, 
%           then for what type of inversion? 
% 
%           0 : Not synthetic, real data
%           1 : Synthetic full
%           2 : Synthetic deviatoric
%           3 : Synthetic double-couple
% 
% savedir : Directory to save our resulting input file, inputted as a 
%           string.
% 
% minmaxfreq : Minimum and maximum frequencies at which we filter our 
%              seismic data. Enter [] if that information is unknown. '
%              Otherwise, enter a 2-element vector, with the first 
%              element being the minimum frequency and the second element
%              being the maximum frequency, both in Hz. 
% 
% OUTPUT
% inputfile : Full path to the resulting input file
% 
% Last Modified : March 16, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack event header information
evtnum=evthdr{1};
evtn=evthdr{2};
evte=evthdr{3};
evtz=evthdr{4};
rho=evthdr{5};

% Frequencies?
if ~isempty(minmaxfreq)
  freqstr=sprintf('_%.2gTO%.2gHz',minmaxfreq(1),minmaxfreq(2));
else
  freqstr='';
end

% Open resulting input file
% Real data
if issynth==0
  inputfile=sprintf('EVT%d.P.Z.vel1dinput%s.txt',evtnum,freqstr);

% Synthetic full
elseif issynth==1
  inputfile=sprintf('EVT%d.P.Z.vel1dinput%s_FULLSYNTH%s.txt',evtnum,...
    freqstr);
  
% Synthetic deviatoric
elseif issynth==2
  inputfile=sprintf('EVT%d.P.Z.vel1dinput%s_DEVSYNTH.txt',evtnum,freqstr);
  
% Synthetic double-Couple
elseif issynth==3
  inputfile=sprintf('EVT%d.P.Z.vel1dinput%s_DCSYNTH.txt',evtnum,freqstr);
end
inputfile=fullfile(savedir,inputfile);
fid=fopen(inputfile,'w');

% Write event header
numsta=length(stanums);
fprintf(fid,'%d %d %g %g %g %g\n',evtnum,numsta,evtn,evte,evtz,rho);


% Iterate through each station and write its information to input file
for s=1:numsta
  % Retrieve info
  nowsta=stanums(s);
  nowomega=omegas(s);
  nowomega=round(nowomega,14,'significant');
  nowe=staes(s);
  nown=stans(s);
  nowz=stazs(s);
    
  % Station number
  if nowsta<10
    stastr=sprintf('S00%d',nowsta);
  elseif nowsta<100
    stastr=sprintf('S0%d',nowsta); 
  else
    stastr=sprintf('S%d',nowsta);
  end
  
  % Construct line and add to input file
  fprintf(fid,'%s Z P %0.13e %g %g %g\n',stastr,nowomega,nown,nowe,nowz);
end

% Close file
fclose(fid);


