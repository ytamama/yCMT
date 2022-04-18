function sacfile=ms2sac(msfile)
% 
% Function to convert a miniseed file to a SAC file, located in the 
% current working directory. 
% 
% INPUT
% msfile : Miniseed file and its full path to it
% 
% OUTPUT
% sacfile : SAC file, located in the current working directory
% 
% Last Modified : December 22, 2021 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert miniseed to SAC
[status,cmdout]=system(sprintf('mseed2sac %s',msfile));
outstrs=strsplit(cmdout);
% Retrieve name of SAC file
for s=1:length(outstrs)
  nowstr=outstrs{s};
  if contains(nowstr,'.SAC')
    sacfile=nowstr;
  end
end

% If SAC file name not retrieved
if ~exist('sacfile','var')
  sacfile='';
end



