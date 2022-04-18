function [rsac,tsac]=rotatesac(esac,nsac)
% 
% Function to rotate SAC files in the North and East directions into
% those of the radial and transverse directions. These SAC files are of
% events recorded across Tajikistan by the TOTAL array, established by 
% SeismoTech Ltd. 
% 
% INPUTS
% esac : SAC file in the East component
% nsac : SAC file in the North component
% 
% OUTPUTS 
% rsac : SAC file in the radial component
% tsac : SAC file in the transverse component
% 
% Note: assume that the East and North component SAC files already
%       have the coordinates of station and event in their headers. 
% 
% Based on rotsac.m in https://github.com/ytamama/GuyotSeismology
% 
% Last Modified: September 13, 2021 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set names for new SAC files
rsac=replace(esac,'.E.','.R.');
tsac=replace(esac,'.E.','.T.');

% Set component inclination to 90 degrees for both files
inccmd_e=sprintf(...
  'echo "r %s ; chnhdr CMPINC 90 ; w %s ; q" | /usr/local/sac/bin/sac',...
  esac,esac);
[status,cmdout]=system(inccmd_e);
%
inccmd_n=sprintf(...
  'echo "r %s ; chnhdr CMPINC 90 ; w %s ; q" | /usr/local/sac/bin/sac',...
  nsac,nsac);
[status,cmdout]=system(inccmd_n);
% Define CMPAZ in both components
azcmdn=sprintf(...
  'echo "r %s ; chnhdr CMPAZ 0 ; w %s ; q" | /usr/local/sac/bin/sac',...
  nsac,nsac);
[status,cmdout]=system(azcmdn);
azcmde=sprintf(...
  'echo "r %s ; chnhdr CMPAZ 90 ; w %s ; q" | /usr/local/sac/bin/sac',...
  esac,esac);
[status,cmdout]=system(azcmde);
% Rotate from XY to RT
rotcmd=sprintf(...
  'echo "r %s %s ; rotate to GCP ; w %s %s ; q" | /usr/local/sac/bin/sac',...
  nsac,esac,rsac,tsac);
[status,cmdout]=system(rotcmd);


% Assign component name
compcmdr=sprintf(...
  'echo "r %s ; chnhdr KCMPNM R ; w %s ; q" | /usr/local/sac/bin/sac',...
  rsac,rsac);
[status,cmdout]=system(compcmdr);
compcmdt=sprintf(...
  'echo "r %s ; chnhdr KCMPNM T ; w %s ; q" | /usr/local/sac/bin/sac',...
  tsac,tsac);
[status,cmdout]=system(compcmdt);


