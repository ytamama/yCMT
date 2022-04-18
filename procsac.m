function newsac=procsac(oldsac,tfreqs,pzfile,ffreqs,remflag,valtype)
% 
% Function to process a SAC file, including removing the mean and trend,
% filtering, and deconvolving
%
% INPUT
% oldsac : Unfiltered, unprocessed SAC file to filter, deconvolve, 
%          detrend, etc. 
% tfreqs : Corner frequencies in Hz, through which we deconvolve our 
%          data. Define a four-element vector with frequencies f1, f2, f3,
%          f4 in Hz, where f2 is at least twice f1, and f4 at least twice 
%          f3. Deconvolve the seismograms if we enter this vector. If we 
%          don't wish to do so, enter an empty vector. 
% pzfile : Full path to the SACPZ file, with which we deconvolve our 
%          SAC file if necessary. Enter an empty string if not
%          deconvolving
% ffreqs : Frequencies through which we filter our data. Define a 2 
%          element vector, where the first one is our lower bound and the
%          second one is our upper bound. Both in Hz. 
% remflag : Remove mean and linear trend? 
%           0 : No
%           1 : Yes - must remove for deconvolution!
% valtype : Do we deconvolve our SAC file into displacement, velocity, or
%           acceleration?
%           0 : Displacement 
%           1 : Velocity 
%           2 : Acceleration 
%           [] : Not deconvolving
% 
% OUTPUT
% newsac : Filtered and processed SAC file, in the same directory as 
%          oldsac. 
% 
% Note: If both tfreqs and ffreqs are empty, return the original SAC file. 
%       If only tfreqs has values, return a deconvolved SAC file
%       If only ffreqs has values, return a filtered SAC file
%       If both tfreqs and ffreqs have values, return a deconvolved and
%       filtered SAC file
%       
% 
% Last Modified: October 20, 2021 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify new file name
part1=oldsac(1:end-4);
newsac=sprintf('%s_PROC.SAC',part1);

% Nothing to do... 
if isempty(tfreqs) && isempty(ffreqs)
  fprintf('Returning original SAC file \n')
  newsac=oldsac;
  return
  
% Only deconvolve
elseif ~isempty(tfreqs)
  % Transfer function command
  transfer=sprintf(...
    'transfer from polezero subtype %s to %s freqlimits %g %g %g %g prewhitening on',...
    pzfile,valtype,tfreqs(1),tfreqs(2),tfreqs(3),tfreqs(4));

  % Deconvolution
  dccmd=sprintf('echo "r %s ; read ; rtr ; rmean ; taper type ; %s ; w %s ; q" | /usr/local/sac/bin/sac',...
    oldsac,transfer,newsac);
  [status,cmdout]=system(dccmd);
  if status>0
    keyboard
  end
    
% Only filter
elseif ~isempty(ffreqs)
  % Filter 
  if remflag==0
    fcmd=sprintf('echo "r %s ; bp c %g %g ; w %s ; q" | /usr/local/sac/bin/sac',...
      oldsac,ffreqs(1),ffreqs(2),newsac);
  else
    fcmd=sprintf('echo "r %s ; rtr ; rmean ; bp c %g %g ; w %s ; q" | /usr/local/sac/bin/sac',...
      oldsac,ffreqs(1),ffreqs(2),newsac);
  end
  [status,cmdout]=system(fcmd);
  if status>0
    keyboard
  end
    
% Deconvolve and filter
else
  % Transfer function command
  transfer=sprintf(...
    'transfer from polezero subtype %s to %s freqlimits %g %g %g %g prewhitening on',...
    pzfile,valtype,tfreqs(1),tfreqs(2),tfreqs(3),tfreqs(4));

  % Deconvolution
  dccmd=sprintf('echo "r %s ; read ; rtr ; rmean ; taper type ; %s ; w %s ; q" | /usr/local/sac/bin/sac',...
    oldsac,transfer,newsac);
  [status,cmdout]=system(dccmd);
  if status>0
    keyboard
  end

  % Filter 
  fcmd=sprintf('echo "r %s ; bp c %g %g ; w %s ; q" | /usr/local/sac/bin/sac',...
    newsac,ffreqs(1),ffreqs(2),newsac);
  [status,cmdout]=system(fcmd);
  if status>0
    keyboard
  end
end


