function [solmat,inptmat,paramat]=focimt_1dvel(inputfile,velmdl,...
  jkpar,resamp)
% 
% Function to construct the output of fociMT, for the P-wave arrivals 
% recorded in the vertical component for the inputted event, 
% identified by SeismoTech Ltd. Note that all seismic data will be in 
% the displacement vs time domain.
% 
% INPUTS
% inputfile : Full path to the file used as the 'input file', in the 1D 
%             velocity format, for fociMT. We save the results of our 
%             inversion in the same path as inputfile. 
% 
% velmdl : 1-D velocity model for the region encompassing the earthquake
%          and stations, inputted as an N x 2 matrix where the first
%          column contains the depths, and the second columns contain the
%          velocity at that depth. 
%
% jkpar : Do we wish to conduct a jacknife test for this inversion, 
%         where we conduct multiple inversions in which we remove one 
%         station at a time? 
% 
%         [] : No [Default]
%         1 : Yes
% 
%         Note that if we are conducting a jacknife test, we will not 
%         conduct resampling, and vice versa. 
% 
%         Note that this flag is for conducting the jacknife test
%         programmed within fociMT. If we wish to conduct our own 
%         jacknife test where we run separate inversions, removing one 
%         station at a time, set jkpar to [] for each inversion. 
% 
% resamp : Do we wish to conduct our inversion using resampling testing, 
%          as coded in focimt.m? Input the same vector that you would 
%          use in focimt.m (see the manual). Otherwise, enter [] [default]. 
% 
%          Note that if we are conducting a jacknife test, we will not 
%          conduct resampling, and vice versa.  
% 
%    
% OUTPUTS
% solmat : Full path to the .mat file, containing the full, deviatoric, 
%          and double-couple solutions to the fociMT inversion. 
% 
% inptmat : Full path to the .mat file, containing the input parameters 
%           to the fociMT inversion. 
% 
% paramat : Full path to the .mat file, containing the parameters for the
%           fociMT inversion. 
% 
% References: 
% Kwiatek, G., Martínez‐Garzón, P. & Bohnhoff, M. (2016) HybridMT: A 
% MATLAB/Shell Environment Package for Seismic Moment Tensor Inversion and 
% Refinement. Seismological Research Letters, 87, 964–976. 
% doi:10.1785/0220150251
% 
% Last Modified : February 27, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default parameters
if ~exist('jkpar','var')
  jkpar=[];
end
if ~exist('resamp','var')
  resamp=[];
end

% Pick only one type of test to do
if ~isempty(jkpar) && ~isempty(resamp)
  jkpar=[];
end

% Directory to save our results
[focdir,~,~]=fileparts(inputfile);


% Run focimt from input file
if ~isempty(resamp) 
  [solution,input,params]=focimt(inputfile,'VelocityModel',velmdl,...
    'BeachBallFormat','PNG','Resample',resamp,'ProjectDir',focdir);

elseif ~isempty(jkpar)
  [solution,input,params]=focimt(inputfile,'VelocityModel',velmdl,...
    'BeachBallFormat','PNG','Jacknife','on','ProjectDir',focdir);

elseif isempty(jkpar) && isempty(resamp)
  [solution,input,params]=focimt(inputfile,'VelocityModel',velmdl,...
    'BeachBallFormat','PNG','ProjectDir',focdir);
end

% Retrieve and rename names of .mat files containing the solution, input, 
% and inversion parameters
solmat=fullfile(focdir,'Solution.mat');
inptmat=fullfile(focdir,'Input.mat');
paramat=fullfile(focdir,'Params.mat');



