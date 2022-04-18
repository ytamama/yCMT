function [sactime,tims]=getsactime(sacfile)
% 
% Function to retrieve the starting time of a SAC file
% 
% INPUT
% sacfile : A SAC file
% 
% OUTPUT
% sactime : Datetime object, representing the starting time of the 
%           inputted SAC file
% tims : Values on the time-axis, should the seismogram be plotted. 
%        Output from csdms-contrib/slepian_oscar/readsac.m
% 
% Last Modified: October 23, 2021 by Yuri Tamama
% 
% Uses readsac.m in csdms-contrib/slepian_oscar
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read SAC file
[~,hdrinfo,~,~,tims]=readsac(sacfile);

% Starting time of SAC file 
startyr=hdrinfo.NZYEAR;
startjday=hdrinfo.NZJDAY;
starthr=hdrinfo.NZHOUR;
startmin=hdrinfo.NZMIN;
startsec=hdrinfo.NZSEC+(hdrinfo.NZMSEC/1000);
startymd=jul2dat(startyr,startjday);
startmon=startymd(1);
startday=startymd(2);
sactime=datetime(startyr,startmon,startday,starthr,startmin,startsec);

