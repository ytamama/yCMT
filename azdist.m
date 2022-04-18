function [az,dist]=azdist(evcoord,stacoord,calcmode)
% 
% Function to compute the azimuth and distance between an event and 
% station.
% 
% INPUTS
% evcoord : Event coordinates
%           [Longitude, Latitude] in degrees, or 
%           [X, Y] in km, with respect to a local origin
%           [Easting, Northing] in m
% 
% stacoord : Station coordinates
%            [Longitude, Latitude] in degrees, or 
%            [X, Y] in km, with respect to a local origin
%            [Easting, Northing] in m
% 
% calcmode : Do we calculate the azimuth and distance in Cartesian 
%            coordinates, or do we use spherical coordinates with a 
%            Great Circle approximation? 
% 
%            1 : Cartesian, with an X-Y coordinate system developed for 
%                events and stations recorded by SeismoTech Ltd. 
% 
%                OR 
% 
%                Cartesian, with Northing and Easting coordinates in m
% 
%            2 : Spherical
% 
% OUTPUTS
% az : Azimuth from event to station, in degrees and clockwise from North
% dist : Distance between event and station, in km or m for calcmode=1 and 
%        in degrees for calcmode=2
% 
% Last Modified : February 2, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian
if calcmode==1 
  % "Unpack" coordinates
  evx=evcoord(1);
  evy=evcoord(2);
  stx=stacoord(1);
  sty=stacoord(2);
    
  % Distance
  dist=sqrt((evx-stx)^2 + (evy-sty)^2);
    
  % Azimuth
  az=atand((stx-evx)/(sty-evy));
  % Quadrant II
  if (az<0) && ((stx-evx)<0)
    az=360+az;
  
  % Quadrant IV
  elseif (az<0) && ((sty-evy)<0)
    az=az+180;
    
  % Quadrant III
  elseif (az>0) && ((sty-evy)<0)
    az=az+180;
    
  % Directly south of the event 
  elseif (az==0) && ((sty-evy)<0)
    az=-180;
  end
    
    
% Spherical
elseif calcmode==2
  % "Unpack" coordinates
  evlo=evcoord(1);
  evla=evcoord(2);
  stlo=stacoord(1);
  stla=stacoord(2);
    
  % Calculate Great Circle distance and azimuth
  [dist,az]=distance(evla,evlo,stla,stlo);

end

