function [M11,M12,M13,M22,M23,M33,mtensor]=ftang2mt(str,dip,rake,M0,...
  convention)
% 
% Function to derive the six independent components, in Cartesian 
% coordinates, of the moment tensor, based on the scalar seismic moment
% and strike, dip, and rake angles. 
% 
% INPUTS
% str : Strike angle, in degrees
% dip : Dip angle, in degrees
% rake : Rake angle, in degrees
% M0 : Scalar seismic moment, in N*m
% convention : 1 for XYZ, 2 for RTF
% 
% In either convention, the moment tensor would be of this form:
% 
% M11 M12 M13
% M12 M22 M23
% M13 M23 M33
% 
% References: 
% Chapter 4 of Aki, K. & Richards, P.G. (2002) Quantitative Seismology, 
% 2nd ed., Mill Valley, CA: University Science Books.
% 
% Note that this formula is only applicable for double-couple sources. 
% 
% Last Modified: April 5, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% XYZ
if convention==1
  M11 = -1*M0*(sind(dip)*cosd(rake)*sind(2*str) +...
    sind(2*dip)*sind(rake)*sind(str)*sind(str));
  M12 = M0*(sind(dip)*cosd(rake)*cosd(2*str) +...
    0.5*sind(2*dip)*sind(rake)*sind(2*str));
  M13 = -1*M0*(cosd(dip)*cosd(rake)*cosd(str) +...
    cosd(2*dip)*sind(rake)*sind(str));
  M22 = M0*(sind(dip)*cosd(rake)*sind(2*str) -...
    sind(2*dip)*sind(rake)*cosd(str)*cosd(str));
  M23 = -1*M0*(cosd(dip)*cosd(rake)*sind(str) -...
    cosd(2*dip)*sind(rake)*cosd(str));
  M33 = M0*sind(2*dip)*sind(rake);

  
% RTF
elseif convention==2
  M11 =  M0*sind(2*dip)*sind(rake);
  M22 =  -M0*(sind(dip)*cosd(rake)*sind(2*str) +...
    sind(2*dip)*sind(rake)*(sind(str))^2 );
  M33 =  M0*(sind(dip)*cosd(rake)*sind(2*str) -...
    sind(2*dip)*sind(rake)*(cosd(str))^2 );
  M12 =  -M0*(cosd(dip)*cosd(rake)*cosd(str)  +...
    cosd(2*dip)*sind(rake)*sind(str) );
  M13 =  M0*(cosd(dip)*cosd(rake)*sind(str)  -...
    cosd(2*dip)*sind(rake)*cosd(str) );
  M23 =  -M0*(sind(dip)*cosd(rake)*cos(2*str) +...
    0.5*sind(2*dip)*sind(rake)*sind(2*str) );
end
    
% Moment Tensor
mtensor=[M11 M12 M13; M12 M22 M23; M13 M23 M33];


