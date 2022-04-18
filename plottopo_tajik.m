function [mapfig,mapax,figname]=plottopo_tajik(maprange,addfaults,...
    addbrdr,makefig,savedir)
% 
% Function to plot the topography across Tajikistan or the Tajik Basin, 
% using data provided by GEBCO 2021
% 
% INPUT    
% maprange : Do we plot across Tajikistan, or just the Tajik Basin? 
%            0 : Across Tajikistan + parts of neighboring countries
%            1 : Across Tajikistan + parts of neighboring countries, but
%                with a box around our area of study
%            2 : Just our area of study
% 
%            [minlon maxlon minlat maxlat] : A custom input of, in this 
%                                            order, minimum longitude, 
%                                            maximum longitude, minimum 
%                                            latitude, maximum latitude. 
% 
% addfaults : Do we wish to add fault lines? 
%             0 : No
%             1 : Yes, but don't color-code
%             2 : Yes, and color-code
%             3 : Yes, but just the major faults
% addbrdr : Add national borders? 
%           0 : No
%           1 : Yes
% makefig : Make new figure? Or plot on an existing figure? 
%           0 : Plot on existing figure
%           1 : Make new figure
% savedir : Full path to the directory where we wish to save our plots. 
%           Enter an empty string if we don't wish to save our plot. 
% 
% OUTPUTS
% mapfig : Figure handle to our topographic map
% mapax : Axis handle to our topographic map
% figname : Full path to our saved figure. Empty string we didn't save
%           our plot. 
% 
% Topography data from GEBCO Compilation Group. (2021) GEBCO 2021 Grid. 
% doi:10.5285/c6612cbe-50b3-0cff-e053-6c86abc09f8f
% 
% Code uses: 
% 1. Code from bathymetry.m to read the elevation CDF file and generate
%    a colormap, in MERMAID_buffer/bathymetry.m by Sirawich (Pete) 
%    Pipatrathanporn
% 2. cax2dem.m in csdms-contrib/slepian_alpha by Frederik Simons
% 3. sergeicol.m in csdms-contrib/slepian_alpha, using a colormap by 
%    Sergei Lebedev
% 
% NetCDF Reference:
% Unidata. (2012) Network Common Data Form (NetCDF) version 4.3.3.1 
% [software], Boulder, CO: UCAR/Unidata. DOI: 10.5065/D6H70CW6
% 
% Last Modified: February 27, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load topography data
topofile=gettopofile();
lonvals=ncread(topofile,'lon');
latvals=ncread(topofile,'lat');
elevals=ncread(topofile,'elevation');
% Minimum and maximum elevation
minele=min(min(elevals));
maxele=max(max(elevals));

% Initialize figure
if makefig==1
  mapfig=figure();
  mapfig.Units='normalized';
  mapfig.Position(3)=0.5;
  mapfig.Position(4)=0.75;
else
  mapfig=gcf;
end

% Plot elevation!
eleplot=imagesc(lonvals,latvals,elevals.');
axis xy;
[cb,~]=cax2dem([minele maxele],'hor');
cb.Location='southoutside';
cb.Label.String='Elevation Above Sea Level (m)';
cb.Label.FontSize=14;
grid on
hold on

% Plot national borders
if addbrdr==1
  brdrlines=plotbrdr_tajik([0 0 0]);
end

% Load earthquake data to impose lon-lat limits
seistbl=getevtinfo();
lonpts=seistbl.LON_dec;
latpts=seistbl.LAT_dec;

% Plot across Tajikistan + parts of neighboring countries
if maprange<2
  minlat=36;
  maxlat=40;
  minlon=66;
  maxlon=72;
  if maprange==1
    minlat_box=min(latpts)-0.1;
    maxlat_box=max(latpts)+0.1;
    minlon_box=min(lonpts)-0.1;
    maxlon_box=max(lonpts)+0.1;
  end
% Plot across Tajik Basin
else
  minlon=min(lonpts)-0.1;
  maxlon=max(lonpts)+0.1;
  minlat=min(latpts)-0.1;
  maxlat=max(latpts)+0.1;
end

% Plot faults
if addfaults>0
  if addfaults==1
    faulttypes={'normal';'thrust';'strike';'subsurface';'uncertain';...
      'uncertain subs';'major'};
    clrcode=0;
  elseif addfaults==2
    faulttypes={'normal';'thrust';'strike';'subsurface';'uncertain';...
      'uncertain subs';'major'};
    clrcode=1;
  else
    faulttypes={'major'};
    clrcode=1;
  end
  plotfaults_tajik(faulttypes,clrcode)
end

% Axes limits and labels
mapax=gca;
mapax.XLim=[minlon maxlon];
mapax.YLim=[minlat maxlat];  
mapax.FontSize=14;
mapax.XLabel.String=sprintf('Longitude (%s)',char(176));
mapax.XLabel.FontSize=15;
mapax.YLabel.String=sprintf('Latitude (%s)',char(176));
mapax.YLabel.FontSize=15;
if addfaults>0
  mapax.Title.String='Faults and Topography of the Tajik Basin';
else
  mapax.Title.String='Topography of the Tajik Basin';  
end
hold on

% Plot box enclosing our area of study, if needed
if maprange==1
  seiswid=maxlon_box-minlon_box;
  seisht=maxlat_box-minlat_box;
  seisboxplt=rectangle('Position',[minlon_box minlat_box seiswid seisht]);
  seisboxplt.FaceColor='none';
  seisboxplt.EdgeColor=[0 0 0];
  seisboxplt.LineWidth=3;
end

% Save figure
if maprange==0
  mapstr='tajikistan.all';
elseif maprange==1
  mapstr='tajikistan.all.box';
else
  mapstr='tajikbasin';
end
%
if addfaults==0
  faultstr='nofault';
elseif addfaults==1
  faultstr='faults';
else
  faultstr='clrcodefaults';
end
%
if ~isempty(savedir)
  figname=sprintf('tajikbasintopo_%s_%s.png',mapstr,faultstr);
  figname=fullfile(savedir,figname);
  saveas(mapfig,figname);
else
  figname='';
end

