function [mapfig,figname,allax]=ploteqs_tajik(evtnums,maprange,addfaults,...
  eqcode,addbrdr,makefig,savedir)
% 
% Function to plot the locations of earthquakes recorded by the TOTAL SA 
% array of SeismoTech Ltd. in 2015, across the Tajik Basin. 
% 
% INPUTS
% evtnums : Do we plot a specific subset of seismic events? If so, enter 
%           a vector containing their ID numbers! If we wish to plot all 
%           of them instead, input nothing (Default)
% maprange : Do we plot across Tajikistan, or just the Tajik Basin? 
%            0 : Across Tajikistan + parts of neighboring countries
%            2 : Just our area of study
% addfaults : Do we wish to add fault lines?
%             0 : No
%             1 : Yes, but don't color-code
%             2 : Yes, and color-code
% eqcode : How do we color-code the earthquakes? Enter a two-element
%          vector, where the first element is how we color-code: 
% 
%          1 : Change in longitude relative to initial estimate
%          2 : Change in latitude relative to initial estimate
%          3 : Change in depth relative to initial estimate
%          4 : Lateral (lat-lon) distance between final and initial
%              estimate
%          5 : 3D Distance between final and initial estimate
%          6 : Depth
%          7 : Magnitude
% 
%          and the second element tells us whether or not we plot all 
%          earthquakes: 
%          0 : Plot all earthquakes
%          1 : Plot only the earthquakes within the top 5% of the value
%              specified in the first element. 
%         -1 : Plot only the earthquakes within the bottom 5% of the value
%              specified in the first element. 
% 
%          Enter [] (empty vector) if we're not color-coding
% 
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
% figname : Full path to our saved figure. Empty string we didn't save
%           our plot. 
% allax : Arrat of all plot axes
% 
% Topography ata from GEBCO Compilation Group. (2021) GEBCO 2021 Grid. 
% doi:10.5285/c6612cbe-50b3-0cff-e053-6c86abc09f8f
% 
% Code uses: 
% 1. Code from bathymetry.m to read the elevation CDF file and generate
%    a colormap, in MERMAID_buffer/bathymetry.m by Sirawich (Pete) 
%    Pipatrathanporn
% 2. cax2dem.m in csdms-contrib/slepian_alpha by Frederik Simons
% 3. sergeicol.m in csdms-contrib/slepian_alpha, using a colormap by 
%    Sergei Lebedev
% 4. figdisp.m in csdms-contrib/slepian_alpha by Frederik Simons
% 5. A modified version of redblueu.m, by Mirko Hrovat (2021)
% 
% References: 
% Mirko Hrovat (2021). 
% RedBlue Colormap Generator with Zero as White or Black 
% (https://www.mathworks.com/matlabcentral/fileexchange/74791-redblue-colormap-generator-with-zero-as-white-or-black), 
% MATLAB Central File Exchange. Retrieved November 18, 2021.
% 
% Unidata. (2012) Network Common Data Form (NetCDF) version 4.3.3.1 
% [software], Boulder, CO: UCAR/Unidata. DOI: 10.5065/D6H70CW6
% 
% Uses defval.m and figdisp.m in csdms-contrib/slepian_alpha
% 
% Last Modified: February 27, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default value
defval('evtnums',[])

% Retrieve values from eqcode
if ~isempty(eqcode)
  codeval=eqcode(1);
  plotall=eqcode(2);
else
  codeval=0;
end

% Plot map of topography and faults
[mapfig,]=plottopo_tajik(maprange,addfaults,addbrdr,makefig,'');
hold on 

% Obtain oordinates of earthquakes to plot
seistbl=getevtinfo();
% Filter by event ID numbers if applicable
if ~isempty(evtnums)
  seistbl=seistbl(ismember(seistbl.EVNT,evtnums),:);
end
% Longitude and latitude coordinates
lonpts=seistbl.LON_dec;
latpts=seistbl.LAT_dec;
numevents=length(lonpts);

% Sort and retrieve variable to color code
if codeval>0
  % Depth
  if codeval==6
    % Sort table by absolute depth
    absdep=table(abs(seistbl.DEP_KM));
    seistbl=[seistbl, absdep];
    seistbl=sortrows(seistbl,22,'ascend');
    
  % Magnitude
  elseif codeval==7
    % Sort table by absolute magnitude
    absmag=table(abs(seistbl.MAG));
    seistbl=[seistbl, absmag];
    seistbl=sortrows(seistbl,22,'ascend');
    
  % Final - Initial Coordinates
  elseif codeval<6
    [diffvals,diffids]=getdiffloc(codeval);
    % Sort differences by ascending event ID
    difftbl=table(diffids,diffvals);
    difftbl=sortrows(difftbl,1,'ascend');
    diffids=difftbl.diffids;
    diffvals=difftbl.diffvals;
    % Sort event table by coordinate or distance differences
    difftbl2=table(diffvals,abs(diffvals));
    seistbl=[seistbl, difftbl2];
    seistbl=sortrows(seistbl,23,'ascend');
  end
  minval=min(seistbl{:,22});
  maxval=max(seistbl{:,22});
  clim=[minval maxval];
  valrange=maxval-minval;
  
  % All values
  allvals=seistbl{:,22};
end


% Lon-lat limits 
if maprange<2
  minlat=36;
  maxlat=40;
  minlon=66;
  maxlon=72;
else
  minlon=min(lonpts)-0.1;
  maxlon=max(lonpts)+0.1;
  minlat=min(latpts)-0.1;
  maxlat=max(latpts)+0.1;
end

% Axes limits
nowax=gca;
allax=[];
allax=horzcat(allax,nowax);
nowax.XLim=[minlon maxlon];
nowax.YLim=[minlat maxlat]; 

% Plot Title
titlestrs={'Change in Longitude';'Change in Latitude';'Change in Depth';...
  'Change in 2D Distance';'Change in 3D Distance';'Depth';'Magnitude'};
if codeval==0
  titlestr1='Earthquakes Recorded by TOTAL Seismic Array';
  titlestr2=sprintf('%d Events',numevents);
  nowax.Title.String={titlestr1;titlestr2};
else
  nowax.Title.String=sprintf('%d Events: Color-Coded by %s',...
    numevents,titlestrs{codeval});
end


% Colors for color-coding earthquake
if codeval>0
  % Use modified redblueu.m
  eqclrs=redblueu_MOD(250,clim);
  [numclr,~]=size(eqclrs);
end


% Plot earthquakes
for i=1:numevents
  nowrow=seistbl(i,:); 
  nowlon=nowrow.LON_dec;
  nowlat=nowrow.LAT_dec;
  if codeval>0 && codeval<6
    nowval=nowrow{:,22};
  elseif codeval==6
    nowval=nowrow.DEP_KM;
  elseif codeval==7
    nowval=nowrow.MAG;
  end
  
  % Check if we should plot this earthquake, if applicable
  if codeval>0
    if plotall==1
      if nowval<quantile(allvals,0.95)
        continue;
      end
    elseif plotall==-1
      if nowval>quantile(allvals,0.05)
        continue;
      end
    end
  end
 
  % Plot Color 
  if codeval>0
    valplace=(nowval-minval)/valrange;
    clrnum=round(valplace*numclr);
    if clrnum==0
      nowclr=eqclrs(1,:);
    else
      nowclr=eqclrs(clrnum,:);
    end
  else
    nowclr=[0.65 0.65 0.65];
  end
  %
  nowplot=scatter(nowlon,nowlat);
  nowplot.MarkerEdgeColor=[0 0 0];
  nowplot.LineWidth=1.5;
  nowplot.MarkerFaceAlpha=0.7;
  if numevents>2000
    nowplot.SizeData=25;
  elseif numevents>1000
    nowplot.SizeData=50;
  elseif numevents>100
    nowplot.SizeData=75;
  else
    nowplot.SizeData=100;
  end
  try
    nowplot.MarkerFaceColor=nowclr;
  catch
    keyboard
  end
  nowplot.LineWidth=1;
  if i==1
    hold on
  end
end

% Color bar if needed!
cbarlbls={'Final Lon - Initial Lon (deg)';...
  'Final Lat - Initial Lat (deg)';'Final Depth - Initial Depth (km)';...
  'Final - Initial Horizontal (deg)';'Final - Initial 3D (km)';...
  'Depth (km)';'Magnitude'};
if codeval>0
  % Adjust position of plot   
  if makefig==1
    mapfig.Position(3)=0.53;
    mapfig.Position(4)=0.77;
    nowax.Position(1)=0.1;
    nowax.Position(3)=0.75;
  end
  % Overlay second axes
  ax2=axes;
  allax=horzcat(allax,ax2);
  ax2.Visible='off';
  % Colormap
  cmap=colormap(ax2,eqclrs);
  % Colorbar
  cbar=colorbar(ax2);
  cbar.Location='eastoutside';
  cbar.Position(1)=0.87;
  cbar.Position(2)=0.268;
  cbar.Position(4)=0.65;
  caxis([minval maxval])
  if codeval==6
    cbar.Ticks=[cbar.Ticks floor(maxval)];
  end
  if codeval==7
    cbar.Ticks=[minval cbar.Ticks];  
  end
  % Label
  cbar.Label.String=cbarlbls{codeval};
  cbar.Label.Rotation=270;
  cbar.Label.Units='normalized';
  cbar.Label.Position(1)=4;
  cbar.FontSize=13;
  if codeval==6
    cbar.Direction='reverse';
  end
end

% Save figure
if maprange==0
  mapstr='tajikistan.all';
else
  mapstr='tajikbasin';
end
%
if addfaults==0
  faultstr='nofault';
elseif addfaults==1 || addfaults==2
  faultstr='faults';
else
  faultstr='majrfaults';
end
%
eqstr=sprintf('%dEVTs',numevents);
%
codestrs={'DeltaLon';'DeltaLat';'DeltaDep';'Delta2D';'Delta3D';'DepCode';...
  'MagCode'};
if ~isempty(savedir)
  if codeval==0
    % No color code
    fname=sprintf('TOTALEQwTOPO_%s_%s_%s.eps',eqstr,...
      mapstr,faultstr);
  else
    if plotall==0
      fname=sprintf('TOTALEQwTOPO_%s_%s_%s_%s.eps',codestrs{codeval},...
        eqstr,mapstr,faultstr);
    elseif plotall==1
      fname=sprintf('TOTALEQwTOPO_%s_95PRC_%s_%s_%s.eps',codestrs{codeval},...
        eqstr,mapstr,faultstr);
    elseif plotall==-1
      fname=sprintf('TOTALEQwTOPO_%s_5PRC_%s_%s_%s.eps',codestrs{codeval},...
        eqstr,mapstr,faultstr);
    end
  end
  setenv('EPS',savedir);
  figname=figdisp(fname,[],[],2,[],'epstopdf');
else
  figname='';
end


