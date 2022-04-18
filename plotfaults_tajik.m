function plotfaults_tajik(faulttypes,clrcode)
% 
% Function to plot on a map the faults from Gagala et al. 2020 and 
% Kufner et al. 2018
% 
% INPUTS
% faulttypes : Which faults do we want to plot? Enter as a cell array, where
%              each type of fault is entered as a string entry in the 
%              array : 
% 
%              'normal' : Normal faults
%              'thrust' : Thrust faults
%              'strike' : Strike-slip faults
%          'subsurface' : Subsurface faults 
%           'uncertain' : Uncertain faults
%      'uncertain subs' : Uncertain subsurface faults
%               'major' : Major, named faults
% 
%              Enter one or more of these strings in a cell array to plot
%              one or more types of fault
% 
% clrcode : Do we wish to color-code our plots? 
%           0 : No 
%           1 : Yes [Default]
% 
% Uses defval.m in csdms-contrib/slepian_alpha
% 
% Last Modified: January 16, 2022 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default variables
defval('clrcode',1);

% Get file(s) containing coordinates of inputted fault(s)
[nrmfile,thrfile,strfile,subfile,uncertfile,unsfile,majrfile]=...
  getfaultfiles(faulttypes);

% Normal fault
if exist(nrmfile)==2
  nrmfaults=shaperead(nrmfile);
  for nfault = 1:length(nrmfaults)
    flt_plot=plot(nrmfaults(nfault).X,nrmfaults(nfault).Y);
    flt_plot.LineWidth=1.5;
    if clrcode==1
      flt_plot.Color=[1 0.5 0]; 
    else
      flt_plot.Color=[0 0 0];
    end
    flt_plot.HandleVisibility='off';
    if nfault==1
      hold on
    end
  end
end
      
% Thrust faults 
if exist(thrfile)==2
  thrfaults=shaperead(thrfile);
  for tfault = 1:length(thrfaults)
    flt_plot=plot(thrfaults(tfault).X,thrfaults(tfault).Y);
    flt_plot.LineWidth=1.5;
    if clrcode==1
      flt_plot.Color=[1 0 0]; 
    else
      flt_plot.Color=[0 0 0];
    end
    flt_plot.HandleVisibility='off';
    if tfault==1
      hold on
    end
  end
end
    
% Strike-slip faults
if exist(strfile)==2
  strfaults=shaperead(strfile);
  for ssfault = 1:length(strfaults)
    flt_plot=plot(strfaults(ssfault).X,strfaults(ssfault).Y);
    flt_plot.LineWidth=1.5;
    if clrcode==1
      flt_plot.Color=[0 1 1]; 
    else
      flt_plot.Color=[0 0 0];
    end
    flt_plot.HandleVisibility='off';
    if ssfault==1
      hold on
    end
  end
end
      
% Subsurface faults 
if exist(subfile)==2
  subsfaults=shaperead(subfile);
  for sfault = 1:length(subsfaults)
    flt_plot=plot(subsfaults(sfault).X,subsfaults(sfault).Y);
    flt_plot.LineWidth=1.5;
    if clrcode==1
      flt_plot.Color=[0 0.5 1]; 
    else
      flt_plot.Color=[0 0 0];
    end
    flt_plot.HandleVisibility='off';
    if sfault==1
      hold on
    end
  end
end
    
% Uncertain faults 
if exist(uncertfile)==2
  unfaults=shaperead(uncertfile);
  for ufault = 1:length(unfaults)
    flt_plot=plot(unfaults(ufault).X,unfaults(ufault).Y);
    flt_plot.LineWidth=1.5;
    if clrcode==1
      flt_plot.Color=[0.3 0.3 0.3]; 
    else
      flt_plot.Color=[0 0 0];
    end
    flt_plot.HandleVisibility='off';
    if ufault==1
      hold on
    end
  end
end
      
% Uncertain subsurface 
if exist(unsfile)==2
  unsubfaults=shaperead(unsfile);
  for usfault = 1:length(unsubfaults)
    flt_plot=plot(unsubfaults(usfault).X,unsubfaults(usfault).Y);
    flt_plot.LineWidth=1.5;

    if clrcode==1
      flt_plot.Color=[0.3 0.3 0.3]; 
      flt_plot.LineStyle=':';
    else
      flt_plot.Color=[0 0 0];
    end
    flt_plot.HandleVisibility='off';
    if usfault==1
      hold on
    end
  end
end
    
% Major faults
if exist(majrfile)==2
  faulttbl=readtable(majrfile,'Format','%f%f%f%f%f%f%f%f%f%f%f%f');
  %
  alburz_lon=faulttbl.ALBURZ_LON;
  alburz_lon=alburz_lon(alburz_lon~=-99);
  alburz_lat=faulttbl.ALBURZ_LAT;
  alburz_lat=alburz_lat(alburz_lat~=-99);
  %
  illiac_lon=faulttbl.ILLIAC_LON;
  illiac_lon=illiac_lon(illiac_lon~=-99);
  illiac_lat=faulttbl.ILLIAC_LAT;
  illiac_lat=illiac_lat(illiac_lat~=-99);
  %
  vakhsh_lon=faulttbl.VAKHSH_LON;
  vakhsh_lon=vakhsh_lon(vakhsh_lon~=-99);
  vakhsh_lat=faulttbl.VAKHSH_LAT;
  vakhsh_lat=vakhsh_lat(vakhsh_lat~=-99);
  %
  darvaz_lon=faulttbl.DARVAZ_LON;
  darvaz_lon=darvaz_lon(darvaz_lon~=-99);
  darvaz_lat=faulttbl.DARVAZ_LAT;
  darvaz_lat=darvaz_lat(darvaz_lat~=-99); 
  %
  lsfknf_lon=faulttbl.LSFKNF_LON;
  lsfknf_lon=lsfknf_lon(lsfknf_lon~=-99);
  lsfknf_lat=faulttbl.LSFKNF_LAT;
  lsfknf_lat=lsfknf_lat(lsfknf_lat~=-99);


  % Plot
  if clrcode==0
    faultclr=[0 0 0];
  else
    faultclr=[1 0.65 1];
  end
  alburz_plot=plot(alburz_lon,alburz_lat);
  alburz_plot.LineWidth=2.5;
  alburz_plot.Color=faultclr;
  hold on 
  %
  illiac_plot=plot(illiac_lon,illiac_lat);
  illiac_plot.LineWidth=2.5;
  illiac_plot.Color=faultclr;
  %
  vakhsh_plot=plot(vakhsh_lon,vakhsh_lat);
  vakhsh_plot.LineWidth=2.5;
  vakhsh_plot.Color=faultclr;
  %
  darvaz_plot=plot(darvaz_lon,darvaz_lat);
  darvaz_plot.LineWidth=2.5;
  darvaz_plot.Color=faultclr;
  %
  lsf_lon=lsfknf_lon(10:end,:);
  lsf_lat=lsfknf_lat(10:end,:);
  lsf_plot=plot(lsf_lon,lsf_lat);
  lsf_plot.LineWidth=2.5;
  lsf_plot.Color=faultclr;
  % 
  knf_lon=lsfknf_lon(1:10,:);
  knf_lat=lsfknf_lat(1:10,:);
  knf_plot=plot(knf_lon,knf_lat);
  knf_plot.LineWidth=2.5;
  knf_plot.Color=faultclr;
end

