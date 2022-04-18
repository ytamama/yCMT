function brdrlines=plotbrdr_tajik(brdrclr)
% 
% Function to plot the borders of Tajikistan and its neighboring countries
% on a map. 
% 
% INPUT
% brdrclr : Color with which to plot the national borders
%           Example: [0.35 0.35 0.35];
% 
% OUTPUTS
% brdrlines : The borders of Tajikistan and its neighboring countries, in 
%             an array and in the following order:
% 
%             Tajikistan, China, Kyrgyzstan, Kazakhstan, Uzbekistan, 
%             Turkmenistan, Afghanistan, Pakistan, India
% 
% Uses the Climate Data Toolbox by Chad Greene
% 
% Last Modified: September 24, 2021 by Yuri Tamama
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot countries, using borders function by Chad Greene
% Tajikistan
plot_tajikistan=borders('Tajikistan');
plot_tajikistan.Color=brdrclr;
plot_tajikistan.LineWidth=0.25;
hold on
% China
plot_china=borders('China');
plot_china.Color=brdrclr;
plot_china.LineWidth=0.25;
% Kyrgyzstan
plot_kyrgyzstan=borders('Kyrgyzstan');
plot_kyrgyzstan.Color=brdrclr;
plot_kyrgyzstan.LineWidth=0.25;
% Kazakhstan
plot_kazakhstan=borders('Kazakhstan');
plot_kazakhstan.Color=brdrclr;
plot_kazakhstan.LineWidth=0.25;
% Uzbekistan
plot_uzbekistan=borders('Uzbekistan');
plot_uzbekistan.Color=brdrclr;
plot_uzbekistan.LineWidth=0.25;
% Turkmenistan
plot_turkmenistan=borders('Turkmenistan');
plot_turkmenistan.Color=brdrclr;
plot_turkmenistan.LineWidth=0.25;
% Afghanistan
plot_afghanistan=borders('Afghanistan');
plot_afghanistan.Color=brdrclr;
plot_afghanistan.LineWidth=0.25;
% Pakistan
plot_pakistan=borders('Pakistan');
plot_pakistan.Color=brdrclr;
plot_pakistan.LineWidth=0.25;
% India
plot_india=borders('India');
plot_india.Color=brdrclr;
plot_india.LineWidth=0.25;
  

brdrlines=[plot_tajikistan; plot_china; plot_kyrgyzstan;... 
  plot_kazakhstan; plot_uzbekistan; plot_turkmenistan;... 
  plot_afghanistan; plot_pakistan; plot_india];
 
  
  