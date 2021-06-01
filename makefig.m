function fig = makefig
%% function fig = makefig
% make standard figure 
% author: Brita Irving <bkirving@alaska.edu>
%
%uses and resets current figure
% fig = figure(1);
% clf('reset');
% make new figure
axiscolor = [0.9 0.9 0.9];
n = 1;
while ishandle(n)
  n = n + 1;
end
fig = figure(n);
fig.Units = 'inches';
%fig.Position = [20.3021 0.5104 17.0312 9.4479];
if ispc
  fig.Position  = [20 0.415 20 10.03125];
else
  fig.Position  = [20 8 15 10];
end
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'portrait';
% axes defaults
set(fig,'DefaultTextFontname',    'Times');
set(fig,'defaultTextFontSize',    20);
set(fig,'DefaultTextFontName',    'Times');
set(fig,'DefaultAxesFontName',    'Times');
set(fig,'DefaultAxesFontSize',    18);
set(fig,'DefaultAxesFontWeight',  'Bold');
set(fig, 'color', 'w');
set(fig,'DefaultAxesLineStyleOrder','-|--|-.|:');

tomLnCmap = [ 0 0 1; 1 0 0; 1 0.6 0; 0.85 0 0.95; 0 .95 0.95; 0.8 0.8 0; 0.8 0.5 0; 0.5 0 0; 0 0 0.5; 0 0.65 0.65];
%               b;     r;   orange;      ~m;        ~bc;          ~y;    dorange;      dr;     db;       dc;
% no greens as Tom is R-G color blind.
set(fig,'DefaultAxesColorOrder',tomLnCmap);
% set axis color
% set(gca,'color',[0.9 0.9 0.9])
set(gca,'YDir','reverse','XAxisLocation','bottom','Layer','top','Color',axiscolor);
end
% % uses and resets current figure
% % fig = figure(1);
% % clf('reset');
% % make new figure
% n = 1;
% while ishandle(n)
%   n = n + 1;
% end
% fig = figure(n);
% fig.Units = 'inches';
% fig.Position = [20.3021 0.5104 17.0312 9.4479];
% fig.PaperPositionMode = 'auto';
% fig.PaperOrientation = 'portrait';
% % axes defaults
% plot_characteristics_standard
% set(fig, 'color', 'w');
% 
% set(fig,'DefaultAxesLineStyleOrder','-|--|-.|:');
% 
% tomLnCmap = [ 0 0 1; 1 0 0; 1 0.6 0; 0.85 0 0.95; 0 .95 0.95; 0.8 0.8 0; 0.8 0.5 0; 0.5 0 0; 0 0 0.5; 0 0.65 0.65];
% %               b;     r;   orange;      ~m;        ~bc;          ~y;    dorange;      dr;     db;       dc;
% % no greens as Tom is R-G color blind.
% set(fig,'DefaultAxesColorOrder',tomLnCmap);
% % set axis color
% % set(gca,'color',[0.9 0.9 0.9])
