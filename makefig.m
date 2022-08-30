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
% if ispc
%   %fig.Position  = [20 0.415 20 10.03125];
%   fig.Position = [11.4583    2.3229    8.5312    8.0312];
% else
%   fig.Position  = [20 8 15 10];
% end

% Set to 1/3 screen size
% Set screen size based on computer
set(0,'units','inches');%Sets the units of your root object (screen) to inches
Inch_SS = get(0,'screensize');%Obtains this inch information
pos_x      = Inch_SS(1);
pos_y      = Inch_SS(2);
plotwidth  = Inch_SS(3)*1/3;
plotheight = Inch_SS(4);
MP = get(0, 'MonitorPositions');
N = size(MP, 1);
if N == 2
  pos_x = pos_x + MP(2,1);
end
fig.Position = [pos_x pos_y plotwidth plotheight];

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
