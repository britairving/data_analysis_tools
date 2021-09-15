function cfg = mooring_calibration_discrete_samples(data_avg,meta_proc,cfg)
%FUNCTION mooring_calibration_discrete_samples
%
%  Syntax:
%     [cfg,discrete] = mooring_calibration_discrete_samples(data_avg,meta_proc,cfg)
%
%  Description:
%     Used to select discrete bottle data as calibration points for
%     determining offset and/or drift for various data.
%
%     Density and salinity are shown from both the bottle and profile data
%     and are useful to indicating if they are sampling the same water mass
%     as the moored sensor itself. I.e. if the bottle was tripped nearby
%     in time and distance, but the water mass has much different
%     charactersitics than reported by the sensor, you know you cannot use
%     that bottle data as a calibration point.
%
%  Inputs:
%     data_avg  | structure containing burst resolution processed data
%     meta_proc | structure containing flags
%     cfg       | structure containing dataset information, table of
%                 calibration casts, necessary paths, parameter, etc.
%
%  Outputs:
%     discrete  | structure containing number_of_samples (# of discrete
%                 bottle samples that can be used for correcting offsets
%                 and/or drift) and detailed information about each
%                 discrete sample.
%     cfg       | structure containing dataset information, table of
%                 calibration casts, necessary paths, parameter, etc.
%
%  Dependencies:
%     GSW Toolbox: https://github.com/TEOS-10/GSW-Matlab
%     SW Toolbox: http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm
%       (I know this is redundant, but sw_dens.m is just so much simpler..)
%     m_map Toolbox: https://www.eoas.ubc.ca/~rich/map.html
%
%  Notes:
%     This script is very interactive and requires the user to carefully
%     examine the plots to determine if specific discrete bottle samples
%     will be useful as calibration points for correcting a dataset.
%
%     Code folding is a very useful way to get an idea of what script
%     entails.
%
%     See Mooring_CalibrationSamples_template.xlsx for template on how to
%     store discrete samples.
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0  | Initialize script variables
if isfield(cfg,'save_figures')
  save_fig = cfg.save_figures;
else
  save_fig = false;
end
close all;

% limits calibration casts to within +/- day_int
day_int = 3;
% Catch if CRUISE field read in oddly
if ismember('x__CRUISE',cfg.calcasts.Properties.VariableNames)
  cfg.calcasts.CRUISE = cfg.calcasts.x__CRUISE;
end
%% 1  | Select which variable to use
% Enter new variables here!
if isfield(cfg,'parameter')
  calvar = cfg.parameter;
elseif isfield(data_avg,'ph_int')     % SeapHOx pH
  inst   = 'SeapHOx';
  calvar = 'pH';      % Variable name that matches corresponding variable in discrete spreadsheet
  if isfield(data_avg,'ph_smo')
    smo_var = 'ph_smo';     % smoothed pH
  else
    smo_var = 'ph_ext_smo'; % smoothed pH from external electrode
  end
  var_string = 'pH';
elseif isfield(data_avg,'NO3_uM') % SUNA Nitrate
  inst    = 'SUNA';
  calvar  = 'NO3_uM';
  smo_var = 'NO3_uM_TCSS_smo';
  var_string = 'Nitrate';
  if isfield(data_avg,'NO3_uM_TCSS_cor')
    cor_var = 'NO3_uM_TCSS_cor';
  end
else
  % Enter information for parameter here following the above model
  % Make sure to add the necessary fields to "ctd" structure in the
  % read_ctd_profiles_by_project.m script.
  fprintf('parameter unknown... set up here\n')
  keyboard
end
smo_str = [inst ' ' strrep(smo_var,'_','\_')];

if exist('cor_var','var')
  plot_cor = 1;
  cor_str  = [inst ' ' strrep(cor_var,'_','\_')];
else
  plot_cor = 0;
end


if strcmp(cfg.project,'CEO_2017_2018') && contains(smo_var,'ph','IgnoreCase',true)
  idx_keep = strcmp(cfg.calcasts.CRUISE,'HLY1702') & strcmp(cfg.calcasts.StationID,'CEO');
  cfg.calcasts = cfg.calcasts(idx_keep,:);
end


%% 2  | Calculate distance in time and space (horizontally and vertically)
% Loops through cfg.calcasts table and calculates the distance in time,
% latitude/longitude, and depth from each row with the nearest point in
% data_avg.
cfg.calcasts.distkm = nan(size(cfg.calcasts.CastID));
cfg.calcasts.distdb = nan(size(cfg.calcasts.CastID));
cfg.calcasts.disthr = nan(size(cfg.calcasts.CastID));
for nsamp = 1:size(cfg.calcasts,1)
  % Find distance horizontally and vertically
  cfg.calcasts.distkm(nsamp) = round(m_lldist([cfg.calcasts.Lon(nsamp) cfg.mooring.longitude],[cfg.calcasts.Lat(nsamp) cfg.mooring.latitude]),2);     % km separating measurements
  % Calculate pressure
  if isnan(cfg.calcasts.Pressure(nsamp))
    cfg.calcasts.Pressure(nsamp) = gsw_p_from_z(-cfg.calcasts.DepSM(nsamp),cfg.calcasts.Lat(nsamp));
  end
  % find distance in time
  idx_gd = find(data_avg.flag <= meta_proc.flag.not_evaluated);
  time_dist = round(cfg.calcasts.dnum(nsamp)-data_avg.datenum(idx_gd),2);
  idx_time_near = nearest(time_dist,0); % distance in days
  idx_time_near = idx_gd(idx_time_near);
  cfg.calcasts.disthr(nsamp) = round(abs(etime(datevec(data_avg.datenum(idx_time_near)),datevec(cfg.calcasts.dnum(nsamp))))./60./60,2);
  %deploy_hours  = round(abs(etime(datevec(cfg.mooring.deploydate),datevec(cfg.discrete_ref.dnum(nsamp))))./60./60,2);  % Hours separating measurements (2 decimal places)
  %recover_hours = round(abs(etime(datevec(cfg.mooring.recoverdate),datevec(cfg.discrete_ref.dnum(nsamp))))./60./60,2); % Hours separating measurements from recovery
  % Calculate the distance vertically (dbar), mean of pressure +/- day_int 
  idx_pressure_range = data_avg.datenum >= data_avg.datenum(idx_time_near)-day_int & data_avg.datenum <= data_avg.datenum(idx_time_near)+day_int;
  mooring_pressure = nanmean(data_avg.pressure(idx_pressure_range));
  
  cfg.calcasts.distdb(nsamp) = round(abs(cfg.calcasts.Pressure(nsamp)-mooring_pressure),2);
  
  % Calculate density
  if isnan(cfg.calcasts.density(nsamp))
    cfg.calcasts.density(nsamp) = sw_dens(cfg.calcasts.Salinity(nsamp),cfg.calcasts.Temp(nsamp), cfg.calcasts.Pressure(nsamp));
  end

end

%% 3  | (ONLY PH) Calculate pH for samples where discrete pH not available
% Insitu pH is calculated later, if applicable
if strcmp(calvar,'pH')
  cfg.calcasts = calculate_co2sys_ph('calc',cfg,cfg.calcasts);
end



%% 4  | Limit to Calibration casts that contain relevant variable
cfg.calcasts.parameter = cell(size(cfg.calcasts.CastID));
% Loop through castID
idx_match = logical(zeros(size(cfg.calcasts.CastID)));
% Loop through each sample
for nsamp = 1:size(cfg.calcasts,1)
  cal = cfg.calcasts(nsamp,:);
  fprintf('Looking for nearby discrete samples | %s %s\n',char(cal.CRUISE), char(cal.StationID))
  % match parameter
  if isfinite(cal.(calvar))
    idx_match(nsamp) = true;
    cal.parameter = {calvar};
  else
    % _____________________________________________________________________
    % _____________________________ NITRATE _______________________________
    if strcmp(calvar,'NO3_uM')
      if any(isfinite(cal.NO3)) && ~any(isfinite(cal.(calvar)))
        fprintf('  NO3 umol/kg values found... converting to umol/L!\n')
        if all(isnan(cal.Pressure))
          cal.Pressure = gsw_p_from_z(-1.0*cal.DepSM,cal.Lat);
          dens = sw_dens(cal.Salinity,cal.Temp, cal.Pressure);
        else
          dens = sw_dens(cal.Salinity,cal.Temp, cal.Pressure);
        end
        %end
        cal.(calvar) = round(cal.NO3 .* 1000 ./dens,3); % [umol/L] * 1000[L/m3] * 1/[kg/m3]
        idx_match(nsamp) = true;
        cal.parameter = {'NO3_uM'};
      elseif any(isfinite(cal.NO3_NO2_uM))
        fprintf('  NO3 + NO2 uM values found\n')
        idx_match(nsamp) = true;
        cal.parameter = {'NO3_NO2_uM'};
      end
    end % NITRATE
    % _____________________________________________________________________
    % ____________________________ pH _____________________________________
    if strcmp(calvar,'pH')
      if isfinite(cal.pH_calc)
        idx_match(nsamp) = true;
        cal.parameter = {'pH_calc'};
        
      end
    end
  end
  try
    cfg.calcasts(nsamp,:) = cal;
  catch
    fprintf('stopped here\n')
    keyboard
  end
end % Loop through unique casts
if ~strcmp(cfg.project,'CEO_2016')
  calcasts = cfg.calcasts(idx_match == 1,:);
else
  calcasts = cfg.calcasts(find(idx_match == 1,1):end,:);
end
% calcasts(calcasts.distdb > 7,:) = [];

%% 5  | Load full depth ctd profile
% This is helpful because it shows the profile behavior with depth so
% allows the user to visually identify if there is something suspicious
% with the bottle data.
ctd = read_ctd_profiles_by_project(cfg);

%% 6  | For each cast, pull out Bottle and CTD indices
mtch = struct();
cruise_station = unique(calcasts.cruise_station,'stable');
for ncast = 1:numel(cruise_station)
  cast_name = cruise_station{ncast};
  mtch.(cast_name) = struct();
  mtch.(cast_name).idx_ctd = find(strcmp(ctd.cruise_station,cruise_station{ncast}));
  mtch.(cast_name).idx_bot = find(strcmp(calcasts.cruise_station,cruise_station{ncast}));
  % only include those with good data
  var_name = calcasts.parameter{ mtch.(cast_name).idx_bot(1)};
  mtch.(cast_name).idx_bot = mtch.(cast_name).idx_bot(isfinite(calcasts.(var_name)(mtch.(cast_name).idx_bot)));
end


%% 3 | Interpolate data to the mooring pressure

for ncast = 1:numel(cruise_station)
idx_cast = strcmp(calcasts.cruise_station,cruise_station{ncast});  
  try
    fields = calcasts.Properties.VariableNames;
    fields = fields(find(strcmp(fields,'Niskin'))+1:end);
    % Loop through fields and calculated interpolated value
    for nfield = 1:numel(fields)
      sfield = fields{nfield};
      % skip these variables
      if ~ismember(sfield,{'dnum' 'distkm' 'distdb' 'disthr'})
        if isnumeric(calcasts.(sfield)) && any(isfinite(calcasts.(sfield))) && ~isdatetime(calcasts.(sfield))
          calcasts.([sfield '_interp'])(idx_cast) = interp1(calcasts.Pressure(idx_cast),calcasts.(sfield)(idx_cast),mooring_pressure);
        end
      end
    end
  catch
    keyboard
  end
end


%% 7  | For each cast, pull out moored sensor data
% now find day_int day range for data
rm_casts = [];
for ncast = 1:numel(cruise_station)
  cast_name = cruise_station{ncast};
  if isempty(mtch.(cast_name).idx_ctd)
    fprintf('No CTD data for cruise_station: %s\n',cast_name)
    fprintf('... maybe just need to add in read_ctd_profiles_by_project.m\n');
    rm_casts = [rm_casts {cast_name}];
  else
    trng = find(data_avg.datenum >= ctd.dnum(mtch.(cast_name).idx_ctd(1))-day_int ...
      & data_avg.datenum <= ctd.dnum(mtch.(cast_name).idx_ctd(1))+day_int ...
      & data_avg.flag    <= meta_proc.flag.not_evaluated);
    if isempty(trng)
      rm_casts = [rm_casts {cast_name}];
      continue
    end
    % now fill mtch struture
    mtch.(cast_name).idx_data     = trng;
    mtch.(cast_name).datenum      = nanmean(data_avg.datenum(trng));
    mtch.(cast_name).pressure     = round(nanmean(data_avg.pressure(trng)),2);
    mtch.(cast_name).temperature  = round(nanmean(data_avg.temperature(trng)),2);
    mtch.(cast_name).salinity     = round(nanmean(data_avg.salinity(trng)),2);
    mtch.(cast_name).density      = round(nanmean(data_avg.density(trng)),2);
    mtch.(cast_name).(smo_var)    = round(nanmean(data_avg.(smo_var)(trng)),2);
    if plot_cor
      mtch.(cast_name).(cor_var) = round(nanmean(data_avg.(cor_var)(trng)),2);
    end
  end
end
% remove casts that do not fit within window
if ~isempty(rm_casts)
  mtch  = rmfield(mtch,rm_casts);
  cruise_station(strcmp(cruise_station,rm_casts)) = [];
end

%% 8  | (ONLY PH) Calculate in-situ pH
if strcmp(calvar,'pH')
  % Calculate insitu pH
  calcasts = calculate_co2sys_ph('insitu',cfg,calcasts,mtch);
end

%% 9  | Plot#1 | CTD cast with locations of discrete samples
% Examine plot to decide whether to use discrete sample as a calibration
% point
% Density and salinity are shown from both the bottle and profile data and
% are useful to indicating if they are sampling the same water mass as the
% moored sensor itself. I.e. if the bottle was tripped nearby in time and
% distance, but the water mass has much different charactersitics than
% reported by the sensor, you know you cannot use that bottle data as a
% calibration point.
makefig; ax = subplot(1,5,1:3); ax2 = subplot(1,5,4); ax3 = subplot(1,5,5);
hold(ax,'on');  grid(ax,'on');  set(ax,'YDir','rev');
hold(ax2,'on'); grid(ax2,'on'); set(ax2,'YDir','rev');
hold(ax3,'on'); grid(ax3,'on'); set(ax3,'YDir','rev');
ax.XLabel.String = 'Density [kg/m^3]';
ax.YLabel.String = 'Pressure [dbar]';
ax2.XLabel.String = var_string;
ax2.YLabel.String = 'Pressure [dbar]';
ax3.XLabel.String = 'Salinity';
ax3.YLabel.String = 'Pressure [dbar]';

clrs = jet(numel(cruise_station));
for ncast = 1:numel(cruise_station)
  cast_name = cruise_station{ncast};
  idx_ctd = mtch.(cast_name).idx_ctd;
  idx_bot = mtch.(cast_name).idx_bot;
  ctd_str = strrep([ctd.Station{idx_ctd(1)} ': Cast#' num2str(ctd.Cast(idx_ctd(1)))],'_','\_');
  bot_str = [ctd_str ' ' strrep(calcasts.parameter{idx_bot(1)},'_','\_') ];
  var_name = calcasts.parameter{idx_bot(1)};
  try
    imin = nearest(calcasts.Pressure(idx_bot),mtch.(cast_name).pressure);
  catch
    keyboard
  end
  imin = idx_bot(imin);
  plot_text = [pad(calcasts.CRUISE{imin},10) ' | Station: ' pad(calcasts.StationID{imin},10) ' Cast#'  num2str(calcasts.CastID(imin)) ...
    ' | ' strrep(var_name,'_','\_') ...
    ...%' @ ' num2str(calcasts.Pressure(imin),3) ' dbar' ...
    ' | ' num2str(calcasts.distkm(imin)) 'km away' ...
    ' & ' num2str(calcasts.disthr(imin)) 'hr off'];
  %% Density vs Pressure
  % Plot full cast density
  plot(ax,ctd.density(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor',clrs(ncast,:),'MarkerSize',8,'DisplayName',plot_text);
  % Plot discrete bottle density
  plot(ax,calcasts.density(idx_bot),calcasts.Pressure(idx_bot),'kh','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' ctd_str]);
  % Plot discrete bottle density interpolated to mooring depth
  plot(ax,calcasts.density_interp(idx_bot),calcasts.Pressure_interp(idx_bot),'k<','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' ctd_str '\_interpolated']);
  % plot sensor density
  plot(ax,mtch.(cast_name).density,mtch.(cast_name).pressure,'kd','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',[inst ' Density @' ctd_str]);

  %% Parameter vs Pressure
  plot(ax2,calcasts.(var_name)(idx_bot),calcasts.Pressure(idx_bot),'k-h','MarkerFaceColor',clrs(ncast,:),'MarkerSize',10,'DisplayName',bot_str);
  if ismember([var_name '_interp'],calcasts.Properties.VariableNames)
    plot(ax2,calcasts.([var_name '_interp'])(imin),mtch.(cast_name).pressure,'k<','MarkerFaceColor',clrs(ncast,:),'MarkerSize',12,'DisplayName',[bot_str '\_interpolated']);
  end
  %Plot insitu values too
  %Plot insitu values too
  if strcmp(calvar,'pH')
    imin = nearest(calcasts.Pressure(mtch.(cast_name).idx_bot),mtch.(cast_name).pressure);
    imin = mtch.(cast_name).idx_bot(imin);
    try
      plot(ax2,calcasts.pH_insitu(imin),calcasts.pH_insitu_pressure(imin),'ko','MarkerFaceColor',clrs(ncast,:),'MarkerSize',8,'DisplayName',strrep(bot_str,'_calc','_insitu(CO2SYS)'));
    catch
      keyboard
    end
  end
  
  plot(ax2,mtch.(cast_name).(smo_var),   mtch.(cast_name).pressure,'kd','MarkerSize',10,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',[smo_str ' @' ctd_str]);
  if plot_cor
    plot(ax2,mtch.(cast_name).(cor_var), mtch.(cast_name).pressure,'kd','MarkerSize',10,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',[cor_str ' @ ' ctd_str]);
  end
  %% Salinity vs Pressure
  % Plot full cast density
  plot(ax3,ctd.salinity(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor',clrs(ncast,:),'MarkerSize',8,'DisplayName',plot_text);
  % Plot discrete bottle salinity
  plot(ax3,calcasts.Salinity(idx_bot),calcasts.Pressure(idx_bot),'kh','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' ctd_str]);
  % Plot discrete bottle salinity interpolated to mooring depth
  plot(ax3,calcasts.Salinity_interp(idx_bot),calcasts.Pressure_interp(idx_bot),'k<','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' ctd_str '\_interpolated']);
  % plot sensor salinity
  plot(ax3,mtch.(cast_name).salinity,mtch.(cast_name).pressure,'kd','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',[inst ' Salinity @' ctd_str]);  

end

% % Plot AMBON CEO16 and CEO17 profiles.... no nitrate but just for
% % comparison sake
% if strcmp(cfg.project,'CEO_2016')
%
%   idx_bot = find(strcmp(calcasts.StationID,'CEO16'),1);
%   idx_ctd = find(strcmp(ctd.Station,'CEO16'));
%   plot_text = [pad(calcasts.CRUISE{idx_bot},10) ' | Station: ' pad(calcasts.StationID{idx_bot},10) ' Cast#'  num2str(calcasts.CastID(idx_bot)) ...
%     ' | ' num2str(calcasts.distkm(idx_bot)) 'km away' ...
%     ' & ' num2str(calcasts.disthr(idx_bot)) 'hr off'];
%   plot(ax,ctd.density(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor','k','MarkerSize',8,'DisplayName',plot_text);
%
%   idx_bot = find(strcmp(calcasts.StationID,'CEO17'),1);
%   idx_ctd = find(strcmp(ctd.Station,'CEO17'));
%   plot_text = [pad(calcasts.CRUISE{idx_bot},10) ' | Station: ' pad(calcasts.StationID{idx_bot},10) ' Cast#'  num2str(calcasts.CastID(idx_bot)) ...
%     ' | ' num2str(calcasts.distkm(idx_bot)) 'km away' ...
%     ' & ' num2str(calcasts.disthr(idx_bot)) 'hr off'];
%   plot(ax,ctd.density(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor','g','MarkerSize',8,'DisplayName',plot_text);
% end

axis(ax2,'tight');
ax2.YLim = ax.YLim;
ax2.Box = 'on';
ax.Title.String  = [strrep(cfg.project,'_','\_') ' | Nearby CTD casts and discrete samples'];
hl = legend(ax,'show','Location','northwest');
hl.FontSize =  12;
hl = legend(ax2,'show','Location','best');

hl.FontSize =  12;
hl.Position(1:2) = [0.53 0.69];% 0.12 0.33];
text(ax,0.01,0.02,'Density and salinity shown to highlight water mass characteristics','Units','normalized','FontWeight','bold','FontSize',12);
if save_fig
  savename = fullfile(cfg.datadir,[cfg.project '_discretebottle_comparison1']);
  standard_printfig_highrespng(savename);
end

%% 10 | Plot#2 | Time vs [parameter] with discrete samples
% Examine plot to decide whether to use discrete sample as a calibration
% point
idx_good = data_avg.flag <= meta_proc.flag.not_evaluated;
% makefig; a = subplot(4,1,1:2); a2 = subplot(4,1,3); a3 = subplot(4,1,4);
if strcmp(calvar,'pH') && isfield(data_avg,'pco2')
  ax = makefig_subplots(1,4); ax = fliplr(ax);
  a = ax(1); a2 = ax(2); a3 = ax(3); a4 = ax(4);
  hold(a,'on');  grid(a,'on');  set(a,'YDir','normal');
  hold(a2,'on'); grid(a2,'on'); set(a2,'YDir','normal');
  hold(a3,'on'); grid(a3,'on'); set(a3,'YDir','normal');
  hold(a4,'on'); grid(a4,'on'); set(a4,'YDir','normal');
else
  ax = makefig_subplots(1,3); ax = fliplr(ax);
  a = ax(1); a2 = ax(2); a3 = ax(3);
  hold(a,'on');  grid(a,'on');  set(a,'YDir','normal');
  hold(a2,'on'); grid(a2,'on'); set(a2,'YDir','normal');
  hold(a3,'on'); grid(a3,'on'); set(a3,'YDir','normal');
end
% basic plot to show corrected values
% makefig; a = gca; a.Color = [0.75 0.75 0.75];hold(a,'on'); grid(a,'on'); a.YDir = 'normal';
plot(a,data_avg.datenum(idx_good),data_avg.(smo_var)(idx_good),'ko','MarkerSize',5,'MarkerFaceColor','k','DisplayName',smo_str)
if plot_cor
  plot(a,data_avg.datenum(idx_good),data_avg.(cor_var)(idx_good),'kd','MarkerSize',5,'MarkerFaceColor','m','DisplayName',cor_str)
end
plot(a2,data_avg.datenum(idx_good),data_avg.density(idx_good),'ko','MarkerSize',5,'MarkerFaceColor','k','DisplayName','SeapHOx Density')
plot(a3,data_avg.datenum(idx_good),data_avg.salinity(idx_good),'ko','MarkerSize',5,'MarkerFaceColor','k','DisplayName','SeapHOx Salinity')
if strcmp(calvar,'pH') && isfield(data_avg,'pco2')
  plot(a4,data_avg.datenum(idx_good),data_avg.pco2(idx_good),'ko','MarkerSize',5,'MarkerFaceColor','k','DisplayName','HydroC pCO_2')
end
% Plot discrete samples
for nclosest = 1:numel(cruise_station)
  cast_name = cruise_station{nclosest};
  imin = nearest(calcasts.Pressure(mtch.(cast_name).idx_bot),mtch.(cast_name).pressure);
  imin = mtch.(cast_name).idx_bot(imin);
  var_name = calcasts.parameter{imin};
  plot_text = [pad(calcasts.CRUISE{imin},12) ' | Station: ' pad(calcasts.StationID{imin},10) ' Cast#'  num2str(calcasts.CastID(imin)) ...
    ' | ' strrep(var_name,'_','\_') ...
    ...%' @ ' num2str(calcasts.Pressure(imin),3) ' dbar' ...
    ' | ' num2str(calcasts.distkm(imin)) 'km away' ...
    ' & ' num2str(calcasts.disthr(imin)) 'hr off'];
  plot(a,calcasts.dnum(imin),calcasts.(var_name)(imin),'kh','MarkerSize',16,'MarkerFaceColor',clrs(nclosest,:),'LineWidth',1,'DisplayName',plot_text);
  plot(a2,calcasts.dnum(imin),calcasts.density(imin),'kh','MarkerSize',16,'MarkerFaceColor',clrs(nclosest,:),'LineWidth',1,'DisplayName',plot_text);
  plot(a3,calcasts.dnum(imin),calcasts.Salinity(imin),'kh','MarkerSize',16,'MarkerFaceColor',clrs(nclosest,:),'LineWidth',1,'DisplayName',plot_text);
  % Plot interpolated values too
  if ismember([var_name '_interp'],calcasts.Properties.VariableNames)
    plot(a,calcasts.dnum(imin),calcasts.([var_name '_interp'])(imin),'k<','MarkerFaceColor',clrs(ncast,:),'MarkerSize',8,'DisplayName',[bot_str '\_interpolated']);
    plot(a2,calcasts.dnum(imin),calcasts.density_interp(imin),'k<','MarkerSize',8,'MarkerFaceColor',clrs(nclosest,:),'LineWidth',1,'DisplayName',[plot_text '\_interpolated']);
    plot(a3,calcasts.dnum(imin),calcasts.Salinity_interp(imin),'k<','MarkerSize',8,'MarkerFaceColor',clrs(nclosest,:),'LineWidth',1,'DisplayName',[plot_text '\_interpolated']);
  end
  %Plot insitu values too
  if strcmp(calvar,'pH')
    plot(a,calcasts.dnum(imin),calcasts.pH_insitu(imin),'ko','MarkerFaceColor',clrs(nclosest,:),'MarkerSize',8,'DisplayName',strrep(plot_text,'_calc','_insitu(CO2SYS)'));
    if isfield(data_avg,'pco2')
      plot(a4,calcasts.dnum(imin),calcasts.pCO2_insitu(imin),'ko','MarkerFaceColor',clrs(nclosest,:),'MarkerSize',8,'DisplayName',strrep(plot_text,'pH_calc','pCO2_insitu(CO2SYS)'));
      a4.YLabel.String = 'pCO_2';
      plot(a4,mtch.(cast_name).datenum,calcasts.pCO2_calc_interp(imin),'k<','MarkerFaceColor',clrs(ncast,:),'MarkerSize',8,'DisplayName',[bot_str '\_interpolated']);
      
    end
  end
end
ylabel(a,var_string); ylim(a,[-10 60]);
axis(a,'tight');
datetick(a,'x','keeplimits'); a.XTickLabel = [];
linkaxes([a a2 a3],'x');
datetick(a2,'x','keeplimits'); a2.XTickLabel = [];
datetick(a3,'x','keeplimits');
a.Title.String  = [strrep(cfg.project,'_','\_') ' | ' inst];
hl = legend(a,'show','Location','north');
hl.FontSize = 10;
a2.YLabel.String = 'Density';
a3.YLabel.String = 'Salinity';
if save_fig
  savename = fullfile(cfg.datadir,[cfg.project '_discretebottle_comparison2']);
  standard_printfig_highrespng(savename);
end

%% 11 | Plot#3 | Density and Salinity regressions
% Examine plot to decide whether to use discrete sample as a calibration
% point

% makefig; a = gca; a.Color = [0.75 0.75 0.75];hold(a,'on'); grid(a,'on'); a.YDir = 'normal';
makefig; a1 = subplot(2,1,1);a1.Color = [0.75 0.75 0.75];hold(a1,'on'); grid(a1,'on'); a1.YDir = 'normal';
a2 = subplot(2,1,2); a2.Color = [0.75 0.75 0.75];hold(a2,'on'); grid(a2,'on'); a2.YDir = 'normal';

% if strcmp(cfg.project,'CEO_2017')
%   casts = casts(2:end);
%   clrs  = clrs(2:end,:);
% end
ypos = 0.96;
for ncast = 1:numel(cruise_station)
  % cast indicies and legend labels
  cast_name = cruise_station{ncast};
  idx_data  = mtch.(cast_name).idx_data;
  idx_bot = mtch.(cast_name).idx_bot;
  idx_bot = idx_bot(calcasts.Pressure(idx_bot) > 25);  %Only look at relationship in lower water column
  var_name = calcasts.parameter{idx_bot(1)};
  bot_str = [strrep(cast_name,'_','\_') ' ' strrep(var_name,'_','\_')];
  
  %  smoothed data
  plot(a1,data_avg.density(idx_data),data_avg.(smo_var)(idx_data),'ko','MarkerSize',8,'MarkerFaceColor',clrs(ncast,:), 'DisplayName',smo_str)
  [p,S] = polyfit(data_avg.density(idx_data),data_avg.(smo_var)(idx_data),1); varfit = polyval(p,data_avg.density(idx_data),S);
  plot(a1,data_avg.density(idx_data),varfit,'k-','LineWidth',2,'DisplayName',[inst ': Density NO_{3 Raw} regression'])
  %text(a1,0.02,ypos,[inst ' ' strrep(smo_var,'_','\_') ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  % corrected data
  if plot_cor
    plot(a1,data_avg.density(idx_data),data_avg.(cor_var)(idx_data),'kd','MarkerSize',8,'MarkerFaceColor',clrs(ncast,:),'DisplayName',cor_str)
    [p,S] = polyfit(data_avg.density(idx_data),data_avg.(cor_var)(idx_data),1); varfit = polyval(p,data_avg.density(idx_data),S);
    plot(a1,data_avg.density(idx_data),varfit,'m-','LineWidth',2,'DisplayName',[inst ': Density NO_{3 TCSS} regression'])
    text(a1,0.02,ypos-0.05,[inst ' ' strrep(cor_var,'_','\_') ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  end
  
  % Discrete Samples
  plot(a1,calcasts.density(idx_bot),calcasts.(var_name)(idx_bot),'kh','MarkerSize',15,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' bot_str]);
  % Discrete sample linear fit
  [p,S] = polyfit(calcasts.density(idx_bot),calcasts.(var_name)(idx_bot),1);  varfit = polyval(p,calcasts.density(idx_bot),S);
  plot(a1,calcasts.density(idx_bot),varfit,'--','MarkerSize',5,'Color',clrs(ncast,:),'LineWidth',2,'DisplayName',['Discrete samples: Density ' strrep(var_string,'_','\_') ' regression'])
  %text(a1,0.02,ypos-0.15,[strrep(cast_name,'_','\_') ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  %ypos = ypos - numel(cruise_station)*0.1;
  
  %% Salinity
  % smoothed data
  plot(a2,data_avg.salinity(idx_data),data_avg.(smo_var)(idx_data),'ko','MarkerSize',8,'MarkerFaceColor',clrs(ncast,:), 'DisplayName',smo_str)
  [p,S] = polyfit(data_avg.salinity(idx_data),data_avg.(smo_var)(idx_data),1); varfit = polyval(p,data_avg.salinity(idx_data),S);
  plot(a2,data_avg.salinity(idx_data),varfit,'k-','LineWidth',2,'DisplayName',[inst ': salinity ' var_string ' regression'])
  %text(a2,0.02,ypos,[inst ' ' strrep(smo_var,'_','\_') ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  % corrected data
  if plot_cor
    plot(a2,data_avg.salinity(idx_data),data_avg.(cor_var)(idx_data),'kd','MarkerSize',8,'MarkerFaceColor',clrs(ncast,:),'DisplayName',cor_str)
    [p,S] = polyfit(data_avg.salinity(idx_data),data_avg.(cor_var)(idx_data),1); varfit = polyval(p,data_avg.salinity(idx_data),S);
    plot(a2,data_avg.salinity(idx_data),varfit,'m-','LineWidth',2,'DisplayName',[inst ': salinity NO_{3 TCSS} regression'])
    text(a2,0.02,ypos-0.05,[inst ' ' strrep(cor_var,'_','\_') ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  end
  
  % Discrete Samples
  plot(a2,calcasts.Salinity(idx_bot),calcasts.(var_name)(idx_bot),'kh','MarkerSize',15,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' bot_str]);
  % Discrete sample linear fit
  [p,S] = polyfit(calcasts.Salinity(idx_bot),calcasts.(var_name)(idx_bot),1);  varfit = polyval(p,calcasts.Salinity(idx_bot),S);
  plot(a2,calcasts.Salinity(idx_bot),varfit,'--','MarkerSize',5,'Color',clrs(ncast,:),'LineWidth',2,'DisplayName',['Discrete samples: Salinity ' strrep(calvar,'_','\_') ' regression'])
  %text(a2,0.02,ypos-0.15,[strrep(cast_name,'_','\_') ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  %ypos = ypos - numel(cruise_station)*0.1;
  
  % plot insitu calculated ph
  if strcmp(calvar,'pH')
    imin = nearest(calcasts.Pressure(mtch.(cast_name).idx_bot),mtch.(cast_name).pressure);
    imin = mtch.(cast_name).idx_bot(imin);
    try
      plot(a2,calcasts.pH_insitu_salinity(imin),calcasts.pH_insitu(imin),'ko','MarkerFaceColor',clrs(ncast,:),'MarkerSize',10,'LineWidth',2,'DisplayName',strrep(bot_str,'_calc','_insitu(CO2SYS)'));
      
    catch
      keyboard
    end
  end
end
a1.XLabel.String = 'Density [kg/m^3]';
a1.YLabel.String =  var_string;
a1.Title.String  = [strrep(cfg.project,'_','\_') ' | Density - ' strrep(calvar,'_','\_') ' linear regression'];
hl = legend(a1,'show');
hl.FontSize = 12;
hl.Location = 'bestoutside';

a2.XLabel.String = 'Salinity';
a2.YLabel.String =  var_string;
a2.Title.String  = ['Salinity - '  strrep(calvar,'_','\_') ' linear regression'];
hl = legend(a2,'show');
hl.FontSize = 12;
hl.Location = 'bestoutside';
if save_fig
  savename = fullfile(cfg.datadir,[cfg.project '_discretebottle_comparison3']);
  standard_printfig_highrespng(savename);
end

%% 12 | FINALLY - CHOOSE WHAT TO DO
% After examining the plots to determine if any discrete samples will be
% useful as calibration points - prompts user to decide which/how to store
% the individual calibration points.
discrete = struct(); % structure that will hold all discrete cal point information
ncal     = 0;        % discrete calibraiton point counter
DONE_CHOOSING = 0;   % switch to exit loop
while ~DONE_CHOOSING
  first_select_done = 0;
  while ~first_select_done
    fprintf('\n')
    fprintf('__________________________________________________________________\n')
    fprintf('Do you want to use discrete sample(s) for processing this %s data?\n',inst)
    fprintf('  <0>  No\n')
    fprintf('  <1>  Yes\n')
    fprintf('  <9>  STOP\n')
    chc = input('  Enter choice: ');
    if isempty(chc) || chc == 0 % Default to NO
      first_select_done = 1;
      DONE_CHOOSING     = 1;
    elseif chc == 9
      fprintf('in keyboard mode... enter "dbcont" to continue\n')
      keyboard
      first_select_done = 0;
    elseif chc == 1
      first_select_done = 1;
    end
  end
  if chc == 1 % Select discrete samples to use as calibration points
    DONE_CHOOSING_CALPOINTS = 0;
    while ~DONE_CHOOSING_CALPOINTS
      %% Select cruise
      fprintf('\nSelect which cruise & cast\n')
      fprintf('  ___________________________________________________________________________\n')
      fprintf('       Cruise       Station       Distance(km)  Distance(hr)\n')
      fprintf('  ---------------------------------------------------------------------------\n')
      num = 0;
      
      for ncast = 1:numel(cruise_station)
        cast_name = cruise_station{ncast};
        idx_bot = mtch.(cast_name).idx_bot;
        num = num + 1;
        fprintf('  <%d>  %s\t%s \t %.1fkm away \t %.1fhr away\n',num,pad(calcasts.CRUISE{idx_bot(1)},10),pad(calcasts.StationID{idx_bot(1)},10),calcasts.distkm(idx_bot(1)), calcasts.disthr(idx_bot(1)));
      end
      fprintf('  <%d>  None - finished selecting discrete calibration points\n',num+1);
      
      chc_cruise = input('  Enter which cruise: ');
      if chc_cruise == num+1
        DONE_CHOOSING_CALPOINTS = 1;
        DONE_CHOOSING           = 1;
        continue
      else
        % Update discrete calibraiton point counter
        ncal = ncal + 1;
      end
      cast_name = cruise_station{chc_cruise};
      idx_bot = mtch.(cast_name).idx_bot;
      %% Select how to handle discrete calibration points (e.g. average, interpolate, single)
      fprintf('  ________________________________________________________________\n')
      fprintf('  Discrete calibration points from cruise_station: %s\n',cast_name)
      for nsamp = 1:numel(idx_bot)
        if strcmp(calvar,'pH')
          var = calcasts.parameter{idx_bot(1)};
          fprintf('    %s\t%d  %.1f \t %.1fkm away \t %.1fhr away \t %.4f \t %.4f@%.1fdbar\n',pad(calcasts.StationID{idx_bot(nsamp)},10),...
            calcasts.CastID(idx_bot(nsamp)), calcasts.Pressure(idx_bot(nsamp)), ...
            calcasts.distkm(idx_bot(nsamp)), calcasts.disthr(idx_bot(nsamp)),...
            calcasts.(var)(idx_bot(nsamp)),  calcasts.pH_insitu(idx_bot(nsamp)), calcasts.pH_insitu_pressure(idx_bot(nsamp)));
        else
          fprintf('    %s\t%d  %.1f \t %.1fkm away \t %.1fhr away\n',calcasts.StationID{idx_bot(nsamp)},...
            calcasts.CastID(idx_bot(nsamp)), calcasts.Pressure(idx_bot(nsamp)), ...
            calcasts.distkm(idx_bot(nsamp)), calcasts.disthr(idx_bot(nsamp)));
        end
      end
      done_method = 0;
      while ~done_method
        fprintf('\n')
        fprintf('How do you want to store the discrete calibration point?\n')
        fprintf('  <S>  SINGLE point\n')
        fprintf('  <A>  AVERAGE points\n')
        fprintf('  <I>  INTERPOLATE points\n')
        if strcmp(calvar,'pH') % don't give the option of interpolating because will be caluclated at in-situ temperature/salinity later
          fprintf('  <P>  INSITU pH\n')
        end
        fprintf('  <99>  STOP\n')
        chc_method = input('  Enter choice: ','s');
        chc_method = upper(chc_method);
        if ismember(chc_method,{'A' 'S' 'I' 'P'})
          done_method = 1;
        elseif strcmp(chc_method,'99')
          keyboard
          done_method = 0;
        else
          fprintf('  .... you input "%s"... not an option, try again\n',chc_method)
          done_method = 0;
        end
      end %% while ~done_method
      switch chc_method
        case 'S'
          method = 'single';
          fprintf('\n  Select discrete sample point \n')
        case 'I'
          method = 'interpolate';
          fprintf('\n  Select discrete sample points to use for interpolation \n')
        case 'A' 
          method = 'average';
          fprintf('\n  Select discrete sample points to use for average \n')
        case 'P' % INSITU pH
          method = 'insitu';
          fprintf('\n  Select discrete sample point\n')
      end
      % print to screen
      if strcmp(calvar,'pH')
        fprintf('          Cruise        Station   Cast Pressure  Distance(km)    Distance(hr)    pH          pH_insitu\n');
      else
        fprintf('          Cruise        Station   Cast Pressure  Distance(km)    Distance(hr)\n')
      end
      for nsamp = 1:numel(idx_bot)
        if strcmp(calvar,'pH')
          var = calcasts.parameter{idx_bot(1)};
          fprintf('    <N%d>  %s\t%s\t%d  %.1f \t %.1fkm away \t %.1fhr away \t %.4f \t %.4f@%.1fdbar\n',nsamp,pad(calcasts.CRUISE{idx_bot(nsamp)},10),pad(calcasts.StationID{idx_bot(nsamp)},10),...
            calcasts.CastID(idx_bot(nsamp)), calcasts.Pressure(idx_bot(nsamp)), ...
            calcasts.distkm(idx_bot(nsamp)), calcasts.disthr(idx_bot(nsamp)),...
            calcasts.(var)(idx_bot(nsamp)),  calcasts.pH_insitu(idx_bot(nsamp)), calcasts.pH_insitu_pressure(idx_bot(nsamp)));
        else
          fprintf('    <N%d>  %s\t%s\t%d  %.1f \t %.1fkm away \t %.1fhr away\n',nsamp,calcasts.CRUISE{idx_bot(nsamp)},calcasts.StationID{idx_bot(nsamp)},...
            calcasts.CastID(idx_bot(nsamp)), calcasts.Pressure(idx_bot(nsamp)), ...
            calcasts.distkm(idx_bot(nsamp)), calcasts.disthr(idx_bot(nsamp)));
        end
      end
      % use full method string for clarity (could use chc_method still)
      switch method
        case {'single' 'insitu'}
          fprintf('\n  Select discrete sample point \n')
          chc_sample = input('  Enter sample choice N#: ');
          discrete_samples = calcasts(idx_bot(chc_sample),:);
        case {'interpolate' 'average'}
          range_choice_done = 0;
          while ~range_choice_done
            chc_samples1 = input('  Enter first N#: ');
            chc_samples2 = input('  Enter last  N#: ');
            if ~isnumeric(chc_samples1) || ~isnumeric(chc_samples1) || isequal(chc_samples1,chc_samples2)
              range_choice_done = 0;
              fprintf('incorrect indice choices, try again\n')
              continue
            end
            try
              if chc_samples2 < chc_samples1
                discrete_samples = calcasts(idx_bot(chc_samples2:chc_samples1),:);
              else
                discrete_samples = calcasts(idx_bot(chc_samples1:chc_samples2),:);
              end
              range_choice_done = 1;
            catch
              fprintf('indice choices did not work, try again\n')
              range_choice_done = 0;
              continue
            end
          end
      end %% CHC_METHOD
      discrete.number_of_samples = ncal;
      
      %% Fill discrete structure with discrete calibration point information
      scal = ['cal' num2str(ncal)];
      discrete.(scal) = struct();
      discrete.(scal).info   = discrete_samples; % table with all cast information
      discrete.(scal).method = method;           % method used to calculate cal point
      var = discrete_samples.parameter{1};
      if strcmp(var,'pH') || strcmp(var,'pH_calc') && strcmp(method,'insitu')
        var = 'pH_insitu';
        pvar = 'pH_insitu_pressure';
        svar = 'pH_insitu_salinity';
        tvar = 'pH_insitu_temperature';
        uvar = 'pH_u';
      elseif strcmp(method,'interpolate')
        var =  [var '_interp'];
        pvar = 'Pressure_interp';
        svar = 'Salinity_interp';
        tvar = 'Temp_interp';
      else
        pvar = 'Pressure';
        svar = 'Salinity';
        tvar = 'Temp';
      end
      discrete.(scal).cal_var = var;
      discrete.(scal).date        = datestr(nanmean(discrete_samples.dnum));
      discrete.(scal).datenum     = nanmean(discrete_samples.dnum);
      discrete.(scal).pressure    = nanmean(discrete_samples.(pvar));
      discrete.(scal).salinity    = nanmean(discrete_samples.(svar));
      discrete.(scal).temperature = nanmean(discrete_samples.(tvar));
      strtitle = [discrete_samples.CRUISE{1} ' ' datestr(nanmean(discrete_samples.dnum),'dd-mmm-yyyy') ' Station:' discrete_samples.StationID{1} ' Cast:' num2str(discrete_samples.CastID(1))];
      shortlab = [discrete_samples.CRUISE{1} '|' datestr(nanmean(discrete_samples.dnum),'dd-mmm-yyyy') '|' discrete_samples.StationID{1} '|cast' num2str(discrete_samples.CastID(1))];
      switch method
        case 'average'
          discrete.(scal).(var)    = round(nanmean(discrete_samples.(var)),5);
          discrete.(scal).header   = [strtitle ' Niskins:' num2str(discrete_samples.Niskin') ' averaged'];
          discrete.(scal).label1   = [discrete.(scal).header ' ' num2str(discrete_samples.distkm(1),'%.1f') 'km ' num2str(discrete_samples.disthr(1),'%.1f') 'hr'];
        case 'interpolate'
          discrete.(scal).pressure = nanmean(mtch.(cast_name).pressure);
          discrete.(scal).(var)    = round(interp1(discrete_samples.Pressure,discrete_samples.(var),mtch.(cast_name).pressure),5);
          discrete.(scal).header   = [strtitle ' Niskins:' num2str(discrete_samples.Niskin') ' interpolated to ' inst ' pressure'];
          discrete.(scal).label1   = [discrete.(scal).header ' ' num2str(discrete_samples.distkm(1),'%.1f') 'km ' num2str(discrete_samples.disthr(1),'%.1f') 'hr'];
        case {'single' 'insitu'}
          discrete.(scal).(var)    = round(discrete_samples.(var),5);
          discrete.(scal).header   = [strtitle ' Niskin:' num2str(discrete_samples.Niskin)];
          discrete.(scal).label1   = [discrete.(scal).header ' ' num2str(discrete_samples.distkm,'%.1f') 'km ' num2str(discrete_samples.disthr,'%.1f') 'hr ' num2str(discrete_samples.distdb,'%.1f') 'dbar'];
      end
      discrete.(scal).label2 = shortlab;
    end %% while ~DONE_CHOOSING_CALPOINTS
  end %% YES/NO SELECT CALPOINTS
end %% while ~DONE_CHOOSING (whether or not to select calpoints

%% 13 | Store discrete data points and save data
%% 13 | Store discrete data points and save data
if isfield(cfg.path,'calfile_step3')
  discrete_filename = cfg.path.calfile_step3;
else
  discrete_filename = fullfile(cfg.datadir,[cfg.project '_discrete_selected.mat']);
end
save(discrete_filename,'discrete');
cfg.cal_discrete = discrete;
end %% MAIN FUNCTION

