function [data, cfg] = manual_timeSeries_QC(cfg,data,vars,backup_file)
%% FUNCTION manual_timeSeries_QC
%
%  SYNTAX: [data, cfg] = manual_timeSeries_QC(cfg,data,vars,backup_file)
%
%  DESCRIPTION: Clean up time series data visually by inspecting and
%    manually flagging for spikes, stuck sensors, drop outs, etc.
%    Developed to QC timeseries data from sensors deployed on a mooring but
%    could be used for any timeseries data.
%
%  INPUT: [data, cfg] = manual_timeSeries_QC(cfg,data,vars,backup_file)
%    cfg         | structure containing metadata, paths, flags etc.
%    data        | structure containing data arrays. Fields must contain
%                   "datenum", as well as all fields defined in "vars"
%    vars        | cell array with field to inspect and plot. First entry
%                   is the field that you are QC'ing.
%    backup_file | optional filename to load, instead of starting from
%                   scratch.
%
%  EXAMPLES:
%    This script can be used for any datatype, as long as data is a
%    structure with fields "datenum" and the entries in cell array "vars"
%
%    SeapHOx ------------------ pH data -----------------------------------
%     QC ph_cor, and also view pCO2, salinity, and density
%     [data_avg, cfg] = manual_timeSeries_QC(cfg,data_avg,{'ph_cor' 'pco2' 'salinity' 'density'})
%     QC ph_cor, & also view raw ph from external electrode, and density
%     [data_avg, cfg] = manual_timeSeries_QC(cfg,data_avg,{'ph_cor' 'ph_ext''density'})
%
%    SUNA  ------------------ Nitrate data --------------------------------
%     QC NO3_cor, and also view raw nitrate, cdom, and density (burst resolution)
%     [data_avg, cfg] = manual_timeSeries_QC(cfg,data_avg,{'NO3_uM_cor' 'NO3_uM' 'cdom' 'density'})
%     QC NO3_cor, and also view raw nitrate, spectra, and density (full resolution)
%        **add "spectra" to cell array if want to plot spectrogram
%     [data_avg, cfg] = manual_timeSeries_QC(cfg,data_proc,{'NO3_uM_cor' 'NO3_uM' 'spectra' 'density'})
%
%  AUTHORS:
%    Brita K Irving  <bkirving@alaska.edu, https://github.com/britairving/>
%    Adapted from a script by Liz Dobbins (UAF CFOS), which was originally
%    adapted from a script by Seth Danielson (UAF CFOS).
%% 0 | Set defaults and script options
% Run script in debug mode - if runs into an error, script will stop with all variables accessible.
dbstop if error
default_flag  = struct();
%----------------------------   %Code  | Value         | Definition
default_flag.good          = 1; % 1    | Good          | Passed documented QC tests
default_flag.not_evaluated = 2; % 2    | Not evaluated | not available, or unknown | Used for data when no QC test was performed, or the information on quality is not available
default_flag.questionable  = 3; % 3    | Questionable  | Failed non-critical documented metric or subjective test
default_flag.bad           = 4; % 4    | Bad           | Failed critical documented QC test(s) or as assigned by the data provider
default_flag.missing_data  = 9; % 9    | Missing Data  | Used as placeholder when data are missing
% Set basic script options
opt.vars          = vars;
opt.zoom_factor   = 2;      % zoom factor, if set to 1, will not zoom
opt.clr_newflag   = 'g';    % color of newly flagged points
opt.clr_oldflag   = 'r';    % color of previously flagged points
opt.prcnt_lims    = [0 100];% adjust to change ylimits on axes
opt.idx_bad       = [];     % array to store flagged points
opt.auto_flag_bad = 1;      % 1 = sets all flags as bad, 0 = asks each time how to flag
opt.var_vs_S      = 1;      % 1 = plots variable vs salinity as well, 0 = does not plot
opt.spectrogram   = 1;      % 1 = plots spectrogram IF nitrate data,  0 = does not plot
% find total time in data so can plot accordingly
tot_sec = etime(datevec(data.datenum(end)),datevec(data.datenum(1)));
tot_day = tot_sec./(60*60*24); % seconds to days
if tot_day < 5
  %  60 minutes in datenum
  opt.time_int = datenum(0000,0,0,1,0,0);
  opt.dntick = 'mm/dd HH:MM';
else
  % 1 day in datenum
  opt.time_int = datenum(0000,0,1,0,0,0);
  opt.dntick = 'mm/dd';
end

%% 1 | Load file, or initialize necessary variables
if nargin == 4 && exist(backup_file,'file')
  fprintf('Loading QC file: %s\n',backup_file)
  load(backup_file);
else
  %% Initialize filenames, variables, etc..
  % first check that variable is a field in data structure
  if ~isfield(data,opt.vars{1})
    fprintf('QC variable "%s" not a field in the data structure!\n',opt.vars{1})
    error('Must pass through "data" as a structure with key variables as fields')
  else
    % short hand for quick access
    var = opt.vars{1};
  end
  % check that datenum is a field in data structure
  if ~isfield(data,'datenum')
    error('"datenum" is a required field in the data structure!')
  end
  % create filename to save qc work to
  if isfield(cfg.path,'qc_name')
    [~,qc_name,~] = fileparts(cfg.path.qc_name);
    qc_file = fullfile([qc_name '_' datestr(now,'yyyymmddHHMM') '.mat']);
  else
    qc_file = fullfile([opt.vars{1} '_qc_' datestr(now,'yyyymmddHHMM') '.mat']);
  end
  qc_file_tmp = strrep(qc_file,'.mat','_temporary.mat');
  if isfield(cfg,'datadir') && exist(cfg.datadir,'dir')
    qc_file = fullfile(cfg.datadir,qc_file);
    qc_file_tmp = fullfile(cfg.datadir,qc_file_tmp);
  end
  fprintf('Data will be automatically saved to %s (as a backup)\n',qc_file_tmp);
  % store original data so backup available if anything is overwritten
  data_original = data;
  
  %% Make sure all passed variables are fields, clear those that are not
  opt.nvars   = numel(opt.vars);
  rm_vars = [];
  % Only applies to SUNA data where spectrogram is useful for QC
  if ismember('spectra',vars) && opt.spectrogram && isfield(data,'spectrum_channels') && isfield(cfg.SUNA.cal1,'caldata')
    opt.spectrogram = 1;
    opt.wavelength = cfg.SUNA.cal1.caldata.Wavelength;
    opt.spec_limits = [0 60000];
    opt.cmap = cmocean('thermal');
  else
    opt.spectrogram = 0;
  end
  
  for nv = 1:opt.nvars
    if ~isfield(data,opt.vars{nv})
      fprintf('%s not available\n',opt.vars{nv})
      rm_vars = [rm_vars; nv];
    end
  end
  % remove
  if ~isempty(rm_vars)
    opt.vars(rm_vars) = [];
    opt.nvars = numel(opt.vars);
  end
  if opt.var_vs_S && ~isfield(data,'salinity')
    opt.var_vs_S = 0;
  end
  
  %% Define default data quality flags, if none exist
  if ~isfield(cfg,'flag')
    % initialize flag structure
    cfg.flag  = default_flag;
  end
  if isfield(data,'flag')
    data.flag_0 = data.flag;
  else
    data.flag = ones(size(data.datenum)) * cfg.flag.not_evaluated;
  end
end

%% 2 | Create figure and set characteristics
if opt.spectrogram
  % Add spectrogram showing spectral intensities measured by the SUNA
  ax = makefig_subplots(1,opt.nvars+2);
  ax = fliplr(ax); % switch order
  ax_spec1 = ax(end-1);
  ax_spec2 = ax(end);
else
  ax = makefig_subplots(1,opt.nvars);
  ax = fliplr(ax); % switch order
end
% Add variable vs salinity plot - readjust axes positions
if opt.var_vs_S
  for nv = 1:opt.nvars
    ax(nv).Position(3) = ax(nv).Position(3) - ax(nv).Position(3)/4;
  end
  % do some arithmatic magic to get axis asthetically aligned
  ax_S = axes('Position',ax(1).Position);
  ax_S.Position(1) = ax(1).Position(1)+ax(1).Position(3)+0.04; % x start
  ax_S.Position(2) = ax(opt.nvars).Position(2)+0.04;           % y start
  ax_S.Position(3) = ax(1).Position(3)/3;                      % width
  ax_S.Position(4) = ax(1).Position(4)*opt.nvars;              % height
end

%% 3 | Plot
% Highlight already flagged data
idx_alreadyflagged = data.flag > cfg.flag.not_evaluated;
for nv = 1:opt.nvars
  idx_good =  data.(opt.vars{nv}) ~= 0;% ignore data that are exactly zero
  plot(ax(nv),data.datenum(idx_good),data.(opt.vars{nv})(idx_good),'ko-','MarkerSize',3,'DisplayName',strrep(opt.vars{nv},'_','\_'));
  hold(ax(nv),'on'); grid(ax(nv),'on'); ax(nv).YDir = 'normal';
  ax(nv).YLabel.String = strrep(opt.vars{nv},'_','\_');
  datetick(ax(nv),'x','mm/yyyy','keeplimits')
  if nv < opt.nvars
    ax(nv).XTickLabel = [];
  end
  hflag = plot(ax(nv),data.datenum(idx_alreadyflagged),data.(opt.vars{nv})(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
  
end
ax(opt.nvars).XAxis.FontSize = 12;

linkaxes(ax,'x'); % links all xaxis so when change limits on one, it will change limits on all
opt.xmin = min(data.datenum);
opt.xmax = max(data.datenum);
hl = legend(ax(1),'show'); hl.FontSize = 12; hl.Location = 'best';
% Initialize Variable VS Salinity plot
if opt.var_vs_S
  hS = plot(ax_S,data.salinity(idx_good),data.(var)(idx_good),'ko','MarkerSize',3,'DisplayName',strrep(opt.vars{nv},'_','\_'));
  ax_S.Title.String = [strrep(var,'_','\_') ' vs salinity'];
  hold(ax_S,'on'); grid(ax_S,'on'); ax_S.YDir = 'normal';
  hS_bad = plot(ax_S,data.salinity(idx_alreadyflagged),data.(var)(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
  hlS = legend(ax_S,'show'); hlS.FontSize = 12; hlS.Location = 'best';
end
% Plot ISUS spectrogram
if opt.spectrogram
  idx_good =  data.dark ~= 0;% % ignore dark
  
  % Plot 2D spectral intensities at 254nm, 350nm, & average 217-240nm (nitrate fit range)
  i217_to_240 = opt.wavelength >= 217 & opt.wavelength <= 240;
  i254        = opt.wavelength >= 254 & opt.wavelength < 255;
  i350        = opt.wavelength >= 350 & opt.wavelength < 351;
  data.spec_254   = nanmean(data.spectrum_channels(:,i254),2);
  data.spec_350   = nanmean(data.spectrum_channels(:,i350),2);
  data.avg_217240 = nanmean(data.spectrum_channels(:,i217_to_240),2);
  plot(ax_spec1,data.datenum(idx_good),data.avg_217240(idx_good),'ko-','MarkerSize',3,'DisplayName','217-240nm');
  hold(ax_spec1,'on'); grid(ax_spec1,'on'); ax_spec1.YDir = 'normal';
  plot(ax_spec1,data.datenum(idx_good),data.spec_254(idx_good),'k<-','MarkerSize',3,'DisplayName','254nm');
  plot(ax_spec1,data.datenum(idx_good),data.spec_350(idx_good),'k>-','MarkerSize',3,'DisplayName','350nm');
  ax_spec1.YLim =  opt.spec_limits;
  hspec_bad1 = plot(ax_spec1,data.datenum(idx_alreadyflagged),data.avg_217240(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
  hspec_bad2 = plot(ax_spec1,data.datenum(idx_alreadyflagged),data.spec_254(idx_alreadyflagged),  '<','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
  hspec_bad3 = plot(ax_spec1,data.datenum(idx_alreadyflagged),data.spec_350(idx_alreadyflagged),  '>','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
  
  ax_spec1.YLabel.String = 'spectra';
  datetick(ax_spec1,'x','mm/yyyy','keeplimits')
  hlspec = legend(ax_spec1,'show'); hlspec.FontSize = 10; hlspec.Location = 'best';
  ax_spec1.Position(3) = ax(opt.nvars).Position(3);
  ax_spec1.XTickLabel = [];
  
  % Now plot spectrogram
  hp = imagesc(ax_spec2,data.datenum(idx_good),opt.wavelength,data.spectrum_channels(idx_good,:)');
  shading(ax_spec2,'flat');
  ax_spec2.YLabel.String = 'Wavelength'; ax_spec2.YDir = 'normal';
  colormap(ax_spec2,opt.cmap);     cb = colorbar(ax_spec2);
  cb.Limits    = opt.spec_limits; ax_spec2.CLim = cb.Limits;
  hold(ax_spec2,'on'); grid(ax_spec2,'on');
  ax_spec2.Position(3) = ax(opt.nvars).Position(3);
  datetick(ax_spec2,'x','mm/yyyy','keeplimits')
  ax(opt.nvars).XTickLabel = [];
  ax_spec2.XAxis.FontSize  = 12;
end

%% 4 | Begin Manual QC

done = 0;
while ~done
  try
    fprintf('           \n')
    fprintf('  < 0 >   Set Zoom Factor\n')
    fprintf('  < 1 >   Zoom In\n')
    fprintf('  < 2 >   Zoom Out\n')
    fprintf('  < 3 >   Move Left\n')
    fprintf('  < 4 >   Move Right\n')
    fprintf('  < 5 >   Select ([X,Y]) Range of Points\n');
    fprintf('  < 6 >   Select (X) Range of Points\n')
    fprintf('  < 7 >   Make High Kill Threshold\n')
    fprintf('  < 8 >   Make Low Kill Threshold\n')
    if opt.var_vs_S
      fprintf('  < 9 >   Select ([X,Y]) Range of Points on SALINITY plot\n')
    end
    fprintf('  <10 >   Finish and return\n')
    fprintf('  <11 >   Save Plot\n')
    fprintf('  <99 >   Stop [default]\n')
    fprintf('  <199>   Reset Badpoints\n')
    s1 = input('Enter choice: ');
    % default to keyboard if empty or incorrect
    if isstring(s1) || isempty(s1)
      if isstring(s1)
        fprintf('\nIncorrect entry - retry (must be numeric)\n')
      end
      s1 = 99;
    end
    xrange = opt.xmax-opt.xmin;
    %% options for s1 selection
    s1_choice_selection
    idx_good = data.flag <= cfg.flag.not_evaluated;
    trng = data.datenum >= opt.xmin & data.datenum <= opt.xmax;
    %% update plots limits to reflect changes
    ax(1).XLim = [opt.xmin opt.xmax];
    if sum(idx_good) > 0 % continue if no good data within range (or if out of range)
      ax(1).YLim = prctile(data.(var)(trng),opt.prcnt_lims);
    end
    xrange = opt.xmax - opt.xmin;
    opt.time_int = fix(xrange/10);
    if opt.time_int == 0
      opt.time_int = round(xrange/10,2);
    end
    if opt.time_int > 5
      opt.dntick = 'mm/dd';
    elseif opt.time_int <= 5
      opt.dntick = 'mm/dd HH:MM';
    end
    set(ax(opt.nvars),'xlim',[opt.xmin opt.xmax],'xtick',[opt.xmin:opt.time_int:opt.xmax]);
    datetick(ax(opt.nvars),'x',opt.dntick,'keeplimits','keepticks');
    ax(opt.nvars).XTickLabelRotation = 45;
    if opt.spectrogram
      for nv = 1:opt.nvars
        ax(nv).XTick = ax(opt.nvars).XTick;
        ax(nv).XTickLabel = [];
      end
      ax_spec1.XTick      = ax(opt.nvars).XTick;
      ax_spec1.XTickLabel = [];
      set(ax_spec2,'xlim',[opt.xmin opt.xmax],'xtick',[opt.xmin:opt.time_int:opt.xmax]);
      datetick(ax_spec2,'x',opt.dntick,'keeplimits','keepticks');
      ax_spec2.XTickLabelRotation = 45;
    else
      for nv = 1:opt.nvars-1
        ax(nv).XTick = ax(opt.nvars).XTick;
        ax(nv).XTickLabel = [];
      end
    end
    
    %% Show new range on salinity plot
    if opt.var_vs_S
      try delete(hS_rng);  end
      hS_rng = plot(ax_S,data.salinity(trng),data.(var)(trng),'o','MarkerSize',2,'Color',[0.8 0.8 0],'MarkerFaceColor',[0.8 0.8 0],'DisplayName','current range');
      delete(hS_bad);
      hS_bad = plot(ax_S,data.salinity(idx_alreadyflagged),data.(var)(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
    end
    %% update plots to show flagged points
    if ~isempty(opt.idx_bad)
      % ignore flagged points that are exactly zero as this usually indicates
      % different measurement interval
      opt.idx_bad = find(opt.idx_bad);
      opt.idx_bad(data.(var)(opt.idx_bad) == 0) = [];
      % ignore flagged points outside of current range
      t1 = datetime(opt.xmin,'ConvertFrom','datenum');
      t2 = datetime(opt.xmax,'ConvertFrom','datenum');
      opt.idx_bad(~isbetween(datetime(data.datenum(opt.idx_bad),'ConvertFrom','datenum'),t1,t2)) = [];
      % ignore flagged points that are already flagged
      idx_alreadyflagged = find(data.flag >  cfg.flag.not_evaluated);
      opt.idx_bad(ismember(opt.idx_bad,idx_alreadyflagged)) = [];
      %% If any flagged points are left, decide to keep or throw out
      if sum(opt.idx_bad) > 0
        for nv = 1:opt.nvars
          hbad.(['a' num2str(nv)]) = plot(ax(nv),data.datenum(opt.idx_bad),data.(opt.vars{nv})(opt.idx_bad),'o','MarkerSize',3,'Color',opt.clr_newflag,'MarkerFaceColor',opt.clr_newflag,'DisplayName','new flags');
        end
        if opt.var_vs_S
          hS_bad_new = plot(ax_S,data.salinity(opt.idx_bad),data.(var)(opt.idx_bad),'o','MarkerSize',3,'Color',opt.clr_newflag,'MarkerFaceColor',opt.clr_newflag,'DisplayName','new flags');
        end
        if opt.spectrogram
          hspec_bad1new = plot(ax_spec1,data.datenum(opt.idx_bad),data.avg_217240(opt.idx_bad),'o','MarkerSize',3,'Color',opt.clr_newflag,'MarkerFaceColor',opt.clr_newflag,'DisplayName','new flags');
          hspec_bad2new = plot(ax_spec1,data.datenum(opt.idx_bad),data.spec_254(opt.idx_bad),  '<','MarkerSize',3,'Color',opt.clr_newflag,'MarkerFaceColor',opt.clr_newflag,'DisplayName','new flags');
          hspec_bad3new = plot(ax_spec1,data.datenum(opt.idx_bad),data.spec_350(opt.idx_bad),  '>','MarkerSize',3,'Color',opt.clr_newflag,'MarkerFaceColor',opt.clr_newflag,'DisplayName','new flags');
        end
        chc_keep = input('keep new flags <1/0>? ');
        if isempty(chc_keep); chc_keep = 0; end % Default to NO
        if chc_keep
          % Let user decide which flag value to apply
          if ~opt.auto_flag_bad
            fprintf(' ...flag as...\n')
            flags = fieldnames(cfg.flag);
            for nflag = 1:numel(flags)
              if strcmp(flags{nflag},'bad')
                fprintf('   <%d> %s       ** default ** <return>\n',cfg.flag.(flags{nflag}),flags{nflag})
              else
                fprintf('   <%d> %s\n',cfg.flag.(flags{nflag}),flags{nflag})
              end
            end
            flag_chc = input('   enter choice: ');
            if isempty(flag_chc)
              flag_val = cfg.flag.bad;
            else
              flag_val = flag_chc;
            end
          else % No user input, just flag as bad
            flag_val = cfg.flag.bad;
          end
          data.flag(opt.idx_bad) = flag_val;
          delete(hflag);
          idx_alreadyflagged = find(data.flag >  cfg.flag.not_evaluated);
          hflag = plot(ax(1),data.datenum(idx_alreadyflagged),data.(var)(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
          % save badpoints so can load if system crash
          fprintf('saving intermediate QC data to %s\n',qc_file_tmp);
          save(qc_file_tmp, 'data','cfg','opt','-v7.3');
        end
        % Clear
        for nv = 1:opt.nvars
          delete( hbad.(['a' num2str(nv)]) );
        end
        if opt.var_vs_S
          delete(hS_bad_new);
          delete(hS_bad);
          hS_bad = plot(ax_S,data.salinity(idx_alreadyflagged),data.(var)(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
        end
        if opt.spectrogram
          delete(hspec_bad1new); delete(hspec_bad1);
          delete(hspec_bad2new); delete(hspec_bad2);
          delete(hspec_bad3new); delete(hspec_bad3);
          hspec_bad1 = plot(ax_spec1,data.datenum(idx_alreadyflagged),data.avg_217240(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
          hspec_bad2 = plot(ax_spec1,data.datenum(idx_alreadyflagged),data.spec_254(idx_alreadyflagged),  '<','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
          hspec_bad3 = plot(ax_spec1,data.datenum(idx_alreadyflagged),data.spec_350(idx_alreadyflagged),  '>','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
        end
        
        % reset opt.idx_bad
        opt.idx_bad = [];
      end
    end
  catch
    fprintf('something went wrong here...\n')
    fprintf('in keyboard mode... enter "dbcont" to continue\n')
    keyboard
    done = 0;
  end
end % while ~done

fprintf('Chose to exit QC script, exiting now\n');

%% function s1_choice_selection
  function s1_choice_selection
    switch s1
      % ---------------------------------------------------------------------
      case 1   % zoom in
        % this can be very slow when large data set. So, add option if
        % initial zoom just automatically zoom to start at beginning and
        % ending at 1/20th data size.
        if isequal(opt.xmin,min(data.datenum)) && isequal(opt.xmax,max(data.datenum))
          opt.xmin = data.datenum(1);
          iend = round(numel(data.datenum)./20);
          opt.xmax = data.datenum(iend);
        else
          fprintf('Click on zoom endpoints \n')
          ep = ginput(2);
          if isempty(ep)
            zoom(ax(1),opt.zoom_factor);
            opt.xmin = ax(1).XLim(1);
            opt.xmax = ax(1).XLim(2);
          elseif min(size(ep)) > 1  % got enough points
            opt.xmin = min([ep(1,1) ep(2,1)]);
            opt.xmax = max([ep(1,1) ep(2,1)]);
            if opt.xmin >= opt.xmax
              opt.xmin = opt.xmax - 0.5;
            end
          end
          
        end
        % ---------------------------------------------------------------------
      case 2   % zoom out
        xmid = (opt.xmax+opt.xmin)/2;
        opt.xmin = xmid - xrange*opt.zoom_factor;
        opt.xmax = xmid + xrange*opt.zoom_factor;
        % make sure range doesn't go beyond limits
        if opt.xmin < min(data.datenum)
          opt.xmin = min(data.datenum);
        end
        if opt.xmax > max(data.datenum)
          opt.xmax = max(data.datenum);
        end
        % ---------------------------------------------------------------------
      case 3   % move left
        opt.xmin = opt.xmin - xrange*3/4;
        opt.xmax = opt.xmax - xrange*3/4;
        if opt.xmin < min(data.datenum)
          opt.xmin = min(data.datenum);
        end
        if opt.xmax > max(data.datenum)
          opt.xmax = max(data.datenum);
        end
        % ---------------------------------------------------------------------
      case 4   % move right
        opt.xmin = opt.xmin + xrange*3/4;
        opt.xmax = opt.xmax + xrange*3/4;
        if opt.xmin < min(data.datenum)
          opt.xmin = min(data.datenum);
        end
        if opt.xmax > max(data.datenum)
          opt.xmax = max(data.datenum);
        end
        % ---------------------------------------------------------------------
      case 5   % Select ([X,Y]) Range of Points
        fprintf('  Draw rectangle around bad points\n');
        x = getrect(ax(1)); %returned as a 4-element numeric vector with the form [opt.xmin ymin width height].
        t1 = x(1);
        t2 = x(1)+x(3);
        rng = data.datenum >= t1 & data.datenum <= t2 & ...
          data.(var)   >= x(2) & data.(var) <= x(2)+x(4);
        opt.idx_bad = rng;
        % ---------------------------------------------------------------------
      case 6   % Select (X) Range of Points
        fprintf('  Select limits around bad points\n');
        s2 = ginput(2);
        s2 = sortrows(s2,1);
        opt.idx_bad = data.datenum >= s2(1,1) & data.datenum < s2(2,1);
        % ---------------------------------------------------------------------
      case 7   % make high kill threshold
        fprintf('  Select high kill threshold\n');
        s2      = ginput(1);
        opt.idx_bad = data.(var) > s2(2);
        % ---------------------------------------------------------------------
      case 8   % make low kill threshold
        fprintf('  Select low kill threshold\n');
        s2 = ginput(1);
        opt.idx_bad = data.(var) < s2(2);
        % ---------------------------------------------------------------------
      case 9
        if opt.var_vs_S
          fprintf('   ** On %s VS salinity plot **\n',var);
          fprintf('  Draw rectangle around bad points \n');
          x = getrect(ax_S); %returned as a 4-element numeric vector with the form [opt.xmin ymin width height].
          rng = data.salinity >= x(1) & data.salinity <= x(1)+x(3) & ...
            data.(var)    >= x(2) & data.(var)    <= x(2)+x(4);
          opt.idx_bad = rng;
          %  reset to full data view
          opt.xmin = min(data.datenum);
          opt.xmax = max(data.datenum);
        end
        % ---------------------------------------------------------------------
      case 10  % finish and return
        fprintf('  Saving data to %s\n',qc_file);
        save(qc_file,'cfg','data','data_original','-v7.3')
        idx_alreadyflagged = find(data.flag >  cfg.flag.not_evaluated);
        delete(hflag); hflag = plot(ax(1),data.datenum(idx_alreadyflagged),data.(var)(idx_alreadyflagged),'o','MarkerSize',3,'Color',opt.clr_oldflag,'MarkerFaceColor',opt.clr_oldflag,'DisplayName','flagged');
        standard_printfig_highrespng(erase(qc_file,'.mat'));
        done = 1;
        % ---------------------------------------------------------------------
      case 11  % save plot
        standard_printfig_highrespng(erase(qc_file_tmp,'.mat'));
        % ---------------------------------------------------------------------
      case 0   % set zoom factor
        figure(fig)
        opt.zoom_factor = input('Enter Zoom Factor = ? ');
        % ---------------------------------------------------------------------
      case 99  % stop
        % stops execution of the file and gives control to the keyboard
        % K appears before prompt
        disp(' KEYBOARD (debug) MODE:')
        disp(' ...................... To terminate debug mode and continue execution, use the dbcont command.')
        disp(' ...................... To terminate debug mode and exit the file without completing execution, use the dbquit command.')
        keyboard
        % ---------------------------------------------------------------------
      case 199 % clear badpoints
        clear_chc = input('Do you want to clear the original flagged points as well? (y/n)','s');
        if strcmp(clear_chc,'y')
          % add extra check
          clear_chc_check = input('Are you sure?? (y/n)','s');
          if strcmp(clear_chc_check,'y')
            data.flag = ones(size(data.flag)) * cfg.flag.not_evaluated; % set all to not evaluated
          else
            % reset to original flags;
            data.flag = data.flag_0;
          end
          save(qc_file_tmp, 'data','cfg','opt','-v7.3');
        end
        % ---------------------------------------------------------------------
      case isempty(s1)
    end
  end

%% close all the figures before exiting function
close all
end % qaqc_bi_glider_hand_cleanup



