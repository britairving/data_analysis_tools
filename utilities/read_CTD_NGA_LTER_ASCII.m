function ctd = read_CTD_NGA_LTER_ASCII(ctd_directory)
%% function read_CTD_NGA_LTER_ASCII
%
%
%
%  Notes:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
files = dir(fullfile(ctd_directory,'*.ascii'));
files(contains({files.name},'._')) = [];
if numel(files) > 1
  fprintf('Select which file to use\n');
  for n = 1:numel(files)
    fprintf('  <%d> %s\n',n,files(n).name)
  end
  chc = input('  Enter choice: ');
  files = files(chc);
end

%% Read header file 
% Contains metadata about each cast
hdr_file = fullfile(ctd_directory, strrep(files.name,'.ascii','.hdr'));
if exist(hdr_file,'file')
  hdr = readtable(hdr_file,'FileType','text','CommentStyle','%');
end

%% Read header from ascii file
% Read this to get column names
fid = fopen(fullfile(ctd_directory,files.name),'r');
cols = {};
desc = {};
done = 0;
nline = 0;
while ~done
  
  str = fgetl(fid);
  nline = nline + 1;
  if contains(str,'Data File Column Contents','IgnoreCase',true)
    fprintf('Reading column names\n')
    str = fgetl(fid);
    nline = nline + 1;
    while ~contains(str,'%%%%%%%%%%')
      def  = strsplit(str,':');
      def  = strtrim(def{2});
      def  = strrep(def,' ','_');
      def  = strrep(def,'/','_');
      def  = strrep(def,'-','_');
      def  = strrep(def,' ','_');
      cols = [cols def];
      desc = [desc; str]; % Store full description
      str = fgetl(fid);
      nline = nline + 1;
    end
  end
  if contains(str,'*END*')
    dataline = nline+1;
    done = 1;
  end
end
fclose(fid);

%% Read data from ascii file
ctd = readtable(fullfile(ctd_directory,files.name),'FileType','text','CommentStyle',{'%'});
ctd.Properties.VariableNames

return
  %% CTD files
  keyboard
  % Read cnv L0 downcasts
  opt.path.ctd_file = fullfile(opt.path.ctddir,'NGA_WSD201807_ctd_L0_v1','downcastupcast.mat');
  %     ctd_dir_downcasts   = fullfile(opt.path.ctddir,'NGA_WSD201807_ctd_L0_v1','downcasts');
  %     ctd = cnv2mat_bi(ctd_dir_downcasts);
  %     vars = {'timeJ' 'temp' 't090c' 'sal00' 'sal11' 'par' 'latitude' 'longitude' 'density00' 'depSM' 'prDM' 'flag'};
  %     % initialize variables in ctd structure
  %     for nv = 1:numel(vars)
  %       ctd.(vars{nv}) = [];
  %     end
  %     ctd.cast_num = [];
  %     cast_names = fieldnames(ctd.data);
  %     % combine all data into single variable
  %     for nf = 1:numel(cast_names)
  %       field = cast_names{nf};
  %       if contains(field,'test')
  %         continue
  %       else
  %         ctd.cast_num = [ctd.cast_num; repmat(str2double(ctd.cast{nf}),size(ctd.data.(field),1),1)];
  %       end
  %       for nv = 1:numel(vars)
  %         ctd.(vars{nv}) = [ctd.(vars{nv}); ctd.data.(field)(:,nv)];
  %       end
  %     end
  %     % convert to datenum
  %     ctd.datenum = datenum(opt.year,1,0) + ctd.timeJ;
  %     % calculate depth
  %     ctd.depth = -gsw_z_from_p(ctd.prDM,ctd.latitude);
  %     ctd.cast  = ctd.cast_num;
  %     save(opt.path.ctd_file,'ctd');
  case 'LISST_sn4025_2018_NGA_SKQ201810S_spring'

    %% CTD files
    % ctd_file = fullfile(opt.path.ctddir,'SKQ201810S.ascii');
    opt.path.ctd_file = fullfile(opt.path.ctddir,'NGA_SKQ201810S_ctd_L3_v2.mat');
    opt.path.ctd_file_unbinned = fullfile(opt.path.ctddir,'NGA_SKQ201810S_ctd_L0_v1','ctd.mat');
    %     ctd_file = fullfile(opt.path.ctddir,'NGA_SKQ201810S_ctd_L3_v2.csv');
    %     fileopts = detectImportOptions(ctd_file,'FileType','text','CommentStyle','%');
    %     % read variable names
    %     % fileopts.VariableNames(1:3) = {'cast' 'pressure' 't090C'};
    %     fileopts.VariableNames(8) = {'cast'};
    %     fileopts.VariableNames(10) = {'pressure'};
    %     fileopts.VariableNames(11) = {'t090C'};
    %     fileopts.VariableNames(5:6) = {'lon' 'lat'};
    %     % fileopts.VariableNames{26}  = 'timeJ';
    %     fileopts.VariableNames{40}  = 'timeJ';
    %
    %     ctd = readtable(ctd_file,fileopts);
    %     % throw out all comment lines because all NaN
    %     firstline = find(isfinite(ctd.cast),1);
    %     ctd = ctd(firstline:end,:);
    %     ctd = table2struct(ctd,'ToScalar',true);
    %     ctd.datenum = datenum(opt.year,1,0) + ctd.timeJ;
    %     ctd.depth = -gsw_z_from_p(ctd.pressure,ctd.lat);
    %     ctd.temp  = ctd.t090C;
    %     save(opt.path.ctd_file,'ctd','-v7.3');
    %
    
    %     % older versions?
    %     %   case 'LISST_sn4025_2019_NGA_SKQ201915S_summer'
    %     %     % ctd_file = fullfile(opt.path.ctddir,'SKQ201810S.ascii');
    %     %     ctd_file = fullfile(opt.path.ctddir,'NGA_SKQ201915S_ctd_L3_v1.csv');
    %     %     fileopts = detectImportOptions(ctd_file)%,'FileType','text','CommentStyle','%');
    %     %     % read variable names
    %     %     % fileopts.VariableNames(1:3) = {'cast' 'pressure' 't090C'};
    %     %     fileopts.VariableNames(8) = {'cast'};
    %     %     fileopts.VariableNames(10) = {'pressure'};
    %     %     fileopts.VariableNames(11) = {'t090C'};
    %     %     fileopts.VariableNames(5:6) = {'lon' 'lat'};
    %     %     % fileopts.VariableNames{26}  = 'timeJ';
    %     %     % fileopts.VariableNames{40}  = 'timeJ';
    %     %
    %     %     ctd = readtable(ctd_file,fileopts);
    %     %     % throw out all comment lines because all NaN
    %     %     firstline = find(isfinite(ctd.cast),1);
    %     %     ctd = ctd(firstline:end,:);
    %     %     ctd = table2struct(ctd,'ToScalar',true);
    %     %     % ctd.datenum = datenum(opt.year,1,0) + ctd.timeJ;
    %     %     ctd.datenum = fileopts.VariableNames(4)
    %     %
    %     %     ctd.depth = -gsw_z_from_p(ctd.pressure,ctd.lat);
    %     %
    %     %     ctd.temp  = ctd.t090C;
    case 'LISST_sn4025_2019_NGA_SKQ201915S_summer'
      %% Background scatterfiles
      % Write in which water types you expect. Do this based on the filenames
      % of the [background_scatterfile_names].asc files. e.g. for this
      % cruise, files are FSW_13Aug18_0207.asc, MilliQ_01Sept2018_1007.asc,
      % etc., so here I chose MilliQ and FSW. This allows you to choose which
      % background scatter files to use for processing.
      % ** NOTE ** You'll have to go into
      % LISST_select_background_scatterfiles.m and add exact string to match
      % for each type in the section "Set up water types expected in files"
      fprintf('STOPPED HERE** Need to define background scatter file types... read comments above\n')
      keyboard
      opt.zscat_watertypes = {'MilliQ' 'FSW'};
      %% CTD files
      fprintf('STOPPED HERE** ALSO need to define matlab datafile where raw/unbinned ctd data is save AND binned ctd data is saved...\n')
      keyboard
      % ctd_file = fullfile(opt.path.ctddir,'SKQ201810S.ascii');
      opt.path.ctd_file = fullfile(opt.path.ctddir,'SKQ201915Sctd.mat');
      %     ctd_file = fullfile(opt.path.ctddir,'SKQ201915Sctd_csv.csv');
      %     fileopts = detectImportOptions(ctd_file,'FileType','text','CommentStyle','%');
      %     % read variable names
      %     fileopts.VariableNames(1:3) = {'cast' 'pressure' 't090C'};
      %     fileopts.VariableNames(24:25) = {'lat' 'lon'};
      %     fileopts.VariableNames{26}  = 'timeJ';
      %
      %     ctd = readtable(ctd_file,fileopts);
      %     % throw out all comment lines because all NaN
      %     firstline = find(isfinite(ctd.cast),1);
      %     ctd = ctd(firstline:end,:);
      %     ctd = table2struct(ctd,'ToScalar',true);
      %     ctd.datenum = datenum(opt.year,1,0) + ctd.timeJ;
      %     ctd.depth = -gsw_z_from_p(ctd.pressure,ctd.lat);
      %     ctd.temp  = ctd.t090C;
      %     save(opt.path.ctd_file,'ctd','-v7.3');
      
      % Data used for LISST_find_time_depth_lat
      %     ctd_file = fullfile(opt.path.ctddir,'SKQ201915Sctd_csv.csv');
      %     fileopts = detectImportOptions(ctd_file,'FileType','text','CommentStyle','%');
      %     % read variable names
      %     fileopts.VariableNames(1:3) = {'cast' 'pressure' 't090C'};
      %     fileopts.VariableNames(24:25) = {'lat' 'lon'};
      %     fileopts.VariableNames{26}  = 'timeJ';
      %     ctd = readtable(ctd_file,fileopts);
      %     % throw out all comment lines because all NaN
      %     firstline = find(isfinite(ctd.cast),1);
      %     ctd = ctd(firstline:end,:);
      %     ctd = table2struct(ctd,'ToScalar',true);
      %     ctd.datenum = datenum(opt.year,1,0) + ctd.timeJ;
      %     ctd.depth = -gsw_z_from_p(ctd.pressure,ctd.lat);
      %     ctd.temp  = ctd.t090C;
      %     save(fullfile(opt.path.ctddir,'NGA_SKQ201915S_ctd_L0_v1','ctd.mat'),'ctd');
      %     flag_cnv = 1;
      
      
      %     % older version?
      %       %     fileopts.VariableNames{24}  = 'timeJ';
      %       %         vars = {'timeJ' 'temp' 't090c' 'sal00' 'sal11' 'par' 'latitude' 'longitude' 'density00' 'depSM' 'prDM' 'flag'};
      %       %         % initialize variables in ctd structure
      %       %         for nv = 1:numel(vars)
      %       %           ctd.(vars{nv}) = [];
      %       %         end
      %       %         ctd.cast_num = [];
      %       %         cast_names = fieldnames(ctd.data);
      %       %         % combine all data into single variable
      %       %         for nf = 1:numel(cast_names)
      %       %           field = cast_names{nf};
      %       %           if contains(field,'test')
      %       %             continue
      %       %           else
      %       %             ctd.cast_num = [ctd.cast_num; repmat(str2double(ctd.cast{nf}),size(ctd.data.(field),1),1)];
      %       %           end
      %       %           for nv = 1:numel(vars)
      %       %             ctd.(vars{nv}) = [ctd.(vars{nv}); ctd.data.(field)(:,nv)];
      %       %           end
      %       %         end
      %       %         % convert to datenum
      %       %         ctd.datenum = datenum(opt.year,1,0) + ctd.timeJ;
      %       %         % calculate depth
      %       %         ctd.depth = -gsw_z_from_p(ctd.prDM,ctd.latitude);
      %       %         ctd.cast  = ctd.cast_num;
      %       %         save(ctd_dir_downcasts,'ctd');
      %
      %
      %       %     ctd_file = fullfile(opt.path.ctddir,'SKQ201915S.ascii');
      %       %     fileopts = detectImportOptions(ctd_file,'FileType','text','CommentStyle','%');
      %       %     % read variable names
      %       %     fileopts.VariableNames(1:3) = {'cast' 'pressure' 't090C'};
      %       %     fileopts.VariableNames(24:25) = {'lat' 'lon'};
      %       %    % fileopts.VariableNames{26}  = 'timeJ';
      %       %
      %       %     ctd = readtable(ctd_file,fileopts);
      %       %     % throw out all comment lines because all NaN
      %       %     firstline = find(isfinite(ctd.cast),1);
      %       %     ctd = ctd(firstline:end,:);
      %       %     ctd = table2struct(ctd,'ToScalar',true);
      %       %    % ctd.datenum = datenum(opt.year,1,0) + ctd.timeJ;
      %       %     ctd.depth = -gsw_z_from_p(ctd.pressure,ctd.lat);
      %       %     ctd.temp  = ctd.t090C;
      %       %     save(opt.path.ctd_file,'ctd','-v7.3');
      
end %% function ctd = read_CTD_NGA_LTER_L0(ctd_directory)
