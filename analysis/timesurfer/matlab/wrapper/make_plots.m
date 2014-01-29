function make_plots(varargin)
% make_plots('event',[1:8],'format','jpg','title','NY70-FWIO')
parms = mmil_args2parms(varargin,...
						{'filepath',pwd,[],...
						 'datafile',[],[],...
						 'tsobject','timefreq_data',[],...
						 'datafield','timefreq',[],...
						 'dataparam','power',[],...
						 'event',[],[],...
						 'event_str','',[],...
                         'conditionfile',[],[],...
                         'cond_key',[],[],...
                         'cond_key_file',[],[],...
						 'nr',8,[],...
						 'nc',8,[],...
						 'cpad',0.01,[],...
						 'rpad',0.01,[],...
						 'cond',1,[],...
						 'zscore',0,[],...
						 'layout',[],[],...
						 'chantype','eeg',[],...
                         'badchans',[],[],...
                         'badchanfile',[],[],...
                         'showbadchans','yes',{'yes','no'},...
						 'multiplot',1,{1,0},...
						 'topoplot',0,{1,0},...
						 'singleplot',0,{1,0},...
						 'channel',[],[],...
						 'showlabels','yes',{'yes','no'},...
						 'zlim',[],[],...
						 'xlim',[-.5 1],[],...
						 'ylim',[],[],...
						 'interactive','no',[],...
						 'stim','yes',{'yes','no'},...
						 'colorbar','yes',[],...
						 'title','',[],...
						 'format','eps',[],...
						 'fontsize',8,[],...
						 'skip_multievent_files',1,{1,0},...
						 'savepath',[],[],...
						 'save',1,{1,0},...
						 'close',1,{1,0},...
						},false);

cwd = pwd;
cd(parms.filepath);
layout_flag = 0;
if strcmp(parms.format,'jpg'), parms.printcmd = {'-djpeg'}; end
if strcmp(parms.format,'eps'), parms.printcmd = {'-depsc','-tiff','-r150'}; end
if strcmp(parms.format,'tif'), parms.printcmd = {'-dtiff'}; end
if strcmp(parms.chantype,'eeg') && isempty(parms.layout), layout_flag = 1; 
elseif strcmp(parms.chantype,'eeg'), layfile{1} = parms.layout;
else layfile{1} = 'meg';
end

% find matfiles
if isempty(parms.datafile)
	parms.datafile = what(parms.filepath);
	parms.datafile = parms.datafile.mat;
end
if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end
if ~iscell(parms.badchans), parms.badchans = {parms.badchans}; end

% load condition info
if ~isempty(parms.cond_key)
    idx = find(ismember([parms.cond_key.event_code],parms.event));
    parms.event_str = [parms.cond_key.name];
    parms.event_str = parms.event_str(idx);
elseif ~isempty(parms.cond_key_file) && exist(parms.cond_key_file,'file')
    load(parms.cond_key_file);
    idx = find(ismember([cond_key.event_code],parms.event));
    parms.event_str = [cond_key.name];
    parms.event_str = parms.event_str(idx);
elseif ~isempty(parms.conditionfile) && exist(parms.conditionfile,'file')
    ts_makecondkey(parms.conditionfile);
    [pathstr,name] = fileparts(parms.conditionfile);
    load(fullfile(pathstr,'cond_key.mat'));
    idx = find(ismember([cond_key.event_code],parms.event));
    parms.event_str = [cond_key.name];
    parms.event_str = parms.event_str(idx);
end 

% loop over files
for f = 1:length(parms.datafile)
	fprintf('scanning file %g of %g: %s\n',f,length(parms.datafile),parms.datafile{f});
	% load timesurfer data
	data = getfield(load(parms.datafile{f},parms.tsobject),parms.tsobject);
	if parms.skip_multievent_files && length(data.(parms.datafield)) > 1, continue; end
    % set parms.badchans to bad
    if ~isempty(parms.badchans{1})
      [bidx dummy] = match_str({data.sensor_info.label},parms.badchans);
      [data.sensor_info(bidx).badchan] = deal(1);
    end
    % read badchanfile
    if ~isempty(parms.badchanfile) && exist(parms.badchanfile,'file')
      bidx = ts_read_txt_badchans(parms.badchanfile,{data.sensor_info.label});
      [data.sensor_info(bidx).badchan] = deal(1);
    end
    % showbadchans?
    if strcmp(parms.showbadchans,'yes')
        bidx = find([data.sensor_info.badchan]);
        [data.sensor_info(bidx).badchan] = deal(0);
    end
	% get event indices
	if isempty(parms.event)
		cond = 1:length(data.(parms.datafield));
	else
		cond = find(ismember([data.(parms.datafield).event_code],parms.event));
	end
	
	% fileparts for output filename
	[pathstr,name,ext,vrsn] = fileparts(parms.datafile{f});
  if isempty(pathstr), pathstr = cwd; end
	if isempty(parms.savepath), parms.savepath = sprintf('%s/images',pathstr); end
	if ~exist(parms.savepath,'file'), mkdir(parms.savepath); end
	
	% loop over events
	for c = 1:length(cond)
		eventcode = data.(parms.datafield)(cond(c)).event_code;
		fprintf('plotting event %g\n',eventcode);
		% add event string to title
		if ~isempty(parms.event_str) && length(parms.event_str) == length(parms.event)
			title_str = sprintf('%s\n%s',parms.title,parms.event_str{parms.event==eventcode});
		else
			title_str = parms.title;
		end		
		% zscore
		if parms.zscore
			data.(parms.datafield)(cond(c)).(parms.dataparam) = ts_zscore(data,'cond',cond(c));
			zstr = '-zscore'; else zstr = '';
		end
		% write layout file
		if layout_flag
			layfile = write_layout(parms.savepath,data.sensor_info,parms.nr,parms.nc,parms.rpad,parms.cpad);
		end
		for p = 1:length(layfile)
			% convert data to fieldtrip
			ftdata = ts_data2fieldtrip(data,'chantype',parms.chantype,'condition',cond(c));
			% fieldtrip configuration
			cfg=[]; 
			if ~strcmp(layfile{p},'meg'), cfg.layout = layfile{p}; end
			cfg.showlabels = parms.showlabels;
      if ischar(parms.xlim) && strcmp(parms.xlim,'all')
        parms.xlim = [ftdata.time(1) ftdata.time(end)];
      end
			if ~isempty(parms.zlim), cfg.zlim = parms.zlim; end
			if ~isempty(parms.xlim), cfg.xlim = parms.xlim; end
			if ~isempty(parms.ylim), cfg.ylim = parms.ylim; end
			cfg.interactive = parms.interactive;
			cfg.colorbar = parms.colorbar;
			cfg.stim = parms.stim;
			cfg.fontsize = parms.fontsize;
			cfg.comment = sprintf('%s\n%s%s\n%s\n',title_str,parms.dataparam,zstr,date);		
			%cfg.title = title_str;	
			% topoplot
			if parms.topoplot && strcmp(parms.chantype,'eeg')
				if p~=2, continue; end  % skip if not grid
				outfile = sprintf('%s/%s.topoplot%s_sec%g-%g_page%g.%s',parms.savepath,name,zstr,parms.xlim(1),parms.xlim(2),p,parms.format);
				if exist(outfile,'file')
					fprintf('%s: skipping figure.  File already exists: %s\n',mfilename,outfile);
					continue;
				end
				if strcmp(parms.tsobject,'timefreq_data')
					figure; topoplotER_iEEG(cfg,ftdata);
				else
					figure; topoplotER(cfg,ftdata);
				end
				if parms.save
					print(gcf,parms.printcmd{:},outfile);				
				end
			end			
			% multiplot
			if parms.multiplot
				outfile = sprintf('%s/%s.multiplot%s_sec%g-%g_page%g.%s',parms.savepath,name,zstr,parms.xlim(1),parms.xlim(2),p,parms.format);
				if exist(outfile,'file')
					fprintf('%s: skipping figure.  File already exists: %s\n',mfilename,outfile);
					continue;
				end
				if strcmp(parms.tsobject,'timefreq_data')
					figure; multiplotTFR(cfg,ftdata);
				else
					figure; multiplotER(cfg,ftdata);
				end
				if parms.save
					print(gcf,parms.printcmd{:},outfile);				
				end
            end			
			% singleplot
			if parms.singleplot
				outfile = sprintf('%s/%s.singleplot.%s_page%g.%s',parms.savepath,name,parms.channel{1},p,parms.format);			
				if exist(outfile,'file')
					fprintf('%s: skipping figure.  File already exists: %s\n',mfilename,outfile);
					continue;
				end
				if ~isempty(parms.channel), cfg.channel = parms.channel; end
				if strcmp(parms.tsobject,'timefreq_data')
					figure; singleplotTFR(cfg,ftdata);
				else
					figure; singleplotER(cfg,ftdata);
				end
				if parms.save
					print(gcf,parms.printcmd{:},outfile);
				end
			end 
		end % end figure loop
		if parms.close, close all; end
	end % end event loop
end % end file loop

cd(cwd);
clear all
