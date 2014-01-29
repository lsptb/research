% control parameters
parms.skip_multievent_files = 1;
parms.tsobject = 'timefreq_data';
parms.datafield = 'timefreq';
parms.dataparam = 'power';
parms.nr = 10;
parms.rpad = 0.01;
parms.nc = 10;
parms.cpad = 0.01;
parms.zscore = 1;
parms.chantype = 'eeg';
parms.layout = [];
parms.multiplot = 0;
parms.singleplot = 0;
parms.topoplot = 1;
parms.channel = [];
parms.save = 0;
parms.format = 'eps';
parms.close = 0;

% display parameters
parms.zlim = [-5 5]; zstr = 'zscore-5-5';
%parms.zlim = [-3 3]; zstr = 'zscore-3-3';
%parms.xlim = [-.2 .6];
parms.xlim = [-.5 1];

parms.stim = 'yes';
parms.colorbar = 'yes';
parms.showlabels = 'yes';
parms.fontsize = 8;

parms.xlabel = 'time (sec)';
parms.ylabel = 'frequency (Hz)';


% batch parameters
subj = {'50'};%{'50','59','62','68','70','87','90','101','133'};
event = [3];%[1:8];
event_str = {'Output New'};% {'Input New','Input Old','Output New','Output Old','NPNW','FF','Target','Output Words'};

args = mmil_parms2args(parms);

for s = 1:length(subj)
	fprintf('plotting %s for NY%s\n',parms.tsobject,subj{s});
	filepath = sprintf('/space/mdeh1/7/halgdev/analysis/iEEG_NYU/NY%s/NY%s_FWIO/matfiles/timefreq',subj{s},subj{s});
	savepath = sprintf('/space/mdeh1/7/halgdev/analysis/iEEG_NYU/NY%s/NY%s_FWIO/images/new/%s',subj{s},subj{s},zstr);
	title_str = sprintf('NY%s-FWIO',subj{s});
	make_plots('event',event,'event_str',event_str,'title',title_str,...
						 'filepath',filepath,'savepath',savepath,args{:});
end

