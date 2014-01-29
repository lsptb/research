%% make averages

event = [1:8];
nchan = 126;
dataname = 'NY50-00_stream2.ica.rej.foi8-200.toi-0.5-1.timefreq.event';
inpath   = '/space/md3/5/halgdev/analysis/iEEG_NYU/NY50/NY50_FWIO/matfiles/timefreq';
outpath  = '/space/mdeh1/7/halgdev/analysis/iEEG_NYU/NY50/NY50_FWIO/matfiles/timefreq';

for i = 1:length(event)
	datapath=sprintf('%s/event%g',inpath,event(i));
	for j = 1:nchan
		inname=sprintf('%s%g.channel_%03i.mat',dataname,i,j);
		fprintf('loading event %g of %g, channel %g of %g\n',i,length(event),j,nchan);
		tfdata = getfield(load(fullfile(datapath,inname),'timefreq_data'),'timefreq_data');
		if j == 1
			timefreq_data = tfdata;
			timefreq_data.num_sensors = nchan;
			timefreq_data.timefreq.power = mean(tfdata.timefreq.power,4);
		else
			timefreq_data.timefreq.power = cat(1,timefreq_data.timefreq.power,mean(tfdata.timefreq.power,4));
			timefreq_data.sensor_info = cat(1,timefreq_data.sensor_info,tfdata.sensor_info);
		end
		clear tfdata;
	end
		outname = sprintf('%s%g.mat',dataname,i);
		fprintf('saving average power for event %g: %s\n',i,fullfile(outpath,outname));
		save(fullfile(outpath,outname),'timefreq_data');
		clear timefreq_data;
end	

