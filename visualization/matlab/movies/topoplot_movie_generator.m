% Slow-wave sleep - interpolated gradiometer movies

files = {...
  '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/matfiles/sleep_s2_2_raw_grad1.mat',...
  '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/matfiles/sleep_s3_4_raw_grad1.mat',...
  '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/matfiles/sleep_s4_18_raw_grad1.mat',...
  '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/s5/matfiles/SL_2_nd01_060726_grad1.mat',...
  '/home/halgdev/projects/jsherfey/sleep/s8/matfiles/SL_2_nb01_060808_grad1.mat',...
  };
subjs  = {'s2_2','s3_4','s4_18','s5_2','s8_2'};  
subjid = [2 3 4 5 8];
% toilim = {[950 1050],[950 1050],[400 600],[600 680],[200 400]};
% toilim = {[900 1200],[900 1200],[350 650],[500 800],[150 450]}; % 5 min
toilim = {[950 1000],[950 1000],[450 500],[600 650],[300 350]}; % 50 sec

%% BEGIN TOPOPLOT MOVIE GENERATOR
% AVI: open figure, {add data, 1 frame at a time} CLOSE
tic
tocs    = [];
tframe  = .01; % time for one frame (sec)
fps     = 1/tframe; % frames per second (dt=100ms ==> fps=1/dt=10)
nframes = 1500;
dsfact  = 1;
s = 0;
for subjnum = subjid
  s = s + 1;
  prefix  = sprintf('sleep_%s_grad1',subjs{s});
  str     = sprintf('%s: grad1 movie',subjs{s});
  outpath = '/mdkm1/2/kmdev/projects/jsherfey/sleep/movies';
  lay     = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
  badchan = sprintf('/home/halgdev/projects/jsherfey/sleep/badchans/s%g_badchans.txt',subjnum);
  avi     = sprintf('%s/%s_SlowOscillations_%gframetarget.avi',outpath,subjs{s},nframes);
  % load data (start w/ 1 file for 1st test; don't need to use SO_combine_data)
  load(files{s});
  data = ts_data_selection(data,'badchanfile',badchan,'toilim',toilim{s});
  data.averages = data.epochs;
  data = rmfield(data,'epochs');
  % create avi object
  aviobj = avifile(avi,'fps',fps);

  % determine zlim based on maxmin of all of data
  zmax = mean([data.averages.data(:)])+std([data.averages.data(:)])*2; % max([data.averages.data(:)])
  zmin = mean([data.averages.data(:)])-std([data.averages.data(:)])*2; % min([data.averages.data(:)])
  zlim = [zmin zmax];
  % zlim    = 'maxmin';
  % determine number of intervals and selection sequence
  t0 = data.averages.time(1);
  tf = data.averages.time(end);
  dt = 1/data.sfreq;
  t  = t0:dt:tf;
  T  = tf-t0;
  n  = floor(T/tframe);
  % % 22 minutes (1 file): 827400 samples
  % % time for 22 minutes with dt=100ms (1370 frames), if 30s/frame: 11.4 hours

  % filtered signals
  filt1 = ts_preproc(data,'bpfilter','yes','bpfreq',[.3 3]);
  if dsfact > 1, filt1 = ts_preproc(filt1,'dsfact',dsfact); end
  zlim1 = [mean([filt1.averages.data(:)])-std([filt1.averages.data(:)])*2 mean([filt1.averages.data(:)])+std([filt1.averages.data(:)])*2];
  filt2 = ts_preproc(data,'bpfilter','yes','bpfreq',[3 8]);
  if dsfact > 1, filt2 = ts_preproc(filt2,'dsfact',dsfact); end
  zlim2 = [mean([filt2.averages.data(:)])-std([filt2.averages.data(:)])*2 mean([filt2.averages.data(:)])+std([filt2.averages.data(:)])*2];
  % filt3 = ts_preproc(data,'lnfilter','yes','lnfreq',[120 180]);
  filt3 = ts_preproc(data,'bpfilter','yes','bpfreq',[30 50]);
  if dsfact > 1, filt3 = ts_preproc(filt3,'dsfact',dsfact); end
  zlim3 = [mean([filt3.averages.data(:)])-std([filt3.averages.data(:)])*2 mean([filt3.averages.data(:)])+std([filt3.averages.data(:)])*2];
  if dsfact > 1, data  = ts_preproc(data,'dsfact',dsfact);  end

  % loop over small time intervals of t = [0 eof]: frames
  saveflag = 0;
  tf       = t0 + tframe;
  sz       = get(0,'ScreenSize');
  figure('Color','w','Position',[.1*sz(3) .1*sz(4) .6*sz(3) .8*sz(4)]);
  tocs  = [tocs toc];
  for k = 1:nframes
    if tf > data.averages.time(end), break; end
    % select data subset
    dat = ts_data_selection(data,'toilim',[t0 tf]);
    subplot(2,2,1)
    ts_ezplot(dat,'showlabels','no','layout',lay,'events',1,'comment','xlim','zlim',zlim,'topoplot',1,'avgovertime','yes','title',[str ''],'newfig',0,'save',saveflag,'outpath',outpath,'prefix',prefix,'overwrite',1,'close',0,'colorbar','no','style','straight','shading','flat','electrodes','on');
    tmp = ts_data_selection(filt1,'toilim',[t0 tf]);
    subplot(2,2,2)
    ts_ezplot(tmp,'showlabels','no','layout',lay,'events',1,'comment','xlim','zlim',zlim1,'topoplot',1,'avgovertime','yes','title',[str ''],'newfig',0,'save',saveflag,'outpath',outpath,'prefix',prefix,'overwrite',1,'close',0,'colorbar','no','style','straight','shading','flat','electrodes','on');
    tmp = ts_data_selection(filt2,'toilim',[t0 tf]);
    subplot(2,2,3)
    ts_ezplot(tmp,'showlabels','no','layout',lay,'events',1,'comment','xlim','zlim',zlim2,'topoplot',1,'avgovertime','yes','title',[str ''],'newfig',0,'save',saveflag,'outpath',outpath,'prefix',prefix,'overwrite',1,'close',0,'colorbar','no','style','straight','shading','flat','electrodes','on');
    tmp = ts_data_selection(filt3,'toilim',[t0 tf]);
    subplot(2,2,4)
    ts_ezplot(tmp,'showlabels','no','layout',lay,'events',1,'comment','xlim','zlim',zlim3,'topoplot',1,'avgovertime','yes','title',[str ''],'newfig',0,'save',saveflag,'outpath',outpath,'prefix',prefix,'overwrite',1,'close',0,'colorbar','no','style','straight','shading','flat','electrodes','on');

    % get and add frame
    F      = getframe(gcf);
    aviobj = addframe(aviobj,F);
    % update limits
    t0 = t0 + tframe;
    tf = tf + tframe;
    clear F tmp dat
    tocs = [tocs toc];
  end
  close
  % close avi object
  aviobj = close(aviobj);
  clear data
  toc

  % mov = aviread(avi);
  % movie(mov,1,fps);

  tocs = [tocs toc];
  toc
  
  clear filt1 filt2 filt3

end % end subject loop

%% MPEG movie maker (function: mpgwrite)
% addpath(genpath('/home/jsherfey/svn/dev/movies'));
% mov = aviread(avi);
% mpgwrite(mov,jet,fullfile(fileparts(avi),'outfile.mpg'));


