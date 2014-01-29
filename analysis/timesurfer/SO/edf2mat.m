% on ip84
addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
addpath(genpath('/home/jsherfey/svn/dev/packages/fieldtrip-20100511'))

parms.datafile = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG19/MG19_sleep1.edf';
parms.datafile = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG25/MG2501Sleep1.edf';
parms.datafile = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG20/MG20_Sleep1.edf';

parms.datafile = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG21/MG21_Sleep1.edf';
epoch_data = ts_loadedf('datafile',parms.datafile);
[fpath,fname] = fileparts(parms.datafile);
outpath = strrep(fpath,'/space/md3/5/halgdev/projects/Spindles/Intracranial_Data','/space/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg');
outname = [fname '.mat'];
if ~exist(outpath),mkdir(outpath); end
save(fullfile(outpath,outname),'data','-v7.3');

parms.datafile = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG22/MG21_Sleep1.edf';
epoch_data = ts_loadedf('datafile',parms.datafile);
[fpath,fname] = fileparts(parms.datafile);
outpath = strrep(fpath,'/space/md3/5/halgdev/projects/Spindles/Intracranial_Data','/space/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg');
outname = [fname '.mat'];
if ~exist(outpath),mkdir(outpath); end
save(fullfile(outpath,outname),'data','-v7.3');

parms.datafile = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG23/MG21_Sleep1.edf';
epoch_data = ts_loadedf('datafile',parms.datafile);
[fpath,fname] = fileparts(parms.datafile);
outpath = strrep(fpath,'/space/md3/5/halgdev/projects/Spindles/Intracranial_Data','/space/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg');
outname = [fname '.mat'];
if ~exist(outpath),mkdir(outpath); end
save(fullfile(outpath,outname),'data','-v7.3');

parms.datafile = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG24/MG21_Sleep1.edf';
epoch_data = ts_loadedf('datafile',parms.datafile);
[fpath,fname] = fileparts(parms.datafile);
outpath = strrep(fpath,'/space/md3/5/halgdev/projects/Spindles/Intracranial_Data','/space/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg');
outname = [fname '.mat'];
if ~exist(outpath),mkdir(outpath); end
save(fullfile(outpath,outname),'data','-v7.3');






% /home/halgdev/incoming/MEG_MGH/MG36_MEG/Sleep01_raw.fif % 1-8, 10-12


