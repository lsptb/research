function [nspike_data,comedi_data]=cap5
ffp_start_record('comedidata', 'nspikedata')
pause(15);
ffp_stop_record
fid = fopen('../nspikedata');
nspike_data = fread(fid, [256,inf], 'int16');
fclose(fid);
fid = fopen('../comedidata');
comedi_data = fread(fid,[7,inf],'ushort');
fclose(fid);
figure;
imagesc(nspike_data);
figure
imagesc(comedi_data);

