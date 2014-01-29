function [nspike_fid,comedi_fid]=cap5m
ffp_start_record('/mnt/raid/adam/comedidata', '/mnt/raid/adam/nspikedata')
pause(300);
ffp_stop_record
nspike_fid = fopen('/mnt/raid/adam/nspikedata');
%nspike_data = fread(fid, [256,inf], 'int16');
%fclose(fid);
comedi_fid = fopen('/mnt/raid/adam/comedidata');
%comedi_data = fread(fid,[7,inf],'ushort');
%fclose(fid);
%figure;
%imagesc(nspike_data);
%figure
%imagesc(comedi_data);

