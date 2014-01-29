function [comedi_data, nspike_data] = compare_ts(secs) 

disp(secs)
nspike_fid = fopen('/mnt/raid/adam/cap5.nspike.dat');
comedi_fid = fopen('/mnt/raid/adam/cap5.comedi.dat');

%fseek(nspike_fid, secs * 3e4 * 128, 'bof');
%fseek(comedi_fid, secs * 3e4 * 8, 'bof');

nspike_data = fread(nspike_fid, [128,5*3e4], 'int16');
comedi_data = fread(comedi_fid,[8,5*3e4],'ushort');

fclose(nspike_fid);
fclose(comedi_fid);

comedi_data2 = comedi_data(1,:) - 32768;
comedi_data2 = comedi_data2 * 2;

figure
plot(nspike_data(1,:));
hold
plot(comedi_data2(1,:),'r');


