function [nspike_data,comedi_data]=cap5
n_start_time = gettime_nspike;
c_start_time = gettime_comedi;
ffp_start_record('/mnt/raid/adam', 'cap5')
pause(5);
ffp_stop_record
disp('nspike time elapsed: ')
disp(gettime_nspike - n_start_time)
disp('comedi time elapsed: ')
disp(gettime_comedi - c_start_time)
fid = fopen('/mnt/raid/adam/cap5.nspike.dat');
nspike_data = fread(fid, [128,inf], 'int16');
fclose(fid);
fid = fopen('/mnt/raid/adam/cap5.comedi.dat');
comedi_data = fread(fid,[8,inf],'ushort');
fclose(fid);

comedi_data2 = comedi_data(1,:) - 32768;
comedi_data2 = comedi_data2 * 2;

figure
plot(nspike_data(255,:));
hold
plot(comedi_data2(1,:),'r');




%figure;
%imagesc(nspike_data);
%figure
%imagesc(comedi_data);

