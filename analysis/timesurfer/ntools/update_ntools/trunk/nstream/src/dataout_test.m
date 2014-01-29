function [nspike_data,comedi_data]=dataout_test
        fid = fopen('nspikedata.out');
        nspike_data = fread(fid, [256,inf], 'int16');
        fclose(fid);
        fid = fopen('comedidata.out');
        comedi_data = fread(fid,[7,inf],'ushort');
        fclose(fid)
%        for i = 1:256
%            data2(i,:) = data(i,:)+i*100;
%        end
%        plot(data2')

