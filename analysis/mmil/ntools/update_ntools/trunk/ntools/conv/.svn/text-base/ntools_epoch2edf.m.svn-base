% ntools_epoch2edf(data, filename)
%   convert data or epoch_data to edf_data and save it.
%
% Parameters
%   data - epoch_data/cont_data structure
%   filename - the location where you want the edf_data to be stored.
%
% example: ntools_epoch2edf(data, filename);
%

function ntools_epoch2edf(data, filename)
    %write the header
    fid = fopen(filename, 'wt');
    num_sensors = data.num_sensors;

    %8 ascii : version of this data format (0)
    fprintf(fid, '%-8s', '0');

    %80 ascii : local patient identification
    %80 ascii : local recording identification
    %8 ascii : startdate of recording (dd.mm.yy)
    %8 ascii : starttime of recording (hh.mm.ss)
    fprintf(fid, '%-176s', '                                                                                Nspike                                                                                          ');

    %8 ascii : number of bytes in header record
    fprintf(fid, '%-8d', 256+256*num_sensors);

    %44 ascii : reserved
    fprintf(fid, '%-44s', 'Exported from Neuroscan SCAN 4 software');

    %8 ascii : number of data records (-1 if unknown)
    fprintf(fid, '%-8s', '-1');

    %8 ascii : duration of a data record, in seconds
    [m n] = size(data.epochs.time);
    duration = data.epochs.time(1, n) - data.epochs.time(1, 1);
    fprintf(fid, '%-8f', duration);

    %4 ascii : number of signals (num_sensors) in data record
    fprintf(fid, '%-4d', num_sensors);

    %num_sensors * 16 ascii : num_sensors * label (e.g. EEG FpzCz or Body temp)num_sensors=128

    for i = 1:num_sensors
        a = data.sensor_info(1, i).label;
        fprintf(fid, '%-16s', a);
    end



    %num_sensors * 80 ascii : num_sensors * transducer type (e.g. AgAgCl electrode)
    %%%can I use repmat?

    fprintf(fid, '%-s', repmat('Unknown. Usually EEG electrode.                                                 ', 1, num_sensors));

    %num_sensors * 8 ascii : num_sensors * physical dimension (e.g. uV or degreeC)
    %char(63, 82, 32, 32, 32, 32, 32);
    fprintf(fid, '%-s', repmat(char(181, 86, 32, 32, 32, 32, 32, 32), 1, num_sensors));

    %num_sensors * 8 ascii : num_sensors * physical minimum (e.g. -500 or 34)
    %fprintf(fid, '%-s', repmat('-32767  ', 1, num_sensors));
    phys_min = round(min(min(data.epochs.data)));
    phys_max = round(max(max(data.epochs.data)));
    for i = 1:num_sensors
        fprintf(fid, '%-8d', phys_min);
    end

    %num_sensors * 8 ascii : num_sensors * physical maximum (e.g. 500 or 40)
    %fprintf(fid, '%-s', repmat('32768   ', 1, num_sensors));
    for i = 1:num_sensors
        fprintf(fid, '%-8d', phys_max);
    end
    %num_sensors * 8 ascii : num_sensors * digital minimum (e.g. -2048)
    fprintf(fid, '%-s', repmat('-32767  ', 1, num_sensors));

    %num_sensors * 8 ascii : num_sensors * digital maximum (e.g. 2047)
    fprintf(fid, '%-s', repmat('32768   ', 1, num_sensors));

    %num_sensors * 80 ascii : num_sensors * prefiltering (e.g. HP:0.1Hz LP:75Hz)
    fprintf(fid, '%-s', repmat('HP:     LP:                                                                     ', 1, num_sensors));

    %num_sensors * 8 ascii : num_sensors * nr of samples in each data record
    nr = n;
    for i=1:num_sensors
        fprintf(fid, '%-8d', nr);
    end


    %num_sensors * 32 ascii : num_sensors * reserved
    fprintf(fid, '%-s', repmat(' ', 1, 32*num_sensors));

    %write the data.
    fclose(fid);
    
    fid = fopen(filename, 'a+b');
    D = double(data.epochs.data)';


    PhysMin = min(D);
    PhysMax = max(D);
    DigMin = repmat(-32767, 1, num_sensors);
    DigMax = repmat(32768, 1, num_sensors);

    b = min(PhysMin);
    a = max(PhysMax);
    if b ~= 32768 || a ~= 32768;
        D = D - repmat(PhysMin(:)', size(D, 1), 1);
        D = D * sparse(1:num_sensors, 1:num_sensors, (DigMax - DigMin) ./ (PhysMax - PhysMin)); % scale Phys->Dig
        D = D + repmat(DigMin(:)', size(D, 1), 1);
    end

    fwrite(fid, D, 'int16');
    fclose(fid);
end

