% ntools_procDecimateIEEG(recording_filename_root)
%   processes decimates raw nspike file to give output files
%
%  Inputs: recording_filename_root - recording root without the suffixes.
%
% E.g.:
%   recording_filename_root = '/mnt/esata/NY140/SP/081030/rec005';
%   ntools_decimate(recording_filename_root);
%
function ntools_procDecimateIEEG(recording_filename_root)

    global experiment

    recording_filename_root = make_filename_root(recording_filename_root); %removes suffix if necess

    SIZEOF_FLOAT = 4;
    SIZEOF_INT = 2;
    CH = length(experiment.channels);
    decimate_factor = round(experiment.recording.sample_rate./experiment.processing.ieeg.sample_rate);
    IEEGFS = experiment.recording.sample_rate./decimate_factor;
    
    T = round(1.*1024.^2/(IEEGFS.*CH*SIZEOF_INT));
    nspike_file = [recording_filename_root '.nspike.dat'];
    if length(dir(nspike_file))
        nspike_fid = fopen(nspike_file, 'r');
        ieeg_fid = fopen([recording_filename_root '.decieeg.dat'], 'w');
        chk=1; N=0;
        ieeg = zeros(1, CH*IEEGFS*T, 'int16');
        while(chk)
            N = N+1;
            disp(['IEEG: Loop ' num2str(N)]);
            ieeg = fread(nspike_fid, CH*IEEGFS*T, [num2str(CH) '*int16=>int16'], (decimate_factor-1)*CH*SIZEOF_INT);
            fwrite(ieeg_fid, ieeg, 'short');
            if (size(ieeg, 1)/CH < IEEGFS*T); chk = 0; end
        end
        fclose(ieeg_fid);
        fclose(nspike_fid);
    else
        disp('No nspike file found.  Skipping.')
    end
end
