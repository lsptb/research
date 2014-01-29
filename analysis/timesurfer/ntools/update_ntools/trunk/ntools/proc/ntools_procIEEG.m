% ntools_procIEEG(recording_filename_root)
%   processes raw nspike file to give iEEG files
%
%  Inputs: recording_filename_root - recording root without the suffixes.
%
% E.g.:
%   recording_filename_root = '/mnt/esata/NY140/SP/081030/rec005';
%   procIEEG(recording_filename_root);
%
function ntools_procIEEG(recording_filename_root)

    global experiment

    recording_filename_root = make_filename_root(recording_filename_root); %removes suffix if necess
    
    figh = get(0, 'CurrentFigure'); %use setappdata to report percent done if in GUI

    FS = experiment.recording.sample_rate;
    CH = length(experiment.channels);
    decimate_factor = round(experiment.recording.sample_rate./experiment.processing.ieeg.sample_rate);

    tic
    T = 5;
    tapers = [0.0025, experiment.processing.ieeg.lowpass];
    n = tapers(1); w = tapers(2); p = n*w; k = floor(2*p-1); tapers = [n, p, k];
    tapers(1) = tapers(1).*FS;
    tapers = dpsschk(tapers);
    filt = mtfilt(tapers, FS, 0);
    filt = single(filt./sum(filt));
    Nf = length(filt);
    nspike_file = [recording_filename_root '.nspike.dat'];

    if isfile(nspike_file)
        nspike_fid = fopen(nspike_file, 'r');
        nspike_size = file_size(nspike_file);
        ieeg_fid = fopen([recording_filename_root '.ieeg.dat'], 'w');

        chk = 1;
        N = 0;
        ieeg = zeros(CH, FS*T+Nf-1, 'single');
        while(chk)
            tic
            N = N+1;
            disp(['IEEG: Loop ' num2str(N)]);
            data = fread(nspike_fid, [CH, FS*T], 'int16=>single');

            if(size(data, 2))
                for j = 1:size(data, 1)
                    ieeg(j, 1:size(data, 2)+Nf-1) = conv(data(j,:), filt);
                end
                ieeg2 = ieeg(:, Nf./2:decimate_factor:size(data, 2)+Nf./2-1);
                fwrite(ieeg_fid, ieeg2, 'float');
            end
            
            if(~isempty(figh)), setappdata(figh, 'percent_done', ftell(nspike_fid)/nspike_size); end
            
            if (size(data, 2) < FS*T)
                chk = 0; 
            end
            toc
        end
        fclose(ieeg_fid);
        fclose(nspike_fid);
    else
        disp('No nspike file found.  Skipping.')
    end
end