% ntools_procComediVideoSync(recording_filename_root)
%   processes raw comedi file to give video sync file
%
%  Inputs: recording_filename_root - recording root without the suffixes.
%
% E.g.:
%   recording_filename_root = '/mnt/esata/NY140/SP/081030/rec005';
%   ntools_procComediVideoSync(recording_filename_root);
%
%  Saves [recording_filename_root '.videosync.txt'] file
%      First column is timestamp in seconds from beginning of comedi file
%      Second column is code.
%

function ret=ntools_procComediVideoSync(recording_filename_root)

    global experiment

    recording_filename_root = make_filename_root(recording_filename_root); %removes suffix if necess
    
    figh = get(0, 'CurrentFigure'); %use setappdata to report percent done if in GUI

    tic
    videosync_channel = experiment.recording.comedi.videosync_channel_index;
    comedi_sample_format = experiment.recording.comedi.sample_format;
    CH = experiment.recording.comedi_num_channels_to_write;
    FS = experiment.recording.comedi.sample_rate;
    T = 0.25;
    comedi_sample_size = experiment.recording.comedi.sample_size;
    THRESHOLD = 1e3;
    new_sample_rate = experiment.processing.comedi.sample_rate;
    DECIMATION_FACTOR = round(FS./new_sample_rate);
    comedi_file = [recording_filename_root '.comedi.dat'];
    prev_stdcomedi = 0;
    prev_code = 0;
    if isfile(comedi_file)
        comedi_fid = fopen(comedi_file, 'r');
        fseek(comedi_fid,comedi_sample_size*(videosync_channel-1),'bof');
        comedi_size = file_size(comedi_file);
        videosync_file = [recording_filename_root '.comedivideosync.txt'];
        videosync_fid = fopen(videosync_file, 'wt');
        
        iCode = 0;  videosync_data= [];
        chk = 1;
        N = 0;
        while(chk)
            N = N+1;
%            disp(['COMEDI: Loop ' num2str(N)]);
            comedi = fread(comedi_fid, FS*T, [comedi_sample_format '=>single'], ...
                comedi_sample_size*(DECIMATION_FACTOR*CH-1));
            comedi = comedi-mean(comedi);

            stdcomedi = std(comedi);
            if stdcomedi > THRESHOLD & prev_stdcomedi < THRESHOLD
                iCode = iCode + 1;
              wintime = min(find(abs(comedi)>3e3));
              timestamp = (N-1)*T + wintime./1e4;  % Timestamp in seconds
              code = demodulate_timecode(comedi,1e4);
              disp(['Code detected: ' num2str(code)]);
              videosync_data(iCode,1) = timestamp;  
              videosync_data(iCode,2) = code;
              if code ~= prev_code + 1;
                %plot(comedi)
                %pause;
              end
              prev_code = code;
            end

            prev_stdcomedi = stdcomedi;

            if(~isempty(figh)),
               setappdata(figh, 'percent_done', ftell(comedi_fid)/comedi_size); 
            end
            
            if (length(comedi) < FS*T)
                chk = 0; 
            end
        end
        
        %  Fix isolated code errors
        a = diff(videosync_data(:,2));
        ind = find(a<1);  %  Timecodes that don't step by one.
        for iInd = 1:length(ind);
            if videosync_data(ind(iInd)+2,2) == videosync_data(ind(iInd),2)+2;
                videosync_data(ind(iInd)+1,2) = videosync_data(ind(iInd)+2,2)-1;
            end
        end
        
        %  Write to disk
        for iCode = 1:size(videosync_data,1)
            timestamp = videosync_data(iCode,1);
            code = videosync_data(iCode,2);
            fprintf(videosync_fid,'%6.3f\t%d\n',timestamp,code);
        end
        fclose(videosync_fid);
        fclose(comedi_fid);

    else
        ret = 0;
        disp('No comedi file found.  Skipping.')
    end
    ret = 1;
end
