% ntools_procGenerateNicoletSync(recording_filename_root)
%  generates Nicolet sync file
%
%  Inputs: recording_filename_root - recording root without the suffixes.

function ret=ntools_procGenerateNicoletSync(recording_filename_root)

    % neural sampling rate and period
    [nFs nCh] = ntools_hdrxml_read([recording_filename_root '.low.nspike.hdrxml']);
    nT = 1/nFs;
    nSize = 2;
    
    % video framerate and period
    vFs = 29.970;
    vT = 1/vFs;
    
    % removes suffix if necess
    recording_filename_root = make_filename_root(recording_filename_root); 

    comedi_ts_file = [recording_filename_root '.comedivideosync.txt'];
    video_ts_file = [recording_filename_root '.wavvideosync.txt'];

    eval(['comedi_ts = ', 'load(''', comedi_ts_file, ''');']);
    eval(['video_ts = ', 'load(''', video_ts_file, ''');']);
    
    % get intersection of timestamps
    % find timecodes shared between video/comedi
    shared_ts = intersect(comedi_ts(:,2), video_ts(:,2));
    ts(:,1) = shared_ts;

    % extract time values for shared indices
    [x, x, ind] = intersect(shared_ts, comedi_ts(:,2)); 
    ts(:,2) = comedi_ts(ind);
    [x, x, ind] = intersect(shared_ts, video_ts(:,2));
    ts(:,3) = video_ts(ind);


    % t_v - time elapsed in the video track from start of recording to 
    %       a given marker
    % t_c - time elapsed in the comedi track from start of recording to 
    %       a given marker
    % calculate t_v - t_c where t_v and t_c are aligned vectors of timestamps, 
    % giving us a vector of time differences.  when comedi leads video, 
    % these qtys are +ive, when the video track leads comedi, these qtys 
    % are -ive.
    ts(:,4) = ts(:,3) - ts(:,2);
    %ts(:,5) = vertcat(ts(1,2) - ts(1,3), diff(ts(:,4)))
    %ts(:,6) = round(ts(:,5)*vFs);
    %ts(:,7) = vertcat(0, diff(ts(:,6)));
    
    % get length of neural data (this could probable use the headers...)
    stat = dir([recording_filename_root '.low.nspike.dat']);
    d_length = stat.bytes/nSize/nCh;
   
    % the edf converter truncates to the last whole second
    d_length = floor(d_length/nFs)*nFs;

    fid = fopen([recording_filename_root '.001.sync'],'w');
    
    f = 1;
    for i = 1:d_length
        t = i*nT;
       
        % find the index of the closest marker to the current value of t 
        [x, ts_index] =  min(abs(ts(:,2) - t));

        % add the time offset from the closest marker this gives us time in
        % in terms of "video time" 
        adj_t = t + ts(ts_index,4);
       
        % if we're at the start and video leads comedi (ie: no video for first
        % data, we'll have a negative time value so hardcode to first frame) 
        if (adj_t < 0)
            f = 1;
        else
            % convert video time to video frame
            f = round(adj_t*vFs);
        end
        
        % this script is really slow because matlab fprintf really sucks
        % if speed is a problem, replace this with a mex function.  i'm
        % leaving it as-is for now as that adds another moving part and
        % performance has not yet shown to be a problem.
        fprintf(fid, '%d\t%d\r\n', i, f);
        
        if mod(t,1) == 0
            display(['Writing sync file, t=' num2str(t)]);
        end

    end
    fclose(fid);
    ret=1;
end
    
    
    
