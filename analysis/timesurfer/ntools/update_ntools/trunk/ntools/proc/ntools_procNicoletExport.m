% ntools_procNicoletExport(recording_filename_root)
%
% Main entry point for Nicolet exporter.
%
% Creates a new directory underneath the recording root named NicoletExport
% containing a copy of the avi, a Nicolet sync file and the neural data
% converted to EDF format.  (Note: The neural data is truncated to the
% last whole second.)
%
% This directory is in the format expected by the Nicolet NSpike module,
% and should be transferred in full.
%
%
%  Inputs: recording_filename_root - recording root without the suffixes.
%
% E.g.:
%  recording_filename_root = '/mnt/esata/NY140/SP/081030/rec005';
%   
%  
% NOTES:
%
% This process is not instantaneous.
%
% The following packages must be installed for this to work:
%
% mplayer - (for demuxing the AVI file)
% libxml2, libxml2-dev (for supporting the C EDF converter)
% 
% additionally, the edfconverter in nsuite/ntools/conv/edfconvert
% must be compiled the edfconvert_binary_path line must be modified below
% 
function ntools_procNicoletExport(recording_filename_root)
    global NYUMCDIR;

    % this sucks, but it will have to do for now as we don't really have
    % any standard place for binaries used in conversion/processing
    % XXX must be edited when the tools are installed on a new machine
    % edfconvert must also be compiled
    % ***************** MODIFY THE FOLLOWING LINE ***************************
    edfconvert_binary_path = [NYUMCDIR '/svn_code/adam/nsuite/ntools/conv/edfconvert/edfconvert'];
 %tic
    recording_filename_root = make_filename_root(recording_filename_root); %removes suffix if necess
    
    figh = get(0, 'CurrentFigure'); %use setappdata to report percent done if in GUI

    % check existence of avi file
    avi_file = [recording_filename_root '.avi'];
    if isempty(dir(avi_file))
        error(['AVI file missing: ' avi_file]);
    end

    % create output directory
    nicolet_export_directory = [fileparts(recording_filename_root) '/NicoletExport'];
    disp(['Making Nicolet export directory: ' nicolet_export_directory]);
    if ~mkdir(nicolet_export_directory);
        error('Unable to create export directory');
    end
    
    % generate edf file (edfconvert)    
    edfconvert_cmd = [edfconvert_binary_path ' ' recording_filename_root ... 
        '.low.nspike.hdrxml'];
    disp(['Converting to EDF format: ' edfconvert_cmd]);
    if system(edfconvert_cmd) < 0
        error('EDF conversion failed');
    end
    
    % move edf file to output directory
    short_filename = get_last_path_element(recording_filename_root);
    source_edf_file = [recording_filename_root '.low.nspike.edf'];
    dest_edf_file = [nicolet_export_directory '/' short_filename '.001.edf'];
    if ~movefile(source_edf_file, dest_edf_file)
        error(['Unable to move EDF file to output directory: src=' source_edf_file ...
            ' dest=' dest_edf_file]);
    end
    
    % demux avi file with mplayer
    demux_cmd = ['mplayer -vo null -vc null -ao pcm:file=' ...
      recording_filename_root '.wav ' avi_file];
    disp(['Extracting audio tracks from video: ' demux_cmd]);
    if system(demux_cmd) < 0
        error('Audio extraction failed');
    end
    
    % run procWavVideoSync
    disp('Extracting and decoding time markers from AVI audio');
    if ~ntools_procWavVideoSync(recording_filename_root,2);
        error('AVI time marker decoding failed');
    end
    
    % run procComediVideoSync
    disp('Extracting and decoding time markers from NI/Comedi recording');
    if ~ntools_procComediVideoSync(recording_filename_root);
        error('NI/Comedi time marker decoding failed');
    end
    
    % generate sync file
    disp('Generating Nicolet synchronization file');
    if ~ntools_procGenerateNicoletSync(recording_filename_root)
        error('Synchronization file generation failed.');
    end
    
    % move sync file to output directory
    dest_sync_file = [nicolet_export_directory '/' short_filename '.001.sync'];
    sync_file = [recording_filename_root '.001.sync'];
    disp(['Moving ' sync_file ' to ' dest_sync_file]);
    if ~movefile(sync_file, dest_sync_file)
        error(['Unable to move sync file to output directory: src=' sync_file ...
            ' dest=' dest_sync_file]);
    end
        
    % copy video to output directory
    dest_avi_file = [nicolet_export_directory '/' short_filename '.001.avi'];
    disp(['Copying ' avi_file ' to ' dest_avi_file]);
    if ~copyfile(avi_file, dest_avi_file)
        error('File copy failed');
    end
     
    %toc
    % ?? 2gb filesize limit?  will this be a problem?
    
    
    function ret=get_last_path_element(path)
        
        [ret,r] = strtok(path,'/');
        
        if isempty(r)
            return;
        end
        
        while ~isempty(r)
            [ret,r] = strtok(r,'/');
        end

    end

end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
