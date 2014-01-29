%
% ntools_audio_load_int(filename, start_samp, stop_samp)
%
% loads data from a given audio datafile starting with start_samp (inclusive)
% and ending at end_samp (inclusive) into a matlab matrix in [channel,sample]
% format.
%
% example: data = ntools_audio_load_int('rec001.audio.dat', 100, 200)

function data = ntools_audio_load_int(varargin)
    error('This will not work at the moment! Does anything even call this?');
    data = ntools_load_samples('audio', varargin{:});
end
