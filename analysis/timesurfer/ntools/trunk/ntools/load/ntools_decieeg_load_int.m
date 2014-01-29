%
% ntools_decieeg_load_int(filename, start_samp, stop_samp)
%
% loads data from a given decieeg datafile starting with start_samp (inclusive)
% and ending at end_samp (inclusive) into a matlab matrix in [channel,sample]
% format.
%
% example: data = ntools_decieeg_load_int('rec001.decieeg.dat', 100, 200)

function data = ntools_decieeg_load_int(varargin)
    data = ntools_load_int('decieeg', varargin{:});
end
