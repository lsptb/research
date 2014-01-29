%
% ntools_cleandecieeg_load_int(filename, start_samp, stop_samp)
%
% loads data from a given cleandecieeg datafile starting with start_samp (inclusive)
% and ending at end_samp (inclusive) into a matlab matrix in [channel,sample]
% format.
%
% example: data = ntools_cleandecieeg_load_int('rec001.cleandecieeg.dat', 100, 200)

function data = ntools_cleandecieeg_load_int(varargin)
    data = ntools_load_int('cleandecieeg', varargin{:});
end
