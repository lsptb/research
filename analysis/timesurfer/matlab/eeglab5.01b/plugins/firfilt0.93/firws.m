%firws() - Designs windowed sinc type I linear phase FIR filter
%
% Usage:
%   >> b = firws(m, f);
%   >> b = firws(m, f, w);
%   >> b = firws(m, f, t);
%   >> b = firws(m, f, w, t);
%
% Inputs:
%   m - filter order (mandatory even)
%   f - vector or scalar of cutoff frequency/ies (-6 dB;
%       pi rad / sample)
%
% Optional inputs:
%   w - vector of length m + 1 defining window {default blackman}
%   t - 'high' for highpass, 'stop' for bandstop filter {default low-/
%       bandpass}
%
% Output:
%   b - filter coefficients
%
% References:
%   Smith, S. W. (1999). The scientist and engineer's guide to digital
%   signal processing (2nd ed.). San Diego, CA: California Technical
%   Publishing.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   pop_firws, pop_firwsord, pop_kaiserbeta, windows

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: firws.m,v $
% Revision 1.1.1.1  2005/12/22 16:07:55  widmann
% initial import into CVS
%

function b = firws(m, f, t, w)

    if nargin < 2
        error('Not enough input arguments');
    end
    if length(m) > 1 || ~isnumeric(m) || ~isreal(m) || mod(m, 2) ~= 0 || m < 2
        error('Filter order must be a real, even, positive integer.');
    end
    f = f / 2;
    if f < 0 | f > 1
        error('Frequencies must fall in range between 0 and 1.');
    end
    if nargin < 3
        t = '';
        w = windows('blackman', (m + 1))';
    end
    if nargin < 4
        if ~ischar(t)
            w = t;
            t = '';
        else
            w = windows('blackman', (m + 1))';
        end
    end
    w = w(:)'; % Make window row vector

    b = fkernel(m, f(1), w);

    if length(f) == 1 && strcmpi(t, 'high')
        b = fspecinv(b);
    end

    if length(f) == 2
        b = b + fspecinv(fkernel(m, f(2), w));
        if isempty(t) || ~strcmpi(t, 'stop')
            b = fspecinv(b);
        end
    end

% Compute filter kernel
function b = fkernel(m, f, w)
    b = [sin(2 * pi * f * [-m / 2:-1]) ./ [-m / 2:-1] ...
        2 * pi * f ...
        sin(2 * pi * f * [1:m / 2]) ./ [1:m / 2]] .* w; % sinc * window
    b = b / sum(b); % Normalization to unity gain at DC

% Spectral inversion
function b = fspecinv(b)
    b = -b;
    b(1, (length(b) - 1) / 2 + 1) = b(1, (length(b) - 1) / 2 + 1) + 1;
