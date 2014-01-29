function r = isempty2(o, ascell)
% r = isempty2(o, [ascell])
%
% Purpose:
%   Behaves as isempty, but
%   has an option to check cell array CONTENTS
%   instead of cell array itself, for empty values.
%
% Input Arguments:
%   o      : input object (can be any type allowable by iscell)
%   ascell : boolean specify whether to evaluate o as a cell
%            {default: false}
%
% Output Arguments:
%   r      : if o is not a cell, or if ascell is not true, same as iscell.
%            otherwise, r is matrix with the same shape as o, 
%            specifying whether each cell is empty.
%
% See also: isempty
%
% Created By:       Ben Cipollini on 06/15/2007
% Last Modified By: Ben Cipollini on 08/09/2007

  % Ascell != true, just do normal
  if (~exist('ascell', 'var') || ~ascell)
    r =  isempty(o);
    
  % Not a cell, so do normal
  elseif (~iscell(o))
    r = isempty(o);

  % Is cell; check each cell's contents for 'isempty',
  % and return the cell indices that are actually empty.
  elseif (isempty(o))
    r = isempty(o);
    
  else
    o_sz = size(o);
    o = reshape(o, [prod(size(o)) 1]); % make into a cell array
    r = [];
    for i=1:length(o)
        r(i) = isempty2(o{i},true);
    end;
    r = reshape(r, o_sz); % put back into original form
  end;
