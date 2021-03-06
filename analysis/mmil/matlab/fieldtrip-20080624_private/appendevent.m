function c = appendevent(a, b)

% APPENDEVENT

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: appendevent.m,v $
% Revision 1.1  2008/01/30 10:35:19  roboos
% moved from subfunction in read_event and write_event into seperate function, also used in bciagent.
% made duration and offset optional, just like timestamp
%

if isempty(a)
  c = b(:);
elseif isempty(b)
  c = a(:);
else
  c = a(:);
  for i=1:numel(b)
    c(end+1).type     = b(i).type;
    c(end  ).value    = b(i).value;
    c(end  ).sample   = b(i).sample;
    if isfield(b, 'timestamp')
      c(end  ).timestamp = b(i).timestamp; % optional
    end
    if isfield(b, 'offset')
      c(end  ).offset    = b(i).offset;    % optional
    end
    if isfield(b, 'duration')
      c(end  ).duration  = b(i).duration;  % optional
    end
  end
end