function [cmd,loadflag,saveflag,funclass] = get_funcalls(varargin)
% format 1: [cmd,loadflag,saveflag,funclass] = get_funcalls(name,itype,otype,iflag,oflag)
% format 2: [cmd,loadflag,saveflag,funclass] = get_funcalls(parms.function(f))
% Columns in funtypes.csv:
% 1-ID
% 2-iflag
% 3-oflag
% 4-loadflag
% 5-saveflag
% 6-instring
% 7-outstring
% 8-funclass
%
% Created by Jason Sherfey on 13-Apr-2009
% 

if nargin==1
  parms = varargin{1};
  fun   = parms.name;
  itype = parms.input;
  otype = parms.output;
  iflag = parms.control_flags.iflag;
  oflag = parms.control_flags.oflag;
elseif nargin==5
  fun = varargin{1};
  itype = varargin{2};
  otype = varargin{3};
  iflag = varargin{4};
  oflag = varargin{5};  
end
[funtypes  ext] = find_masterfile('funtypes');
[data, result]  = mmil_readtext(funtypes, ',','','','textual-empty2NaN');

% find loadflag & saveflag given iflag & oflag
ind1 = find(ismember(data(:,2),num2str(iflag)));
ind2 = find(ismember(data(:,3),num2str(oflag)));
ind  = intersect(ind1,ind2);
if length(ind) > 1
  error('%s: redundant function definition in funtypes.csv\n',mfilename);
end

loadflag = str2num(data{ind,4});
saveflag = str2num(data{ind,5});
funclass = data{ind,8};

istr = data{ind,6};
istr(regexp(istr,'\w\s\w')+1)=',';

ostr = data{ind,7};
if strcmp(ostr,'NaN') || isempty(ostr) || ~ischar(ostr), ostr = ''; end

if length(regexp(istr,';'))>1
  ind = regexp(istr,';');
  ind = ind(end-1);
  cmd = [istr(1:ind+1) ostr istr(ind+2:end)];
else
  cmd = [ostr ' ' istr];
end

cmd = strrep(cmd,'fun',fun);
if ischar(itype), cmd = strrep(cmd,'itype',itype); end
if ischar(otype), cmd = strrep(cmd,'otype',otype); end
  
  





