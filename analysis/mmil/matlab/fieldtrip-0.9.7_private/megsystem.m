function [type] = megsystem(grad)

% MEGSYSTEM returns a string that describes the manufacturer of the
% system and the type of MEG sensors. The heuristic approach is to test
% the gradiometer definition on a few ascpects and score points for each
% of them. The system that is most "similar" wins.
% 
% Use as
%   [str] = megsystem(grad)
% where the input is a gradiometer structure and the output will be any of
%   'ctf151'
%   'ctf275'
%   'ctf151_planar'
%   'ctf275_planar'
%   'neuromag122'
%   'neuromag306'

% Copyright (C) 2004-2006, Robert Oostenveld
%
% $Log: megsystem.m,v $
% Revision 1.5  2006/01/30 14:06:04  roboos
% added square brackets around output variable in function definition
% added copyrights and log
% cleaned up help documentation
%

description = {
  'ctf151'
  'ctf275'
  'ctf151_planar'
  'ctf275_planar'
  'neuromag122'
  'neuromag306'
};

% start with an empty counter for each system
similar = zeros(size(description));

% look at the number of channels
Nchan = length(grad.label);
similar(1) = similar(1) + (abs(Nchan-151)   <  20);
similar(2) = similar(2) + (abs(Nchan-275)   <  20);
similar(3) = similar(3) + (abs(Nchan-151*2) <  20);
similar(4) = similar(4) + (abs(Nchan-275*2) <  20);
similar(5) = similar(5) + (abs(Nchan-122)   <  20);
similar(6) = similar(6) + (abs(Nchan-306)   <  20);
similar(1) = similar(1) + (abs(Nchan-151)   <= abs(Nchan-275));
similar(2) = similar(2) + (abs(Nchan-275)   <= abs(Nchan-151));
similar(3) = similar(3) + (abs(Nchan-151*2) <= abs(Nchan-2*275));
similar(4) = similar(4) + (abs(Nchan-275*2) <= abs(Nchan-2*151));
similar(5) = similar(5) + (abs(Nchan-122)   <= abs(Nchan-306));
similar(6) = similar(6) + (abs(Nchan-306)   <= abs(Nchan-122));

% look at the beginning of Neuromag channel names
similar(5) = similar(5) + length(find(strmatch('MEG', grad.label)));
similar(6) = similar(6) + length(find(strmatch('MEG', grad.label)));

% look at the beginning of CTF channel names
similar(1) = similar(1) + length(find(strmatch('MZ', grad.label)));
similar(2) = similar(2) + length(find(strmatch('MZ', grad.label)));
similar(3) = similar(3) + length(find(strmatch('MZ', grad.label)));
similar(4) = similar(4) + length(find(strmatch('MZ', grad.label)));

% look at whether the names contain _dH and _dV
similar(3) = similar(3) + length(find(a_strmatch('_dH', grad.label)));
similar(3) = similar(3) + length(find(a_strmatch('_dV', grad.label)));
similar(4) = similar(4) + length(find(a_strmatch('_dH', grad.label)));
similar(4) = similar(4) + length(find(a_strmatch('_dV', grad.label)));

% if they do not contain any _dH or _dV, the planar CTF systems is less likely that their axial counterpart
similar(3) = similar(3) - (length(find(a_strmatch('_dH', grad.label)))==0);
similar(3) = similar(3) - (length(find(a_strmatch('_dV', grad.label)))==0);
similar(4) = similar(4) - (length(find(a_strmatch('_dH', grad.label)))==0);
similar(4) = similar(4) - (length(find(a_strmatch('_dV', grad.label)))==0);

% determine to which MEG ssytem the input data is the most similar
[m, i] = max(similar);
if m==0 || length(find(similar==m))>1
  error('could not detect the type of MEG system');
else
  type = description{i};
end

% search for substring in each element of a cell-array
function a = a_strmatch(str, strs)
a = zeros(length(strs),1);
for i=1:length(strs)
  a(i) = ~isempty(strfind(strs{i}, str));
end
a = find(a);

