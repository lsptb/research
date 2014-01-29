function [peaks,sel] = select_peakpairs(peaks,tau)
% peaks = select_peakpairs(peaks,tau)
% Purpose: return only pos & neg peaks separated by less than tau [sec]
% Inputs:
%   peaks - output of SO_detection()
%   tau   - max distance allowed between peak pairs
% Output:
%   peaks structure containing peak pairs
%
% Created by JSS on 18-Jun-2010
tic

if nargin < 2, tau = 1; end
Fs      = peaks(1).sfreq;
t       = peaks(1).tstart:1/Fs:peaks(1).tstop;
[sel(1:length(peaks)).pospeak] = deal([]);
[sel(1:length(peaks)).negpeak] = deal([]);
for ref = 1:length(peaks)
  pos   = peaks(ref).pospeak;
  neg   = peaks(ref).negpeak;
  npos  = length(pos);
  nneg  = length(neg);
  tpos  = t(pos);
  tneg  = t(neg);
  % look for tneg within tau of each tpos
  newpos = [];
  newneg = [];
  cnt   = 1;
  for k = 1:npos
    tk  = tpos(k);
    ind = nearest(tneg,tk);
    tmp = tneg(ind);
    if isempty(tmp), continue; end
    if abs(tk-tmp)<=tau && ~ismember(pos(k),newpos) && ~ismember(neg(ind),newneg)
      sel(ref).pospeak(cnt) = k;
      sel(ref).negpeak(cnt) = ind;
      newpos(cnt) = pos(k);
      newneg(cnt) = neg(ind);
      cnt        = cnt + 1;
    end
  end
  peaks(ref).pospeak = newpos;
  peaks(ref).negpeak = newneg;
end

toc