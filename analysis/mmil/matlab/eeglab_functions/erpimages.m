% erpimages() - image multiple single-trial ERPs optionally sorted on 
%               an input variable and smoothed by moving-average 
%               Plot selected single channels in separate panels.
% Usage:
%   >> [outdata,outvar,outx] = erpimages(data,chans,rt,times,'title',...
%                                           avewidth,decimate,option(s));
%
% Inputs:
%   data     - (channels,ntrials*frames) data matrix
%   chans    - vector of channel numbers to plot in multiple panels
%   sortvar  - vector variable to sort trials on (ntrials = length(sortvar))
%   times    - vector of times (frames = length(times)) {def|0->[0:frames-1]}
%  'title'   - title string {default none}
%   avewidth - ntrials to perform moving average on (can be non-int) {def|0->1}
%   decimate - ntrials to decimate output by (can be non-int) {def|0->1}
%
% Optional Inputs:
%  'noplot'  - don't plot sortvar {default: plot if in times range}
%  'nosort'  - don't sort data on sortvar {default: sort}
%  'caxis'   - [lo hi] -> set color axis limits {default: data bounds}
%  'weights' - [weights] ICA weight matrix
%  'comps'   - [component_numbers] add contour lines for specified components
%
% Outputs:
%   outdata  - (times,pointsout) data matrix (after moving average)
%   outvar   - (1,pointsout) sortvar vector (after moving average)
%   outx     - (1,pointsout) smoothed trial numbers (after moving average)
%
% Authors: Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, La Jolla, 3/1998 
%
% See also: erpimage()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Tzyy-Ping Jung & Scott Makeig, CNL / Salk Institute, La Jolla 3-2-98
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

% $Log: erpimages.m,v $
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 3/5/98 added nosort option -sm
% 01-25-02 reformated help & license, added links -ad 

function [data,sortvar,outx] = erpimage(data,chans,sortvar,times,titl,avewidth,decfactor,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)

YES = 1;  % logical argument
NO  = 0;
DEFAULT_NOSHOW = NO; % show sortvar by default
DEFAULT_NOSORT = NO; % sort on sortvar by default
DEFAULT_DECIMATION_FACTOR = NO; % -> 1
DEFAULT_SMOOTHING_WIDTH = NO;   % -> 1
LABELFONT = 16;
TICKFONT  = 14;

if nargin < 4
  help erpimage
  return
end

if size(data,1) < max(chans)
  help erpimage
  fprintf('erpimage(): data and chans args do not agree\n');
  return
end

framestot = length(data);
ntrials = length(sortvar);
if ntrials < 2
  help erpimage
  return
end

frames = floor(framestot/ntrials);
if frames*ntrials ~= framestot
  fprintf('erpimage(); length of sortvar does not divide data length.\n')
  help erpimage
  return
end

noshow = DEFAULT_NOSHOW;
nosort = DEFAULT_NOSORT;

Caxis = [];
base = [];
icaweights = [];
icacomps = [];
Caxflag  = NO;
baseflag = NO;
icaflag  = NO;
compflag = NO;
if nargin > 7
  for n=8:nargin
    eval(['arg = arg' int2str(n-7) ';']);
    if Caxflag == YES
      Caxis = arg;
      if size(Caxis) ~= [1 2]
        help erpimage
        return
      end
      Caxflag == NO;
    elseif baseflag == YES
      base = arg;
      baseflag = NO;
    elseif icaflag  == YES
      icaweights = arg;
      icaflag = NO;
    elseif compflag  == YES
      icacomps = arg;
      compflag = NO;
    elseif strcmp(arg,'base') | strcmp(arg,'baseline')
      baseflag = YES;
    elseif strcmp(arg,'noshow')
      noshow = YES;
    elseif strcmp(arg,'nosort')
      nosort = YES;
    elseif strcmp(arg,'caxis')
      Caxflag = YES;
    elseif strcmp(arg,'weights')
      icaflag = YES;
    elseif strcmp(arg,'comps')
      compflag = YES;
      fprintf(' comps argument not yet implemented... \n');
    end
  end
end

if nargin < 7
  decfactor = DEFAULT_DECIMATION_FACTOR;
end
if nargin < 6
  avewidth = DEFAULT_SMOOTHING_WIDTH;
end
if nargin < 5
  titl = '';
end
if nargin < 4
  times = NO;
end
if times ==NO,
   times = NO:frames-1;
end

if avewidth == NO,
  avewidth = 1;
end
if decfactor == NO,
  decfactor = 1;
end
if avewidth < 1
  help erpimage
  fprintf('Variable avewidth cannot be < 1.\n')
  return
end
if avewidth > ntrials
  fprintf('Setting variable avewidth to max %d.\n',ntrials)
  avewidth = ntrials;  
end
if decfactor < 1
  help erpimage
  fprintf('Variable decfactor cannot be < 1.\n')
  return
end
if decfactor > ntrials
  fprintf('Setting variable decfactor to max %d.\n',ntrials)
  decfactor = ntrials;  
end

if ~isempty(icaweights) 
  act = weights*data;
  if size(weights,1) == size(weights,2)
     winv = inv(weights);
  else
     winv= pinv(weights);
  end
end

clf
npanes = length(chans);
ht = floor(sqrt(npanes));
wdth=ht;
while wdth*ht<npanes
    ht = ht+1;
end
   
for s = 1:npanes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chan = chans(s);
fprintf('Chan %d:\n',chan);
dat=reshape(data(chan,:),frames,ntrials);
fprintf('Plotting input data as %d epochs of %d frames.\n',...
                             ntrials,frames);
if nosort == YES
  fprintf('Not sorting data on input sortvar.\n');
else
  fprintf('Sorting data on input sortvar.\n');
if s==1
  [sortvar,ix] = sort(sortvar);
end
  dat = dat(:,ix);
end

if avewidth > 1 | decfactor > 1
  if nosort == YES
    fprintf('Smoothing the data using a window width of %g epochs\n',avewidth);
  else
    fprintf('Smoothing the sorted data using a window width of %g epochs\n',...
                       avewidth);
  end
  fprintf('  and a decimation factor of %g\n',decfactor);
  [dat,outx] = movav(dat,1:ntrials,avewidth,decfactor); % use square window
if s==1
  [sortvar,outx] = movav(sortvar,1:ntrials,avewidth,decfactor); 
end
  fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                          frames,length(outx));
else
  outx = 1:ntrials;
end

if ~isempty(Caxis)
  mindat = Caxis(1);
  maxdat = Caxis(2);
  fprintf('Using the specified caxis range of [%g,%g].\n',mindat,maxdat);
else
  mindat = min(min(dat)); % dat is single-channel
  maxdat = max(max(dat));

  maxdat =  max(abs([mindat maxdat])); % make symmetrical about 0
  mindat = -max(abs([mindat maxdat]));
  fprintf('The caxis range will be the data range in the specified channels [%g,%g].\n',...
                   mindat,maxdat);
end
subplot(ht,wdth,s)
if ~isempty(base)
  size(base)
  basevec = mean(dat(base,:));
  dat = dat - ones(length(times),1)*basevec;
end
imagesc(outx,times,dat,[mindat,maxdat]); 
                                   % imagesc() = image() with scaling

if ~isempty(icaweights) 
  proj = zeros(size(dat));
  for c=icacomps
     proj = proj + winv(c,chan)*act(c,:);
  end
  proj = reshape(proj,length(times),ntrials);
  if avewidth > 1 | decfactor > 1
    [proj,outx] = movav(proj,1:ntrials,avewidth,decfactor); % use square window
  end
  hold on
  cval = 0.75*std(abs(dat))*sign(max(dat)); 
  contour(outx,times,dat,cval,'k','LineWidth',2);
end

if noshow == NO & (min(sortvar) < min(times) | max(sortvar) > max(times))
  fprintf('Variable sortvar not in range of times: will not be plotted.\n');
  noshow = YES;
end

if noshow == YES
  fprintf('Not overplotting sorted sortvar on data.\n');
else
  fprintf('Overplotting sorted sortvar on data.\n');
  e = 0:ntrials-1;
  hold on; plot(outx,sortvar,'k','LineWidth',2); hold off % overplot sortvar
end

if s/wdth>ht-1
 if s==npanes
   if nosort == YES
     l=xlabel('Epoch Number');
   else
     l=xlabel('Sorted Epoch Number');
   end
   set(l,'FontSize',LABELFONT);
 end
else
 set(gca,'XTickLabels',[]);
end
if rem(s,wdth)==1
 if s==1
   l=ylabel('Time (msec)');
   set(l,'FontSize',LABELFONT);
 end
else
 set(gca,'YTickLabels',[]);
end

ttl = [titl ', chan ' int2str(chan)];
t=title(ttl);
set(t,'FontSize',LABELFONT);
set(gca,'FontSize',TICKFONT)

% if s == npanes
  colormap(jet(256))
  caxis([mindat maxdat]);
  cax = colorbar('vertical')
  set(cax,'FontSize',TICKFONT)
% end

fprintf('\n');
end % subplot s
