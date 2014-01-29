function [flip] = calc_flip_matrix(init_dat,init_peaks,varargin)
% Create flip matrix
% flip(i,j) = -1 (flip) or 1 (don't flip)
% [flip] = nrefchan x nchan
%   row i: whether to flip each channel j for a particular reference channel i
%   col j: whether to flip a particular channel j for each reference channel i
% flipped_data{ref} = data.*repmat(flip(ref,:)',1,size(data,2));

%%
%% MAY NOT BE WORKING PROPERLY (12-Aug-2010, JSS)
%%

parms = mmil_args2parms( varargin, ...
                         { 'blcwindow',[],[],...
                           'FlipWindow',.2,[],...
                           'EpochPad',.5,[],...
                           'toilim',[],[],...
                           'refindex',[],[],...
                           'blc','yes',[],...
                         }, ...
                         false );
                       
% data = SO_combine_matfiles(parms.matfiles);
% data = ts_data_selection(data,'toilim',parms.toilim);
if ~isfield(parms,'blcwindow'), parms.blcwindow = []; end
if ~isfield(parms,'FlipWindow'), parms.FlipWindow = .2; end
if ~isfield(parms,'EpochPad'), parms.EpochPad = .5; end
if ~isfield(parms,'toilim'), parms.toilim = []; end
t       = init_dat.epochs.time;
Fs      = init_dat.sfreq;
nchan   = init_dat.num_sensors;
flip_pad    = parms.FlipWindow / 2;
epoch_pad   = parms.EpochPad*1000;
I           = round(flip_pad*Fs);
[flip(1:2).sensor_info] = deal(init_dat.sensor_info);
[flip(1:2).matrix]      = deal(ones(nchan,nchan));
peaktime = init_peaks(1).tstart:1/init_peaks(1).sfreq:init_peaks(1).tstop;
if isfield(parms,'refindex') && ~isempty(parms.refindex)
  refindex = parms.refindex;
else
  refindex = 1:nchan;
end
tstart = tic;
fprintf('Calculating flip matrix:\n');
% Loop over peak types
peaktypes     = {'pospeak','negpeak'};
for type      = 1:length(peaktypes)
  refpeaktype = peaktypes{type};
  flip(type).peaktype = refpeaktype;
  flip(type).parms    = parms;
  % Loop over reference channels
  for rr = 1:length(refindex)
    ref = refindex(rr);
    fprintf('Reference %g of %g: %s (%s). Time elapsed: %g min\n',ref,nchan,init_peaks(ref).label,refpeaktype,toc(tstart)/60);
    refpeaks = init_peaks(ref).(refpeaktype);
    if ~isempty(parms.toilim)
      refpeaks = refpeaks(peaktime(refpeaks)>=parms.toilim(1)&peaktime(refpeaks)<=parms.toilim(2));
    end
    % Flip init_dat wrt polarity of ref-averaged SO waves over refpeaktype (all grads) => flip_dat
    % find simultaneous detections
    % - convert reference peaks to a cell array
    cellrefpeaks = num2cell(refpeaks);  % one element per refpeak
%     [cellepochs{1:nchan}] = deal([]);   % will have one element per channel
    % - for each chan, find the refpeaks with which it is invovled
    
% %     %% Option 1: Average chans on their own peaks near ref peaks
% %     bin = zeros(1,length(peaktime));
% %     ind = cellfun(@(x)[x-I:x+I],cellrefpeaks,'uniformoutput',false);
% %     ind = unique([ind{:}]);
% %     bin(ind) = 1;
% %     a = {init_peaks.pospeak};
% %     b = {init_peaks.negpeak};
% %     c = cellfun(@(x,y)[x y],a,b,'uniformoutput',false);
% %     for k = 1:nchan
% %       sel           = bin(c{k});
% %       cellepochs{k} = c{k}(sel==1);
% %     end
%     %% Option 2: Average chans on ref peaks near which they also have peaks
%     for k = 1:nchan
%       % combine pos & neg peaks in this channel & compare all to refpeaks
%       tk = []; tk = init_peaks(k).pospeak; tk = sort([tk init_peaks(k).negpeak]);
%       % find detections tk in k that occur near detections in cellrefpeaks
%       involved      = cellfun(@(y)(any((tk>=y-flip_pad*Fs)&(tk<=y+flip_pad*Fs))),cellrefpeaks);
%       cellepochs{k} = refpeaks(involved); % = indices to refpeaks with which chan k is involved
%     end    
    % % Option 3: Average chans on all ref peaks regardless of their own state
    [cellepochs{1:nchan}] = deal(refpeaks);
    % generate averages based on simultaneous detections    
    npeak           = cellfun(@length,cellepochs);
    epochdata       = SO_epochs(init_dat,cellepochs,epoch_pad);
    if strcmp(parms.blc,'yes'), epochdata = ts_preproc(epochdata,'blc',parms.blc,'verbose',0,'blcwindow',parms.blcwindow); end
    averagedata     = SO_average(epochdata,npeak);
    clear epochdata
    % create flip vector & flip_dat based on average over epochs from simultaneous detections
    tt        = averagedata.averages.time;
    tix       = nearest(tt,-flip_pad):nearest(tt,flip_pad); % flip window
    tmpdat    = averagedata.averages.data(:,tix);
    ref0      = averagedata.averages.data(ref,nearest(tt,0)); % reference value            
    refpol    = ref0 > 0; % 1 if pos, 0 if neg, reference polarity
    % determine whether to flip each channel for this reference
    for k  = 1:nchan
      % select window tk+/-fpad for this channel around {tk}ref
      tmp  = tmpdat(k,:);
      % find peaks in this window in this channel
      ix   = crossing(diff(tmp));
      % find the max peak
      amp  = max(abs(tmp(ix)));
      % threshold at 25% of the max
      thsh = .25*amp;
      % get indices to peaks above the threshold
      ix   = ix(abs(tmp(ix))>thsh);
      % find the remaining peak closest to t=0
      ix   = ix(nearest(tt(tix(ix)),0));
      x0   = averagedata.averages.data(k,tix(ix));
      xpol = x0 > 0;
      if xpol ~= refpol
        % opposite polarity => flip the channel
        flip(type).matrix(ref,k)  = -1;
      end
      clear tmp
    end
    clear tmpdat
  end
end
fprintf('done\n');