% addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
% addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/scripts
% 
% cd /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s2
% load s2_SO_flip_matrix_filt0.01-4Hz_toi0-3609.31_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_pairedpeaks_05-Jul-2010.mat
load s2_SO_flip_matrix_filt0.01-4Hz_toi0-8855.14_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_07-Jul-2010.mat


%% Method 1
%  consider consistency of relative polarity b/w refchan-locked averages

th  = 90; % threshold
clear badchans x D n per LHS RHS eqs badgrad1 badgrad2 badchanlabels pktype ch
for pktype = 1:2 % pktype = 2;
  x   = flip(pktype).matrix;
  D   = diag(x);
  n   = length(D);
  per = zeros(n,1);
  for ch = 1:n
    LHS = D==x(:,ch);
    RHS = D(ch)==x(ch,:);
    eqs = LHS==RHS';
    per(ch) = 100*sum(eqs)/length(D);
  end
  badchans{pktype}  = match_str({flip(pktype).sensor_info.label},{flip(pktype).sensor_info(per<th).label});
  figure; highlight = badchans{pktype}; dewar; view(0,90); title(sprintf('%s: threshold @ %g%',flip(pktype).peaktype,th))
end
badchans = intersect(badchans{1},badchans{2});
badgrad1 = strmatch('grad1',{flip(pktype).sensor_info(badchans).typestring}); nbadgrad1 = length(badgrad1); sen1 = {flip(1).sensor_info(badgrad1).label};
badgrad2 = strmatch('grad2',{flip(pktype).sensor_info(badchans).typestring}); nbadgrad2 = length(badgrad2); sen2 = {flip(1).sensor_info(badgrad2).label};
fprintf('Number of bad grad1 sensors: %g (%s)\nNumber of bad grad2 sensors: %g (%s)\n',nbadgrad1,[sen1{:}],nbadgrad2,[sen2{:}]);

badchanlabels = {flip(1).sensor_info(badchans).label};



%% Method 2
%  consider mean cross-chan sum after channel-dependent normalizations

pktype     = 1; % pospeak
for pktype = 1:2
  flipmat = flip(pktype).matrix;
  nchan   = size(flipmat,1);
  rowsums = zeros(nchan,1);
  medians = zeros(nchan,1);
  medper  = zeros(nchan,1);
  for k = 1:nchan
    thisrow = flipmat(k,:);
    flipix  = find(thisrow==-1);
    thismat = flipmat;
    thismat(:,flipix) = -thismat(:,flipix);   % flip if -1 in this row
    rowsums(:,k)      = sum(thismat,2);       % sum across row
    medians(:,k)      = median(thismat,2);
    med1          = thismat == repmat(medians(:,k),1,nchan);
    medper(:,k)   = 100*sum(med1,2) / nchan;
  end
  totsums = sum(abs(rowsums),1);
  figure; 
  subplot(2,1,1),plot(totsums,'.'); axis tight; hline(0,'k');
  subplot(2,1,2),plot(mean(medper,1),'.'); axis tight; title('mean median'); xlabel('chan');
  highlight = find(totsums == max(totsums)); figure; dewar; view(0,90);
  % [cidx,ctrs] = kmeans(totsums,8);
  % highlight = find(cidx==cidx(ctrs==min(ctrs))); figure; dewar; view(0,90);

  tmp       = mean(medper,1);
  highlight = find(tmp > 65); figure; dewar; view(0,90);
  if pktype == 1, poslist = {flip(1).sensor_info(highlight).label}'; end
  if pktype == 2, neglist = {flip(2).sensor_info(highlight).label}'; end  

nbin  = 50;
[N,X] = hist(tmp,nbin);
[cidx,ctrs] = kmeans(N.^2,2);
subplot(2,1,1),plot(X(cidx==1),N(cidx==1),'b.',X(cidx==2),N(cidx==2),'r.')
% subplot(2,1,1),plot(X,N,'.-')
N     = smooth(N,round(sqrt(nbin)),'lowess');
[cidx,ctrs] = kmeans(N.^2,2);
subplot(2,1,2),plot(X(cidx==1),N(cidx==1),'b.',X(cidx==2),N(cidx==2),'r.')





