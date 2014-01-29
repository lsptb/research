function [res,peaks] = find_SlowOscillations(data)
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts

% MEG data
% lay      = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
% badchans = '/home/jsherfey/projects/sleep/badchans/s8_badchans.txt';
% matfiles = {...
% '/home/halgdev/projects/jsherfey/sleep/s8/matfiles/SL_1_nb01_060808_grad1.mat' ...
% '/home/halgdev/projects/jsherfey/sleep/s8/matfiles/SL_2_nb01_060808_grad1.mat' ...  
% '/home/halgdev/projects/jsherfey/sleep/s8/matfiles/SL_3_nb01_060808_grad1.mat' ...
% '/home/halgdev/projects/jsherfey/sleep/s8/matfiles/SL_4_nb01_060808_grad1.mat' ...
% '/home/halgdev/projects/jsherfey/sleep/s8/matfiles/SL_5_nb01_060808_grad1.mat' ...
% '/home/halgdev/projects/jsherfey/sleep/s8/matfiles/SL_6_nb01_060808_grad1.mat' ...
% };
% 
% % iEEG data
% badchans = '/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg/MG19/MG19_badchans.txt';
% matfiles = '/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg/MG19/MG19_sleep1.mat';
% 
% data = SO_combine_matfiles(matfiles); 
% data = ts_data_selection(data,'badchanfile',badchans);

% SO
T2T  = 1000;
fc   = [.3 2];
chan = [1:data.num_sensors];% chan=51:53;
% chan = 53:59;
dn   = ceil(T2T * data.sfreq / 1000);
thsh = fc(1);

cnt = 0;
res = [];
for ch = chan
  cnt  = cnt + 1;  
  proc = ts_data_selection(data,'channels',ch);
  proc = ts_preproc(proc,'bpfilter','yes','bpfreq',fc,'blc','no');
  
  x    = proc.epochs.data;
  t    = proc.epochs.time;

  % calculate instantaneous phase and frequency
  Xa   = hilbert(x); % analytic signal
  phi  = angle(Xa);  % instantaneous phase
  env  = abs(Xa);    % amplitude envelope
  frq  = (1/(2*pi)) * diff(unwrap(phi))./diff(t); % instantaneous frequency

  % troughs: phi > pi-pi/12 | phi < -pi+pi/12
  ix1    = find(phi > pi-pi/12);
  ix1    = ix1(diff(ix1) > 1);
  ix1(phi(ix1+1) > -pi+pi/12) = [];
  ix2    = ix1+1;
  ii     = bsxfun(@lt,pi-phi(ix1),abs(-pi-phi(ix2)));
  negix  = sort([ix1(ii) ix2(~ii)]);
  nkeep  = 1:length(negix);
  allneg = negix;

  % peaks (phi zero-crossings): phi < pi/12 & phi > -pi/12
  ix     = find((phi > -pi/12) & (phi < pi/12));
  if max(ix) == length(phi), ix = ix(1:end-1); end
  if min(ix) == 1, ix = ix(2:end); end
  ix     = ix(phi(ix)<0 & phi(ix+1)>0);
  ii     = bsxfun(@lt,abs(phi(ix)),phi(ix+1));
  posix  = sort([ix(ii) ix(~ii)+1]);
  pkeep  = 1:length(posix);
  allpos = posix;

  prej   = [];
  nrej   = [];
  % require monotonic phi between troughs: 
  % only keep n such that frq(t) > thresh (default: fc1) for all t=[t(n),t(n+1)]
  ind  = find(frq < thsh);
  ind  = ind(find(diff(ind)>dn)+1);
  for i = 1:length(ind)
      tmp = find(ind(i)-negix > 0);
      if ~isempty(tmp)
        nrej = [nrej tmp(end)]; 
      end
      if phi(ind(i)) <= 0
        tmp = find(ind(i)-posix < 0); 
        if ~isempty(tmp), tmp = tmp(1); end
      else
        tmp = find(ind(i)-posix > 0); 
        if ~isempty(tmp), tmp = tmp(end); end
      end
      if ~isempty(tmp)
        prej = [prej tmp]; 
      end
  end
  posix(prej) = [];
  negix(nrej) = [];
  pkeep(prej) = [];
  nkeep(nrej) = [];
  
  % amplitude criterion
  E     = median(x(posix));% + std(x(posix));
  ind   = x(posix) > E;
  posix = posix(ind);
  pkeep = pkeep(ind);
  E     = median(x(negix));% - std(x(negix));
  ind   = x(negix) < E;
  negix = negix(ind);
  nkeep = nkeep(ind);
  
  % get zero-crossings from waveform
  % cutoff = median(x(x>0))/10;
  u = posix; npos = length(u);
  v = negix; nneg = length(v);
  I = 1:length(x);
  N = I(x<0);
  P = I(x>0);
  Xlo = ones(1,npos);
  Xhi = ones(1,npos)*length(x);
  Ylo = ones(1,nneg);
  Yhi = ones(1,nneg)*length(x);
  % up
  tic
%   progress('init')
  for k = 1:npos
%     progress(k/npos)
    ix  = N(N<u(k));
    if isempty(ix)
      Xlo(k) = 1;
    else
      Xlo(k) = ix(end);
    end
    ix  = N(N>u(k));
    if isempty(ix)
      Xhi(k) = length(x);
    else
      Xhi(k) = ix(1);
    end
  %   disp(num2str(k))
  end
%   progress('close')
  toc
  % down
%   progress('init')
  for k = 1:length(v)
%     progress(k/nneg)
    ix  = P(P<v(k));
    if isempty(ix)
      Ylo(k) = 1;
    else
      Ylo(k) = ix(end);
    end
    ix  = P(P>v(k));
    if isempty(ix)
      Yhi(k) = length(x);
    else
      Yhi(k) = ix(1);
    end
  %   disp(num2str(k))
  end
%   progress('close')
  toc

  % collect results
  res(cnt).label        = proc.sensor_info.label;
  res(cnt).peak_index   = posix;
  res(cnt).trough_index = negix;
  res(cnt).peakstart_index   = Xlo;
  res(cnt).peakend_index     = Xhi;
  res(cnt).troughstart_index = Ylo;
  res(cnt).troughend_index   = Yhi;
  
  peaks(cnt).label = proc.sensor_info.label;
  peaks(cnt).time  = proc.epochs.time(posix);
  peaks(cnt).type  = ones(1,length(posix));
  peaks(cnt).time  = [peaks(cnt).time proc.epochs.time(negix)];
  peaks(cnt).type  = [peaks(cnt).type 2*ones(1,length(negix))];
end

% 
% log = '/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg/MG19/MG19_detection.log';
% fid = fopen(log,'w+');
% for ch = 1:length(res)
%   a = res(ch).peak_index;
%   b = res(ch).peakstart_index;
%   c = res(ch).peakend_index;
%   d = res(ch).trough_index;
%   e = res(ch).troughstart_index;
%   f = res(ch).troughend_index;
%   
%   ta = t(a);
%   tb = t(b);
%   tc = t(c);
%   td = t(d);
%   te = t(e);
%   tf = t(f);
%   xa = x(a);
%   xb = x(b);
%   xc = x(c);
%   xd = x(d);
%   xe = x(e);
%   xf = x(f);
%   data_a = data.epochs.data(ch,a);
%   data_b = data.epochs.data(ch,b);
%   data_c = data.epochs.data(ch,c);
%   data_d = data.epochs.data(ch,d);
%   data_e = data.epochs.data(ch,e);
%   data_f = data.epochs.data(ch,f);
%   % slopes
%   s1 = (xa-xb)./(ta-tb);
%   s2 = (xa-xc)./(tc-ta);
%   s3 = (xe-xd)./(td-te);
%   s4 = (xf-xd)./(tf-td);
%   s5 = (data_a-data_b)./(ta-tb);
%   s6 = (data_c-data_a)./(tc-ta);
%   s7 = (data_e-data_d)./(td-te);
%   s8 = (data_f-data_d)./(tf-td);
%   s12 = [s1' s2'];
%   s34 = [s3' s4'];
%   s1256 = [s1' s2' s5' s6'];
%   s56 = [s5' s6'];
%   s78 = [s7' s8'];
%   s3478 = [s3' s4' s7' s8'];
%   sfilt = [mean(s1) mean(s2) mean(s3) mean(s4) mean([s12(:)]) mean([s34(:)]) mean([s1256(:)])];
%   sraw  = [mean(s5) mean(s6) mean(s7) mean(s8) mean([s56(:)]) mean([s78(:)]) mean([s3478(:)])];
%   
% %   Xc  = t(res(ch).peak_index);
% %   Xlo = t(res(ch).peakstart_index);
% %   Xhi = t(res(ch).peakend_index);
% %   Yc  = t(res(ch).trough_index);
% %   Ylo = t(res(ch).troughstart_index);
% %   Yhi = t(res(ch).troughend_index);
% % 
% %   ampXc  = x(res(ch).peak_index);
% %   ampXlo = x(res(ch).peakstart_index);
% %   ampXhi = x(res(ch).peakend_index);
% %   ampYc  = x(res(ch).trough_index);
% %   ampYlo = x(res(ch).troughstart_index);
% %   ampYhi = x(res(ch).troughend_index);
%     
%   Xzz = Xhi-Xlo;
%   Yzz = Yhi-Ylo;
%   Wzz = [Xzz Yzz];
% 
%   % figure; hist(Xzz',50); set(gca,'xlim',[0 2])
%   % figure; hist(Yzz',50); set(gca,'xlim',[0 2])
%   figure; hist(Wzz',100); set(gca,'xlim',[0 2])
% 
%   [n1,x1] = hist(Xzz);
%   cntrd1  = sum(n1.*x1)/sum(n1);
%   [n2,x2] = hist(Yzz);
%   cntrd2  = sum(n2.*x2)/sum(n2);
%   [n3,x3] = hist(Wzz);
%   cntrd3  = sum(n3.*x3)/sum(n3);
%   fprintf(fid,'%s: COUNTS; CENTROID DURATION & FREQUENCY\n', res(ch).label);
%   fprintf(fid,'peak,    trough, both,   trough/peak (zero-to-zero)\n');
%   fprintf(fid,'%g     %g     %g     %g\n',[sum(n1) sum(n2) sum(n3) 100*sum(n2)/sum(n1)]);
%   fprintf(fid,'%g  %g  %g  %g\n',[2*[cntrd1 cntrd2 cntrd3]*1000 100*cntrd2/cntrd1]);
%   fprintf(fid,'%g  %g  %g  %g\n\n',[1./(2*[cntrd1 cntrd2 cntrd3]) 100*cntrd1/cntrd2]);
%   fprintf(fid,'%s: SLOPES\n',res(ch).label)
%   fprintf(fid,'filt: pos1,    pos2,    neg1,    neg2,    pos,    neg,    all)\n');
%   fprintf(fid,'      %g, %g, %g, %g, %g, %g, %g\n',sfilt);
%   fprintf(fid,'raw:  pos1,    pos2,    neg1,     neg2,     pos,     neg,     all)\n');
%   fprintf(fid,'      %g, %g, %g, %g, %g, %g, %g\n\n',sraw);  
%   fprintf(fid,'-----------------------------------------------------\n');
% end
% fclose(fid);
% disp('finished');
% 
% 


