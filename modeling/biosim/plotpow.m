function [fig,pows,f]=plotpow(data,spec,varargin)
% Purpose: plot power and spectrogram
% Builds on: PowerSpecTA and SpecGramTA

parms = mmil_args2parms( varargin, ...
                   {  'plot_flag',1,[],...
                      'NFFT',[],[],...     % Wind, 8000
                      'WINDOW',[],[],...   % Wind, 8000
                      'NOVERLAP',[],[],... % WindOL, 7800
                      'FreqRange',[10 80],[],...
                      'NormAbs','Normalized',[],...
                      'Notch',[],[],...
                      'var',[],[],...
                      'varlabel','V',[],...
                      'spectrogram_flag',1,[],...
                   }, false);
fig=[];
var=parms.var;

t = data(1).epochs.time;
if isempty(parms.NFFT)
  NFFT=2^(nextpow2(length(t)-1)-1);
else
  NFFT=parms.NFFT;
end
if isempty(parms.WINDOW)
  WINDOW=NFFT;
else
  WINDOW=parms.WINDOW;
end
if isempty(parms.NOVERLAP)
  NOVERLAP = WINDOW-round(.02*NFFT);
else
  NOVERLAP = parms.NOVERLAP;
end

if ~isfield(spec,'entities') && isfield(spec,'cells')
  spec.entities=spec.cells;
end
npop = length(spec.entities);
Fs = fix(data(1).sfreq);

for pop=1:npop
  labels = {data(pop).sensor_info.label};
  if isempty(parms.var)
    var = find(~cellfun(@isempty,regexp(labels,['_' parms.varlabel '$'])),1,'first');
  else
    var = parms.var;
  end
  n = spec.entities(pop).multiplicity;
  t = data(pop).epochs.time;
  dat = double(squeeze(data(pop).epochs.data(var,:,1:n))');
  lfp = mean(dat,1)';
  res = PowerSpecTA(t,lfp,parms.FreqRange,WINDOW,parms.NormAbs,parms.Notch);
      % res.AreaPower (in FreqRange)
      % res.PeakPower
      % res.OscFreq
      % res.Pxx_HzPerBin = Bins/Fs
  f = res.f;
  if pop==1
    lfps = zeros(npop,length(t));
    pows = zeros(npop,length(f)); 
  end
  lfps(pop,:) = lfp;
  pows(pop,:) = res.Pxx;
  if parms.spectrogram_flag
    [yo,fo,to] = specmw(detrend(lfp),NFFT,Fs,WINDOW,NOVERLAP);
    %[yo,fo,to] = mtmspecTA2(detrend(dat),NFFT,Fs,WINDOW,NOVERLAP);
    %[yo,fo,to] = mtmspecTA(detrend(dat),NFFT,Fs,WINDOW,NOVERLAP,2.5);
    if pop==1
      tfpows = zeros(npop,length(to),length(fo));
    end
    tfpows(pop,:,:) = yo';
  end
end

% Fill screen with subplots:
nrows = npop;
ncols = 2;
nplots = nrows*ncols;

% ------------------------------
dx = 1 / ncols;
dy = 1 / nrows;
xstep = mod((1:ncols)-1,ncols);
ystep = mod((1:nrows)-1,nrows)+1;
xpos = .02+xstep*dx;
ypos = .02+1-ystep*dy;
plims=[0 200];
zlims=[-5 5];

if parms.plot_flag
  screensize = get(0,'screensize');
  fig=figure('position',screensize);
  set(gca,'Units','normalized','Position',[0 0 1 1]);
  for i = 1:nplots
    xi = mod(i-1,ncols)+1;
    yi = floor((i-1)./ncols)+1;
    subplot('Position',[xpos(xi) ypos(yi) .9/ncols .9/nrows]); set(gca,'units','normalized');
    pop = ceil(i/ncols); % index to this population
    labels = {data(pop).sensor_info.label};
    if isempty(parms.var)
      var = find(~cellfun(@isempty,regexp(labels,['_' parms.varlabel '$'])),1,'first');
    else
      var = parms.var;
    end    
    lab = data(pop).sensor_info(var).label;
    if mod(i,ncols)==1 % plot power spectrum
      plot(f,log10(pows(pop,:)));
      ylabel('power [normalized]');
      title([strrep(lab,'_','\_') ' power spectrum'],'fontsize',14,'fontweight','bold');    
      xlim(plims);
      if pop==npop
        xlabel('freq [Hz]'); 
      else
        set(gca,'xtick',[],'xticklabel',[]);
      end
   elseif parms.spectrogram_flag % plot spectrogram
     yo = squeeze(tfpows(pop,:,:));
     yy = (abs(yo)+eps).^2;
     yz = (yy-repmat(mean(yy,1),[size(yy,1) 1]))./(repmat(std(yy,0,1),[size(yy,1) 1]));
     imagesc(to,fo,yz'); axis xy; axis tight;
     caxis(zlims);
     ylim(plims); xlim([to(1) to(end)]);
     colorbar; %colormap(flipud(colormap('hot')));
     ylabel('freq [Hz]')   
     title([strrep(lab,'_','\_') ' spectrogram (power z-score)'],'fontsize',14,'fontweight','bold');    
      if pop==npop
        xlabel('time [s]')
      else
        set(gca,'xtick',[],'xticklabel',[]);
      end     
     %surf(to,fo,yo), shading interp, set(gca,'YScale','log','YTick',[[1:10] [15:5:50] [100 200]])
     %view(2), axis tight, ylim([1 250]), colormap(flipud(colormap('hot'))), colorbar
    end
  end
end

% res = PowerSpecTA(t,y,parms.FreqRange,WINDOW,NormAbs,parms.Notch);
% %res = PowerSpecTA(x,y,[10 80],min(8000,round(numel(T)/2)),'Normalized',[]);
% % PowerSpecTA(x,y,FreqRange,Bins,NormAbs,Notch)
% % plot(res.f,log10(res.Pxx));
%     
% res = SpecGramTA(x,y,'mtm',[],8000,7800,2.5);
% %Wind=parms.WINDOW; WindOL=parms.NOVERLAP;
% %   % SpecGramTA(x,y,Smooth,Notch,Wind,WindOL,SmoothWindow)
% %   Fs = fix(1/(x(2)-x(1))); %Fs = roundn(1/(data{1}(2,1)-data{1}(1,1)),3);
% %   d = y.*1e6;
% %   case {'standard','std'}
% %     [yo,fo,to] = specmw(detrend(d),Wind,Fs,Wind,WindOL);
% %   case {'pwelch','pw'}
% %       [yo,fo,to] = mtmspecTA2(detrend(d),Wind,Fs,Wind,WindOL);
% %   case {'multitaper','mtm'}
% %       [yo,fo,to] = mtmspecTA(detrend(d),Wind,Fs,Wind,WindOL,SmoothWindow);
% %   surf(to,fo,yo), shading interp, set(gca,'YScale','log','YTick',[[1:10] [15:5:50] [100 200]])
% %   view(2), axis tight, ylim([1 250]), colormap(flipud(colormap('hot'))), colorbar
% %   %imagesc(to,fo,(abs(yo)+eps).^2); axis xy; colormap(flipud(colormap('hot'))), axis tight, ylim([0 250]), colorbar
% %   xlabel('Time (secs)'), ylabel('Frequency (Hz)')
% 
