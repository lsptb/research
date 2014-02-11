function plotv(data,spec)
addpath(genpath('/project/crc-nak/sherfey/code/research/external/tallie'));

keyboard

npop = length(spec.entities);

% Minimize white space
set(gca,'Units','normalized','Position',[0 0 1 1]);

% Fill screen with subplots:
nrows = npop;%ceil(sqrt(nplots));
ncols = 2;%nrows;
nplots = nrows*ncols;
 
% ------------------------------
dx = 1 / ncols;
dy = 1 / nrows;
 
xstep = mod((1:ncols)-1,ncols);
ystep = mod((1:nrows)-1,nrows)+1;
xpos = .02+xstep*dx;
ypos = .02+1-ystep*dy;
 
screensize = get(0,'screensize');
figure('position',screensize);
for i = 1:nplots
  xi = mod(i-1,ncols)+1;%xi = ncols-mod(i,ncols);
  yi = floor((i-1)./ncols)+1;
  subplot('Position',[xpos(xi) ypos(yi) 1/ncols 1/nrows]); set(gca,'units','normalized');
  pop = ceil(i/ncols); % index to this population
  var = 1;
  if mod(i,2)==1 % odd, plot trace
    T = data(pop).epochs.time;
    n = data(pop).epochs.num_trials;
    dat = squeeze(data(pop).epochs.data(var,:,:))';
    lab = data(pop).sensor_info(var).label;
    imagesc(T,1:n,dat); axis xy;
    % overlay LFP
    lfp = sum(dat,1)';
    % ...
    text(min(xlim)+.2*diff(xlim),min(ylim)+.8*diff(ylim),strrep(lab,'_','\_'),'fontsize',14);
  else % even, plot spectrum
    fh = @(x,y) PowerSpecTA(x,y,[10 80],min(8000,round(numel(T)/2)),'Normalized',[]);
    res = feval(fh,T,lfp);
    set(gca,'xdata',res.f,'ydata',log10(res.Pxx));
  end
end

