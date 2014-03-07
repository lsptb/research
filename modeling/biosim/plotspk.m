function [fig,rates,tmins,spiketimes,spikeinds]=plotspk(data,spec,varargin)
% Purpose: plot firing rate(t) and FRH

parms = mmil_args2parms( varargin, ...
                   {  'plot_flag',1,[],... % whether to plot w/ this function
                      'doPlot',0,[],... % whether to plot w/ Grant's function
                      'window_size',1,[],...
                      'dW',1/4,[],...
                      'spikethreshold',0,[],...
                      'smooth',[],[],... % integer span for moving average
                      'var',1,[],...
                      'dt',1/data(1).sfreq,[],...
                   }, false);
                 
if ~isfield(spec,'entities') && isfield(spec,'cells')
  spec.entities=spec.cells;
end
threshold = parms.spikethreshold; % used to create point process {dn}
npop = length(spec.entities);
var = parms.var;
fig = [];

% prepare point process data
spiketimes={}; spikeinds={}; rates={}; alldn={}; spikes={}; orders={};
for pop=1:npop
  t=data(pop).epochs.time;
  ncells=spec.entities(pop).multiplicity; %data(pop).epochs.num_trials;
  if pop==1
    tmins = t(1):parms.dW:t(end)-parms.window_size;
    tmaxs = tmins+parms.window_size;
    %rates{pop}=zeros(ncells,length(tmins)); 
    parms.tmin=t(1);
    parms.tmax=t(end);
  end
  % create point process [dn] = cells x samples
  %   -- increment process (0|1 for point event at each time point)
  dn = zeros(ncells,length(t));
  for cell = 1:ncells
    dat = data(pop).epochs.data(var,:,cell);
    if ~isempty(parms.smooth)
      dat = smooth(dat,parms.smooth);
    end
    [~,ind] = findpeaks(double(dat),'MinPeakHeight',threshold);
    dn(cell,ind) = 1;
    spiketimes{pop}{cell}=t(ind);
    spikeinds{pop}{cell}=ind;
  end
  % calc firing rate from point process
  [frates, jnk, ordered_elec] = frate(dn, parms);
  rates{pop} = frates;
  alldn{pop} = dn;
  spikes{pop} = sum(dn,1);
  orders{pop} = ordered_elec;
end

if parms.plot_flag
  nrows = npop;
  ncols = 3;
  nplots = nrows*ncols;
  %maxtraces=10;
  % ------------------------------
  dx = 1 / ncols;
  dy = 1 / nrows;
  xstep = mod((1:ncols)-1,ncols);
  ystep = mod((1:nrows)-1,nrows)+1;
  xpos = .02+xstep*dx;
  ypos = .02+1-ystep*dy;

  screensize = get(0,'screensize');
  fig=figure('position',screensize);
  set(gca,'Units','normalized','Position',[0 0 1 1]);
  for i = 1:nplots
    xi = mod(i-1,ncols)+1;%xi = ncols-mod(i,ncols);
    yi = floor((i-1)./ncols)+1;
    subplot('Position',[xpos(xi) ypos(yi) .9/ncols .9/nrows]); set(gca,'units','normalized');
    pop = ceil(i/ncols); % index to this population
    lab = data(pop).sensor_info(var).label;   
    r=rates{pop};
    rmin=min(r(:)); rmax=max(r(:)); dr=rmax-rmin;    
    ncells=size(r,1);
    if mod(i,ncols)==1 % plot firing rates
      plot(tmins,r); hold on
      if ncells==1
        % add ticks at spike times
        spikes=t(alldn{pop}==1);
        M=rmax+.1*dr;
        for k=1:length(spikes)
          text(double(spikes(k)),M,'|','color',plot_color);
        end
      else
        % overlay mean
        plot(tmins,mean(r,1), 'LineWidth', 3);
      end
      lims = [max(0,rmin-.2*dr) rmax+.2*dr];
      %lims = [max(0,min(r(:)) 1.3*max(r(:))]
      ylim(lims);
      xlim([tmins(1) tmins(end)]); 
      %xlim([t(1) t(end)]);
      xlabel('time [s]');
      ylabel('firing rate [Hz]');
      %text(min(xlim)+.2*diff(xlim),min(ylim)+.8*diff(ylim),[strrep(lab,'_','\_') ' spike rate (threshold=' num2str(threshold) ')'],'fontsize',14,'fontweight','bold');
      title([strrep(lab,'_','\_') ' spike rate (threshold=' num2str(threshold) ')'],'fontsize',14,'fontweight','bold');
    elseif mod(i,ncols)==2 % plot heat map
      imagesc(tmins, 1:ncells,r(orders{pop}, :)); axis xy
      %text(min(xlim)+.2*diff(xlim),min(ylim)+.8*diff(ylim),[strrep(lab,'_','\_') ' sorted heat map'],'fontsize',14,'fontweight','bold');
      title([strrep(lab,'_','\_') ' sorted heat map'],'fontsize',14,'fontweight','bold');
      ylabel('cell','fontsize',12); caxis(lims);
      colorbar
    else % plot psth      
      psth = spikes{pop}/ncells;
      plot(t,psth,'LineWidth',3);
      ylim([0, 1]); xlim([t(1) t(end)]);
      %text(min(xlim)+.2*diff(xlim),min(ylim)+.8*diff(ylim),[strrep(lab,'_','\_') ' PSTH'],'fontsize',14,'fontweight','bold');
      title([strrep(lab,'_','\_') ' PSTH'],'fontsize',14,'fontweight','bold');
    end
  end
end

