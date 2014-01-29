function plot_ClusterHistograms(clusters,alpha)
if nargin < 2, alpha = .05; end

cnt         = 0; % subplot count
n           = 1;
m           = 5; % type of plots to show for each mod/ref
nbin        = 20;

figure('color','w');
   
sens      = clusters.sensor_info;
CF        = clusters.deg2meters;
clusters  = clusters.epochs;
nchan     = length(sens);
labels    = {sens.label};
R   = [clusters.R];
p   = [clusters.p];
N   = [clusters.N];
str = '';

subplot(n,m,cnt+1),try hist(R,nbin); end; axis tight;
set(gca,'xlim',[-1 1]); ylim=get(gca,'ylim'); xlim=get(gca,'xlim');
xlabel('R'); vline(0,'k');
text(xlim(1)+.1*diff(xlim),ylim(2)-.1*diff(ylim),sprintf('%g SOs',length(R)));
ylabel(str);

subplot(n,m,cnt+2),try hist([clusters.R_RandomLoc],nbin); end; axis tight;
set(gca,'xlim',[-1 1]); ylim=get(gca,'ylim'); xlim=get(gca,'xlim');
xlabel('R (rand loc)'); vline(0,'k');
text(xlim(1)+.1*diff(xlim),ylim(2)-.1*diff(ylim),sprintf('%g SOs',length(R)));

for c = 1:length(clusters)
  rch = strmatch(clusters(c).RefChan,clusters(c).InvolvedChans,'exact');
  tk  = clusters(c).DetectionTimes;
  x   = clusters(c).Distance * CF; % distance in meters
  tau = tk - tk(rch);
  del = abs(tau);
  s(c)= mean( x(tau~=0) ./ del(tau~=0) );
end
s = s(p<alpha);

subplot(n,m,cnt+3),try hist(s,nbin); end; axis tight;
ylim=get(gca,'ylim'); xlim=get(gca,'xlim');
xlabel('speed (m/s)');
str = sprintf('%2.2g+/-%2.2gm/s\n\n%g SOs\np<%g',mean(s),std(s),length(s),alpha);
text(xlim(1)+.5*diff(xlim),ylim(2)-.3*diff(ylim),str);

origins     = arrayfun(@(x)x.RefChan,clusters,'uniformoutput',false);
AllInvChn   = {clusters.InvolvedChans};
AllInvChn   = [AllInvChn{:}];
OriginCount = zeros(nchan,1);
InvolvCount = zeros(nchan,1);
tmp     = {clusters.InvolvedChans}; 
invchan = cellfun(@(x)[x{:}],tmp,'uniformoutput',false);
for k   = 1:nchan
  label = labels{k};
  OriginCount(k) = numel(find(ismember(origins,label)));
  InvolvCount(k) = numel(find(ismember(AllInvChn,label)));
end
ClusterSize = arrayfun(@(x)length(x.InvolvedChans),clusters);
InvolvCount = InvolvCount(InvolvCount~=0);

subplot(n,m,cnt+4),try hist(InvolvCount,nbin); end; axis tight;
ylim=get(gca,'ylim'); xlim=get(gca,'xlim');
xlabel('# SOs');
str = sprintf('%g chans',length(InvolvCount));
text(xlim(1)+.1*diff(xlim),ylim(2)-.2*diff(ylim),str);

subplot(n,m,cnt+5),try hist(ClusterSize,nbin); end; axis tight;
ylim=get(gca,'ylim'); xlim=get(gca,'xlim');
xlabel('# chans'); title('cluster size');
str = sprintf('%g clusters',length(ClusterSize));
text(xlim(1)+.7*diff(xlim),ylim(2)-.2*diff(ylim),str);
    
subplot(1,m,1); % make subplot 1 the current object for setting title later