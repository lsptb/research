function SpikeFieldCoh_RoseAll_TA(Ang_peak,DataE_st,DataE_en,varargin)


cmap1 = bone(length(Ang_peak)+3);
if length(varargin)==2
    figure, set(gcf,'Position',[625 45 675 940])
    subplot(3,1,3), [~,h1,h2] = plotyy(1:length(varargin{1}),varargin{1},1:length(varargin{2}),varargin{2});
    set(h1,'LineStyle','-','Marker','o'), set(h2,'LineStyle',':','Marker','o'), legend('Spike rate','Number of spikes'), xlabel('Depol level')
    hs = subplot(3,1,[1 2]);
else hs = figure;
end

for k = 1:length(Ang_peak)
    hp(k) = polar(Ang_peak{k},10-k,'-o');
    set(hp(k),'LineWidth',5,'Color',cmap1(k,:))
    hold on
    sdi{k} = linspace(DataE_st{k},DataE_en{k},100);
    hp2(k) = polar(sdi{k},(10-k)*ones(1,100),'-k');
    set(hp2(k),'LineWidth',2,'Color',cmap1(k,:))
end
pt = findall(hs,'type','text');
pl = findall(hs,'type','line');
for k = 1:length(Ang_peak)
    pl(pl==hp(k)) = []; pl(pl==hp2(k)) = [];
end
delete(pt), delete(pl)
pp = polar(0,0,'+k'); set(pp,'MarkerSize',10)


end

