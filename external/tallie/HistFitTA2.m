function [X,Y,PD,Hist,Ang_peak,Ang_std,Distrib_kurt,Peak_Density,DataE_st,DataE_en] = HistFitTA2(bins,BandWidth,x1,varargin)


if nargin>3, PlotsYN = varargin{1}; end 

% shift the distribution centrally to make a sensibly shaped fit later
D = x1(:);
Xhist1 = linspace(-pi,pi,72);
Dhist1 = histc(D,Xhist1); % 5 degree resolution on angle shift
[~,Dhist1_maxloc] = max(Dhist1);
ShiftAng = Xhist1(Dhist1_maxloc);
D2 = D - ShiftAng;
D2(D2<(-pi)) = D2(D2<(-pi)) + (2*pi);
D2(D2>pi) = D2(D2>pi) - (2*pi);

% create and plot the density histogram
[CdfF,CdfX] = ecdf(D2,'Function','cdf');
BinInfo.rule = 3; BinInfo.nbins = bins;
[~,BinEdge] = internal.stats.histbins(D,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);

if strcmp('Y',PlotsYN)
    figure, maximize(gcf), subplot(2,3,1), hold on
    hb = bar(BinCenter,BinHeight,'hist');
    xlabel('Data'), ylabel('Density'), set(hb,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9])
end

% overlay the probability density of the non-parametric fit
PD = fitdist(D2,'kernel','kernel','epanechnikov','support','unbounded','width',BandWidth);
XGrid = linspace(-2*pi,2*pi,1000);
YPlot = pdf(PD,XGrid);
if strcmp('Y',PlotsYN), plot(XGrid,YPlot,'-k','LineWidth',2), axis tight, end

% rose the original data
if strcmp('Y',PlotsYN), subplot(2,3,2), hr = rose(x1,bins); patch(get(hr,'XData'),get(hr,'YData'),[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]); end

% polar plot the fit
XGrid2 = XGrid + ShiftAng;
XGst = findnearest(-pi+ShiftAng,XGrid2);
XGen = findnearest(pi+ShiftAng,XGrid2);
if strcmp('Y',PlotsYN), subplot(2,3,4), hp = polar(XGrid2(XGst:XGen),YPlot(XGst:XGen)); set(hp,'LineWidth',3,'Color',[0 0.85 0.85]), end

% polar plot the peak and standard error
[~,FitPeak_loc] = max(YPlot);
Ang_peak = XGrid2(FitPeak_loc);
Ang_std = std(D2);
Distrib_kurt = kurtosis(YPlot(XGst:XGen));
Peak_Density = max(YPlot);
s1 = std(D2)/2; %semTA(D2);
DataE_st = Ang_peak - s1;
DataE_en = Ang_peak + s1;

if strcmp('Y',PlotsYN)
    subplot(2,3,5)
    hp = polar(Ang_peak,9,'-o'); set(hp,'LineWidth',5,'Color',[0.25 0.25 0.25])
    hold on
    sdi = linspace(DataE_st,DataE_en,100);
    hp2 = polar(sdi,9*ones(1,100));
    set(hp2,'LineWidth',2,'Color',[0.25 0.25 0.25])
    pt = findall(gca,'type','text');
    pl = findall(gca,'type','line');
    pl(pl==hp) = []; pl(pl==hp2) = [];
    delete(pt), delete(pl)
    pp = polar(0,0,'+'); set(pp,'MarkerSize',10)
end

X = XGrid2;
Y = YPlot;
Hist = BinHeight;


end

