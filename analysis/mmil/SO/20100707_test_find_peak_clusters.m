load /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s2/s2_SO_init_peaks_filt0.01-4Hz_toi0-8855.14_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_07-Jul-2010.mat % peaks

method            = 'histogram';
thresh            = 'meanstd';    % threshold for defining histogram peaks
StepSize          = .01;          % step size in sliding aggregate count
IntegrationWindow = .1;           % size of sliding window in aggregate count
MinSeparation     = .1;           % combine peaks closer than this
ClusterWindow     = 1;            % cluster window size
peaktype          = 'both';       % peaktype to use for aggregate count
cluster_peaks = find_peak_clusters(peaks,'method',method,'thresh',thresh,...
  'StepSize',StepSize,'IntegrationWindow',IntegrationWindow,'MinSeparation',MinSeparation,...
  'count',count,'ClusterWindow',ClusterWindow);

peaks    = cluster_peaks;
events   = [];
for k   = 1:length(peaks)
  events(k).label = data.sensor_info(k).label;
  events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
  events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
end
save('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s2/s2_cluster_peaks_20100707.mat','events','peaks');

% load s2 data
visualizer(data); % load s2_cluster_peaks

