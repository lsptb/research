function detections = Constrain_ClusterSize(detections,MinChansPerCluster)
% Purpose: constrain cluster size
% Created by Jason Sherfey on 08-Nov-2010

t         = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
typestr   = '';
cnumfield = sprintf('%scluster_number',typestr);
cindfield = sprintf('%scluster_time_index',typestr);
clusterID = arrayfun(@(x)(x.(cnumfield)),detections,'UniformOutput',false);
tmpID     = [clusterID{:}];
clusterID = unique(tmpID);
clusterIX = arrayfun(@(x)(x.(cindfield)),detections,'UniformOutput',false);
clusterIX = [clusterIX{:}];
clusterIX = cellfun(@(x)unique(clusterIX(tmpID==x)),num2cell(unique(tmpID)));
nclusters = length(clusterID);
cnum      = {detections.(cnumfield)};
cnum      = [cnum{:}];
Ninvolved = cellfun(@(x)sum(cnum==x),num2cell(clusterID));
Nthresh   = MinChansPerCluster;
ckeep     = find(Ninvolved > Nthresh);
Ninvolved = Ninvolved(ckeep);
clusterID = clusterID(ckeep);
clusterIX = clusterIX(ckeep);
for k = 1:length(detections)
  pkeep   = find(ismember(detections(k).pospeak_cluster_number,clusterID));
  nkeep   = find(ismember(detections(k).negpeak_cluster_number,clusterID));
  bkeep   = find(ismember(detections(k).cluster_number,clusterID));
  detections(k).pospeak                     = detections(k).pospeak(pkeep);
  detections(k).pospeak_cluster_number      = detections(k).pospeak_cluster_number(pkeep);
  detections(k).pospeak_cluster_time_index  = detections(k).pospeak_cluster_time_index(pkeep);
  detections(k).negpeak                     = detections(k).negpeak(nkeep);
  detections(k).negpeak_cluster_number      = detections(k).negpeak_cluster_number(nkeep);
  detections(k).negpeak_cluster_time_index  = detections(k).negpeak_cluster_time_index(nkeep);
  detections(k).cluster_number              = detections(k).cluster_number(bkeep);
  detections(k).cluster_time_index          = detections(k).cluster_time_index(bkeep);
end
