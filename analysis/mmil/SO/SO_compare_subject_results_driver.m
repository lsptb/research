subjects        = [1 2 4 5 6 8];
outfile         = [];%'/space/emc2/1/halgdev/projects/sleep/MEG/SO/images/noflip/alltrials/workspace_AllSubjects_Summary_EarliestDetection_noflip_alltrials.mat';
outpath         = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/images';

FileParams(1,:) = {'SO_clustered_detections.mat',           'SO_clusters_noflip.mat',     'clusters_Near_all',1,'noflip-alltrials-earliest'};
FileParams(2,:) = {'SO_clustered_detections.mat',           'SO_clusters_noflip.mat',     'clusters_Near_pks',1,'noflip-postrials-earliest'};
FileParams(3,:) = {'SO_clustered_detections.mat',           'SO_clusters_noflip.mat',     'clusters_Near_pks',2,'noflip-negtrials-earliest'};
FileParams(4,:) = {'SO_clustered_consistent_detections.mat','SO_clusters_consistency.mat','clusters_Cear_all',1,'consistency-alltrials-earliest'};
FileParams(5,:) = {'SO_clustered_consistent_detections.mat','SO_clusters_consistency.mat','clusters_Cear_pks',1,'consistency-postrials-earliest'};
FileParams(6,:) = {'SO_clustered_consistent_detections.mat','SO_clusters_consistency.mat','clusters_Cear_pks',2,'consistency-negtrials-earliest'};
FileParams(7,:) = {'SO_clustered_flipmatrix_detections.mat','SO_clusters_flipmatrix.mat', 'clusters_Mear_all',1,'flipmatrix-alltrials-earliest'};
FileParams(8,:) = {'SO_clustered_flipmatrix_detections.mat','SO_clusters_flipmatrix.mat', 'clusters_Mear_pks',1,'flipmatrix-postrials-earliest'};
FileParams(9,:) = {'SO_clustered_flipmatrix_detections.mat','SO_clusters_flipmatrix.mat', 'clusters_Mear_pks',2,'flipmatrix-negtrials-earliest'};
FileParams(10,:) = {'SO_clustered_detections.mat',           'SO_clusters_noflip.mat',     'clusters_Nmax_all',1,'noflip-alltrials-maxR'};
FileParams(11,:) = {'SO_clustered_detections.mat',           'SO_clusters_noflip.mat',     'clusters_Nmax_pks',1,'noflip-postrials-maxR'};
FileParams(12,:) = {'SO_clustered_detections.mat',           'SO_clusters_noflip.mat',     'clusters_Nmax_pks',2,'noflip-negtrials-maxR'};
FileParams(13,:) = {'SO_clustered_consistent_detections.mat','SO_clusters_consistency.mat','clusters_Cmax_all',1,'consistency-alltrials-maxR'};
FileParams(14,:) = {'SO_clustered_consistent_detections.mat','SO_clusters_consistency.mat','clusters_Cmax_pks',1,'consistency-postrials-maxR'};
FileParams(15,:) = {'SO_clustered_consistent_detections.mat','SO_clusters_consistency.mat','clusters_Cmax_pks',2,'consistency-negtrials-maxR'};
FileParams(16,:) = {'SO_clustered_flipmatrix_detections.mat','SO_clusters_flipmatrix.mat', 'clusters_Mmax_all',1,'flipmatrix-alltrials-maxR'};
FileParams(17,:) = {'SO_clustered_flipmatrix_detections.mat','SO_clusters_flipmatrix.mat', 'clusters_Mmax_pks',1,'flipmatrix-postrials-maxR'};
FileParams(18,:) = {'SO_clustered_flipmatrix_detections.mat','SO_clusters_flipmatrix.mat', 'clusters_Mmax_pks',2,'flipmatrix-negtrials-maxR'};

selset = 1;%4:6;

clear AllSets
for s = 1:length(selset)
  f   = selset(s);
  DetectionFile   = FileParams{f,1};
  ClusterFile     = FileParams{f,2};
  StructSelection = FileParams{f,3};
  ArrayIndex      = FileParams{f,4};
  SetID           = FileParams{f,5};
  [res,grand]     = SO_compare_subject_results(subjects,ClusterFile,DetectionFile,StructSelection,ArrayIndex,outfile,SetID);
  AllSets(s).SetID    = SetID;
  AllSets(s).params   = [FileParams{f,:} {subjects}];
  AllSets(s).subjects = res;
  AllSets(s).grandavg = grand;
%   figfile = sprintf('%s/AllSubjects_Summary_%s.jpg',outpath,SetID);
%   print('-djpeg',figfile);
%   figfile = sprintf('%s/AllSubjects_Summary_%s.eps',outpath,SetID);
%   print('-depsc','-tiff','-r300',figfile);
end
% for k=1:7,subplot(6,7,k); set(gca,'clim',[0 .5]); end;set(gca,'clim',[0 1])
% for k=8:14,subplot(6,7,k); set(gca,'clim',[0 .03]); end
% for k=15:21,subplot(6,7,k); set(gca,'clim',[.08 .10]); end
% % subplot(6,7,22); set(gca,'ylim',[0 4000]); subplot(6,7,29); set(gca,'ylim',[0 300]);
% % subplot(6,7,23); set(gca,'ylim',[0 3000]); subplot(6,7,30); set(gca,'ylim',[0 250]);
% % subplot(6,7,24); set(gca,'ylim',[0 6000]); subplot(6,7,31); set(gca,'ylim',[0 500]);
% % subplot(6,7,25); set(gca,'ylim',[0 1750]); subplot(6,7,32); set(gca,'ylim',[0 150]);
% % subplot(6,7,26); set(gca,'ylim',[0 1750]); subplot(6,7,33); set(gca,'ylim',[0 150]);
% % subplot(6,7,27); set(gca,'ylim',[0 2500]); subplot(6,7,34); set(gca,'ylim',[0 200]);
% % subplot(6,7,28); set(gca,'ylim',[0 20000]); subplot(6,7,35); set(gca,'ylim',[0 1500]);
% for k=36:41,subplot(6,7,k); set(gca,'ylim',[-.03 .03]); end
% subplot(6,7,42),set(gca,'ylim',[-.003 .003]);

  set(gcf,'PaperOrientation' ,'portrait',...
          'PaperType'        ,'usletter',...
          'PaperUnits'       ,'inches'  ,...
          'PaperPositionMode','manual'  ,...
          'InvertHardcopy'   ,'off'     );
  ppsize = get(gcf,'PaperSize');
  width  = ppsize(2) - 2;
  height = ppsize(1) - 2;
  left   = (ppsize(1) - width ) / 2;
  bottom = (ppsize(2) - height) / 2;
  set(gcf,'PaperPosition',[left bottom width height]);
  pos    = get(gcf,'Position');
  if (pos(3)/pos(4)) > (11/8.5)
    pos(3) = (11/8.5)*pos(4);
    set(gcf,'position',pos);
  end
  file = sprintf('%s/AllSubjects_Summary_%s_Fig4_percents.eps',outpath,SetID);
  if exist(file,'file')
    fprintf('not saving, figure already exists:%s\n',file);
  else
    print(gcf,'-depsc2','-painters','-loose',file);
  end
end

% save results
matfile = sprintf('%s/AllSubjects_Summary_sets%s.mat',outpath,strrep(num2str(selset),'  ','_'));
save(matfile,'subjects','FileParams','AllSets');



% %% Format the figures for CorelDraw and Illustrator
% outpath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/images/OriginDetectionDensity/eps';
% for k   = 1:length(selset)
%   file  = sprintf('%s/AllSubjects_Summary_%s.eps',outpath,FileParams{selset(k),5});
%   figure(k);
%   % Format subplots
%   for s = 1:21
%     subplot(7,7,s);
%     h = get(gca,'Children');
%     delete(h(1));
%     set(h(2:5),'LineWidth',.2);
%     set(h(6),'Marker','.','MarkerSize',.1);
%   end
%   for s = 22:49
%     subplot(7,7,s);
%     set(gca,'FontSize',5);
%     set(get(gca,'xlabel'),'fontsize',5);
%     set(get(gca,'ylabel'),'fontsize',5);
%     if s >= 22 && s <= 28
%       h = get(gca,'Children');
%       str = get(h(1),'String'); set(h(1),'String',[str(1:4) str(end-2:end)]);
%       ylim=get(gca,'ylim'); set(h(1),'Position',[7 .8*ylim(2)],'FontSize',6);      
%     end
%     if s >= 29 && s <= 35
%       h = get(gca,'Children');
%       str = get(h(1),'String'); 
%       set(h(1),'String',[str(1:5) str(end-3:end)]);
%       ylim=get(gca,'ylim'); set(h(1),'Position',[-.9 .8*ylim(2)],'FontSize',6);
%       delete(h(2));
%       set(h(3),'MarkerSize',5);
%       vline(0,'k');
%     end
%     if s >= 36 && s <= 42
%       h = legend;
%       str = get(h,'String');
%       ix = findstr(',',str{1});
%       str{1} = [str{1}(1:16) str{1}(ix:ix+7)];
%       str{2} = [str{2}(1:16) str{2}(ix:ix+7)];
%       set(h,'String',str,'FontSize',4,'Location','NorthOutside');
%       set(gca,'xlim',[0 .1]);
%       h = get(gca,'Children');
%       set(h(1:2),'MarkerSize',2);
%     end
%     if s >= 43
%       set(gca,'xlim',[0 12]);
%       h = get(gca,'Children');
%       set(h(2),'MarkerSize',2);
%       if s < 49
%         set(h(3),'MarkerSize',2);
%       end
%     end
%   end
%   set(gcf,'PaperOrientation' ,'portrait',...
%           'PaperType'        ,'usletter',...
%           'PaperUnits'       ,'inches'  ,...
%           'PaperPositionMode','manual'  ,...
%           'InvertHardcopy'   ,'off'     );
%   ppsize = get(gcf,'PaperSize');
%   width  = ppsize(2) - 2;
%   height = ppsize(1) - 2;
%   left   = (ppsize(1) - width ) / 2;
%   bottom = (ppsize(2) - height) / 2;
%   set(gcf,'PaperPosition',[left bottom width height]);
%   pos    = get(gcf,'Position');
%   if (pos(3)/pos(4)) > (11/8.5)
%     pos(3) = (11/8.5)*pos(4);
%     set(gcf,'position',pos);
%   end
%   if exist(file,'file')
%     fprintf('not saving, figure already exists:%s\n',file);
%   else
%     print(gcf,'-depsc2','-painters','-loose',file);
%   end
% end

% subplots 1-21, Children:
% 1 = grad-label  => delete
% 2 = left ear    => set(_,'LineWidth',.2);
% 3 = right ear   => set(_,'LineWidth',.2);
% 4 = nose        => set(_,'LineWidth',.2);
% 5 = head        => set(_,'LineWidth',.2);
% 6 = elec (o)    => set(_,'Marker','.','MarkerSize',.1)

% subplots 22-49
% set(gca,'FontSize',5)

% subplots 29-35, Children:
% 1 = text => str = get(h(1),'String'); set(h(1),'String',[str(1:5) str(end-3:end)]);
% 2 = vline => delete
% 3 = data => set(h(3),'MarkerSize',5);
% vline(0,'k');

% subplots 36-42
% set(gca,'xlim',[0 .1]);

% subplots 43-48, Children:
% 2 = averages => set(_,'MarkerSize',2);
% 3 = points => set(_,'MarkerSize',2)

% subplot 49
% 2 = averages => set(_,'MarkerSize',2);
