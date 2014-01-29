% Find example for Figure 3 (18-Aug-2010, JSS)
cd ~jsherfey/projects/sleep3/
load s1/matfiles/proc_epoch_data_ICA.mat
load s1/matfiles/SO_clusters_consistency.mat
R=[clusters_Cmax_all(1).epochs.R];
p=[clusters_Cmax_all(1).epochs.p];
N=[clusters_Cmax_all(1).epochs.N];
T=[clusters_Cmax_all(1).epochs.RefTime];
s=[clusters_Cmax_all(1).epochs.Speed];

Rsig  = p < .05;
Rhi   = R > .7;
Rlo   = R < .3;
Nth   = 40;%N > 40;

sel   = find((Rhi(2:end)==1 & Rsig(2:end)) & (Rlo(1:end-1)==1 & Rsig(1:end-1)));
% sel   = find((Rhi(1:end-1)==1 & Rsig(1:end-1)) & (Rlo(2:end)==1 & Rsig(2:end)));
Tlo   = T(sel);
Rvlo  = R(sel);
Nvlo  = N(sel);
Thi   = T(sel+1);
Rvhi  = R(sel+1);
Nvhi  = N(sel+1);

% [Rvlo' Nvlo' Tlo' Rvhi' Nvhi' Thi']

sel   = find((Nvlo >= Nth) | (Nvhi >= Nth));
Tlo   = Tlo(sel);
Rvlo  = Rvlo(sel);
Nvlo  = Nvlo(sel);
Thi   = Thi(sel);
Rvhi  = Rvhi(sel);
Nvhi  = Nvhi(sel);

[Rvlo' Nvlo' Tlo' Rvhi' Nvhi' Thi']



% sync example
trl = 6498;
cst = clusters_Cmax_all(1).epochs(trl);
tmp = ts_matrix2avg([cst.Delays;cst.Delays]');
[s1,s2] = match_str({clusters_Cmax_all(1).sensor_info.label},cst.InvolvedChans);
tmp.sensor_info   = clusters_Cmax_all(1).sensor_info(s1);
tmp.averages.data = tmp.averages.data(:,1);
tmp.averages.time = 0;

tmpdat=ts_data_selection(tmp,'chantype','grad1'); highlight=24;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight'); colorbar
tmpdat=ts_data_selection(tmp,'chantype','grad2'); highlight=18;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight'); colorbar

% traveling wave example
trl = 3096;
cst = clusters_Cmax_all(1).epochs(trl);
tmp = ts_matrix2avg([cst.Delays;cst.Delays]');
[s1,s2] = match_str({clusters_Cmax_all(1).sensor_info.label},cst.InvolvedChans);
tmp.sensor_info   = clusters_Cmax_all(1).sensor_info(s1);
tmp.averages.data = tmp.averages.data(:,1);
tmp.averages.time = 0;

tmpdat=ts_data_selection(tmp,'chantype','grad1'); highlight=11;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight'); colorbar
tmpdat=ts_data_selection(tmp,'chantype','grad2'); highlight=5;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight'); colorbar

%%
sens = clusters_Cmax_all(1).sensor_info;

% sync example
trl = 6498;
cst = clusters_Cmax_all(1).epochs(trl);
tmp = ts_matrix2avg([cst.Delays;cst.Delays]');
[s1,s2] = match_str({clusters_Cmax_all(1).sensor_info.label},cst.InvolvedChans);
tmp.sensor_info   = clusters_Cmax_all(1).sensor_info(s1);
tmp.averages.data = tmp.averages.data(:,1);
tmp.averages.time = 0;

tmpdat=ts_data_selection(tmp,'chantype','grad1'); highlight=24;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight','newfig',1); colorbar
tmpdat=ts_data_selection(tmp,'chantype','grad2'); highlight=18;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight','newfig',0); colorbar

figure
plot(cst.Delays,cst.Theta3D,'.'); axis tight; lsline;
figure; highlight = s1; dewar

% sync example
trl = 3096;
cst = clusters_Cmax_all(1).epochs(trl);
tmp = ts_matrix2avg([cst.Delays;cst.Delays]');
[s1,s2] = match_str({clusters_Cmax_all(1).sensor_info.label},cst.InvolvedChans);
tmp.sensor_info   = clusters_Cmax_all(1).sensor_info(s1);
tmp.averages.data = tmp.averages.data(:,1);
tmp.averages.time = 0;

tmpdat=ts_data_selection(tmp,'chantype','grad1'); highlight=11;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight','newfig',1); colorbar
tmpdat=ts_data_selection(tmp,'chantype','grad2'); highlight=5;
ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'electrodes','on','style','straight','newfig',0); colorbar

figure
plot(cst.Delays,cst.Theta3D,'.'); axis tight; lsline;

figure; highlight = s1; dewar

