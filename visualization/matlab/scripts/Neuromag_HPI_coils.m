fif='/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_19_raw.fif';
megmodel('head',fif,'head'); 
[n,s,t]=chaninfo;
% n = Chan ID
% s = Chan Name
% t = Coil coordinate transformations
[ty,na] = chaninfo('type');
% ty = coil type (0=mag,1=planar grad)
% na = name

{s(1,:) data.sensor_info(1).label}

devicecoord = data.sensor_info(1).loc
device2head = data.coor_trans.device2head
headcoord   = devicecoord * device2head

coiltransformation = t{1};
ts_extract_trans_from_fiff
ts_extract_hpts_from_fiff

[coords,kind,num] = hpipoints(datafile); % kind: 1=cardinal,2=HPI,3=EEG,4=EXTRA

coords(:,(kind==2))

ref_EEG_coords = ts_extract_hpts_from_fiff(fif);


