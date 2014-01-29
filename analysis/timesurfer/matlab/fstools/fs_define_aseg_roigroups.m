function roigroups = fs_define_aseg_roigroups()
%function roigroups = fs_define_aseg_roigroups()
%
%
% created:  04/01/09 by Don Hagler
% last mod: 04/01/09 by Don Hagler
%

roigroups(1).roiname = 'WholeBrain';
roigroups(1).roicode = 10001;
roigroups(1).roicodes = [2 3 7 8 10 11 12 13 17 18 26 28 41 42 46 47 ...
                         49 50 51 52 53 54 58 60 78 79];

roigroups(2).roiname = 'LatVentricles';
roigroups(2).roicode = 10002;
roigroups(2).roicodes = [4,5,43,44];

roigroups(2).roiname = 'AllVentricles';
roigroups(2).roicode = 10003;
roigroups(2).roicodes = [4,5,14,15,43,44,72,75,76,213];

