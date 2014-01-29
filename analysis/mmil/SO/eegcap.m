% {'EEG 001','EEG 002','EEG 003','EEG 004','EEG 005','EEG 006','EEG 007','EEG 008','EEG 009','EEG 010','EEG 011','EEG 012','EEG 013','EEG 014','EEG 015','EEG 016','EEG 017','EEG 018','EEG 019','EEG 020','EEG 021','EEG 022','EEG 023','EEG 024','EEG 025','EEG 026','EEG 027','EEG 028','EEG 029','EEG 030','EEG 031','EEG 032','EEG 033','EEG 034','EEG 035','EEG 036','EEG 037','EEG 038','EEG 039','EEG 040','EEG 041','EEG 042','EEG 043','EEG 044','EEG 045','EEG 046','EEG 047','EEG 048','EEG 049','EEG 050','EEG 051','EEG 052','EEG 053','EEG 054','EEG 055','EEG 056','EEG 057','EEG 058','EEG 059','EEG 060','EEG 061','EEG 062','EEG 063','EEG 064','EEG 065','EEG 066','EEG 067','EEG 068','EEG 069','EEG 070','EEG 071','EEG 072','EEG 073','EEG 074'}

eeglabels = ...
  {'FP1' 'FPz' 'FP2' 'AF7' 'AF3' 'AFz' 'AF4' 'AF8' 'F7' 'F5' 'F3' 'F1' 'Fz' ...
   'F2' 'F4' 'F6' 'F8' 'FT9' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'FC6' ...
   'FT8' 'FT10' 'T9' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'T10' 'TP9' ...
   'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10' 'P9' 'P7' 'P5' ...
   'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' 'EOG061' 'EOG062' 'ECG063' 'EMG064' ...
   'P10' 'PO7' 'PO3' 'POz' 'PO4' 'PO8' 'O1' 'Oz' 'O2' 'Iz'};
ii = [1:60 65:74]; eeglabels = eeglabels(ii);

xx = [.407 .5 .593 .325 .405 .5 .595 .675 .25 .31 .37 .43 .5 .57 .63 .69 .75 .19 ...
      .26 .325 .38 .44 .5 .56 .62 .675 .74 .81 .18 .245 .31 .38 .44 .5 .56 ...
      .62 .69 .755 .82 .185 .25 .315 .38 .44 .5 .56 .62 .685 .75 .815 .195 ...
      .26 .325 .38 .44 .5 .56 .62 .675 .74 .805 .31 .41 .5 .59 .69 .41 .5 .59 .5];
 
yy = [.89 .89 .89 .85 .81 .78 .81 .85 .78 .76 .74 .72 .71 .72 .74 .76 .78 .67 ...
      .66 .65 .64 .63 .64 .63 .64 .65 .66 .67 .55 .55 .55 .55 .55 .55 .55 ...
      .55 .55 .55 .55 .43 .44 .45 .46 .47 .46 .47 .46 .45 .44 .43 .3 .32 ...
      .35 .37 .39 .39 .39 .37 .35 .32 .3 .25 .29 .31 .29 .25 .21 .23 .21 .15];  
 
figure('Name','Standard Positions');
for k = 1:length(eeglabels)
  subplot('position',[xx(k) yy(k) .04 .04]); axis off
  title(sprintf('%s (%g)',eeglabels{k},ii(k)));
end

%% From Nima
% %subplot('position',[0.5 0.1 0.04 0.04]);  %
% %subplot('position',[0.5 0.17 0.04 0.04]); %
% subplot('position',[0.5 0.15 0.04 0.04]); %74
% subplot('position',[0.5 0.23 0.04 0.04]); %72
% subplot('position',[0.5 0.31 0.04 0.04]); %68
% subplot('position',[0.5 0.39 0.04 0.04]); %56
% subplot('position',[0.5 0.46 0.04 0.04]); %45
% subplot('position',[0.5 0.55 0.04 0.04]); %34
% subplot('position',[0.5 0.64 0.04 0.04]); %23 
% subplot('position',[0.5 0.71 0.04 0.04]); %13
% subplot('position',[0.5 0.78 0.04 0.04]);  %6
% %subplot('position',[0.5 0.87 0.04 0.04]); %
% subplot('position',[0.5 0.89 0.04 0.04]); %2
% 
% %subplot('position',[0.57 0.94 0.04 0.04]);
% subplot('position',[0.56 0.55 0.04 0.04]);%35
% subplot('position',[0.62 0.55 0.04 0.04]);%36
% subplot('position',[0.69 0.55 0.04 0.04]);%37
% subplot('position',[0.755 0.55 0.04 0.04]);%38
% subplot('position',[0.82 0.55 0.04 0.04]);%39
% 
% 
% subplot('position',[0.56 0.47 0.04 0.04]);%46
% %subplot('position',[0.57 0.44 0.04 0.04]);
% subplot('position',[0.56 0.63 0.04 0.04]);%24
% %subplot('position',[0.58 0.76 0.04 0.04]);
% %subplot('position',[0.62 0.84 0.04 0.04]);
% 
% subplot('position',[0.62 0.64 0.04 0.04]);%25
% subplot('position',[0.62 0.46 0.04 0.04]);%47
% subplot('position',[0.675 0.65 0.04 0.04]);%26
% subplot('position',[0.685 0.45 0.04 0.04]);%48
% subplot('position',[0.74 0.66 0.04 0.04]);%27
% subplot('position',[0.75 0.44 0.04 0.04]);%49
% subplot('position',[0.81 0.67 0.04 0.04]);%28
% subplot('position',[0.815 0.43 0.04 0.04]);%50
% 
% subplot('position',[0.56 0.39 0.04 0.04]);%57
% subplot('position',[0.62 0.37 0.04 0.04]);%58
% subplot('position',[0.675 0.35 0.04 0.04]);%59
% subplot('position',[0.74 0.32 0.04 0.04]);%60
% subplot('position',[0.805 0.3 0.04 0.04]);%65
% 
% subplot('position',[0.59 0.29 0.04 0.04]);%69
% subplot('position',[0.59 0.21 0.04 0.04]);%73
% subplot('position',[0.69 0.25 0.04 0.04]);%70
% 
% subplot('position',[0.57 0.72 0.04 0.04]);%14
% subplot('position',[0.63 0.74 0.04 0.04]);%15
% subplot('position',[0.69 0.76 0.04 0.04]);%16
% subplot('position',[0.75 0.78 0.04 0.04]);%17
% 
% subplot('position',[0.595 0.81 0.04 0.04]);%7
% subplot('position',[0.675 0.85 0.04 0.04]);%8
% 
% subplot('position',[0.593 0.89 0.04 0.04]); %3
% 
% 
% %------------------------
% 
% 
% %subplot('position',[0.57 0.94 0.04 0.04]);
% subplot('position',[0.44 0.55 0.04 0.04]);%33
% subplot('position',[0.38 0.55 0.04 0.04]);%32
% subplot('position',[0.31 0.55 0.04 0.04]);%31
% subplot('position',[0.245 0.55 0.04 0.04]);%30
% subplot('position',[0.18 0.55 0.04 0.04]);%29
% 
% 
% subplot('position',[0.44 0.47 0.04 0.04]);%44
% %subplot('position',[0.57 0.44 0.04 0.04]);
% subplot('position',[0.44 0.63 0.04 0.04]);%22
% %subplot('position',[0.58 0.76 0.04 0.04]);
% %subplot('position',[0.62 0.84 0.04 0.04]);
% 
% subplot('position',[0.38 0.64 0.04 0.04]);%21
% subplot('position',[0.38 0.46 0.04 0.04]);%43
% subplot('position',[0.325 0.65 0.04 0.04]);%20
% subplot('position',[0.315 0.45 0.04 0.04]);%42
% subplot('position',[0.26 0.66 0.04 0.04]);%19
% subplot('position',[0.25 0.44 0.04 0.04]);%41
% subplot('position',[0.19 0.67 0.04 0.04]);%18
% subplot('position',[0.185 0.43 0.04 0.04]);%40
% 
% subplot('position',[0.44 0.39 0.04 0.04]);%55
% subplot('position',[0.38 0.37 0.04 0.04]);%54
% subplot('position',[0.325 0.35 0.04 0.04]);%53
% subplot('position',[0.26 0.32 0.04 0.04]);%52
% subplot('position',[0.195 0.3 0.04 0.04]);%51
% 
% subplot('position',[0.41 0.29 0.04 0.04]);%67
% subplot('position',[0.41 0.21 0.04 0.04]);%71
% subplot('position',[0.31 0.25 0.04 0.04]);%66
% 
% subplot('position',[0.43 0.72 0.04 0.04]);%12
% subplot('position',[0.37 0.74 0.04 0.04]);%11
% subplot('position',[0.31 0.76 0.04 0.04]);%10
% subplot('position',[0.25 0.78 0.04 0.04]);%9
% 
% subplot('position',[0.405 0.81 0.04 0.04]);%5
% subplot('position',[0.325 0.85 0.04 0.04]);%4
% 
% subplot('position',[0.407 0.89 0.04 0.04]); %1


