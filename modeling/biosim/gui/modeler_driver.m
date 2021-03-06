%% Cell spiking in modeler()
cell = []; net=[];
cell.label = 'E';        % --> cell name
cell.multiplicity = 2;   % --> number of cells
cell.mechanisms = {'itonic' 'ileak','iK','iNa','noise'};
cell.parameters = {'Cm',1};
cell.dynamics = 'V''=(current)/Cm';
net.cells=cell;
net.connections.label = 'E-E';
net.connections.mechanisms='iSYN';
net.connections.parameters=[];
modeler(net);

%% Cell Modeler (1-compartment): driver to define initial cell model and launch the cell modeler
%addpath /space/mdeh3/9/halgdev/projects/jsherfey/code/modeler;
addpath('C:\Users\jsherfey\Desktop\My World\Code\modelers');
addpath('C:\Users\jsherfey\Desktop\My World\Code\modelers\biosim');

cell = [];
cell.label = 'E';        % --> cell name
cell.multiplicity = 1;   % --> number of cells
cell.mechanisms = {'itonic','ileak','iK','iNa','noise'}; % predefined: get_mechlist
cell.parameters = {'stim',10,'E_l',-54.4,'g_l',.3,'Cm',1,'gNa',120,'EKf',-77,'V_noise',3,'IC_noise',1};
cell.dynamics = 'V''=(current)';

cellmodeler(cell); % rt_biosim(cell); data=biosim(cell); biosim_plots(data);

%% cell characterization
clear all
cell = [];
cell.label = 'E';        % --> cell name
cell.multiplicity = 5;   % --> number of cells
cell.mechanisms = {'iStepProtocol' 'ileak','iK','iNa','noise'}; % predefined: get_mechlist
cell.parameters = {'E_l',-54.4,'g_l',.3,'Cm',1,'ENa',50,'gNa',120,'EKf',-77,'V_noise',2,...
  'isi',1000,'nsteps',3,'steptime',400,'nsections',3,'membranearea',2500,'tonictime',6000,'dt',.01};
cell.dynamics = 'V''=(current)'; % note: set timelimits=[0 25000]
%cellmodeler(cell);
% 
% [X,t,Iinj]=SimCharacterizeCells(cell,120*30);
% figure; 
% subplot(2,1,1); plot(t,X(:,1)); xlim([min(t) max(t)]); 
% subplot(2,1,2); plot(t,Iinj); xlim([min(t) max(t)]);

cells.cells=cell;
cells.connections.mechanisms=[];
cells.connections.parameters=[];
cells.connections.label = [];
netmodeler(cells);

% E_iStepProtocol_I    = getStepProtocolStim((0.02),(500),(2),(200),(100),(1500),(3),(60000));
% figure; V=X(:,1);
% subplot(2,1,1); plot(V); xl=xlim; 
% subplot(2,1,2); plot(E_iStepProtocol_I); xlim(xl);

%% Cell spiking
clear all
cell = [];
cell.label = 'E';        % --> cell name
cell.multiplicity = 5;   % --> number of cells
cell.mechanisms = {'itonic' 'ileak','iK','iNa','noise'}; % predefined: get_mechlist
cell.parameters = {'E_l',-54.4,'g_l',.3,'Cm',1,'ENa',50,'gNa',120,'EKf',-77,'V_noise',2,...
  'stim',25,'dt',.01};
cell.dynamics = 'V''=(current)'; % note: set timelimits=[0 25000]
cells.cells=cell;
cells.connections.mechanisms=[];
cells.connections.parameters=[];
cells.connections.label = [];
netmodeler(cells);

%% Cell Modeler (2-compartments): driver to define initial cell model and launch the cell modeler
cd('C:\Users\jsherfey\Desktop\My World\Code\research\modeling\biosim\gui');
addpath(genpath('C:\Users\jsherfey\Desktop\My World\Code\modelers'));
close all; clc;

cell=[]; cells=[];
cell.label = 'soma';        % --> cell name
cell.multiplicity = 1;   % --> number of cells
cell.mechanisms = {'itonic','ileak','iK','iNa','noise'}; % predefined: get_mechlist
cell.parameters = {'stim',10,'E_l',-54.4,'g_l',.3,'Cm',1,'gNa',120,'EKf',-77,'V_noise',3,'IC_noise',1};
cell.dynamics = 'V''=(current)';
cells.cells(1) = cell;

cell=[];
cell.label = 'dend';        % --> cell name
cell.multiplicity = 1;   % --> number of cells
cell.mechanisms = {'itonic','ileak','iK','iNa','noise'}; % predefined: get_mechlist
cell.parameters = {'stim',10,'E_l',-54.4,'g_l',.3,'Cm',1,'gNa',120,'EKf',-77,'V_noise',3,'IC_noise',1};
cell.dynamics = 'V''=(current)';
cells.cells(2) = cell;

% Connect cells
from=2; to=1; % Ed->Es
cells.connections(from,to).label = [cells.cells(from).label '-' cells.cells(to).label];
cells.connections(from,to).mechanisms = 'iCOM';
cells.connections(from,to).parameters = {'g_COM',.2};

cellmodeler(cells); % rt_biosim(cells); data=biosim(cell); biosim_plots(data);

% - update global entity param in model: change param in spec.cells(1).parameters = {key,val} and pass through buildmodel2()
% - update mech param in model: change param in spec.cells(1).mechs(1).params.(key) and pass through buildmodel2()
% - update mech list in model: 
% 	add: add mech label to spec.cells(1).mechanisms{end+1} and file to spec.files{end+1}, then pass through buildmodel2()
% 	remove: remove mech label from spec.cells(1).mechanisms and model from spec.cells(1).mechs
% - update compartments/cells in model:
% 	add: add spec.cells(end+1) and give a unique spec.cells(#).label; pass through buildmodel2()
% 	remove: spec.cells(#)=[], spec.connections(:,#)=[], spec.connections(#,:)=[]

%% Network Modeler: driver to define initial network model and launch the network modeler
%addpath /space/mdeh3/9/halgdev/projects/jsherfey/code/modelers;
%cd('C:\Users\jsherfey\Desktop\My World\Code');
%addpath(genpath('C:\Users\jsherfey\Desktop\My World\Code\research\modeling'));
%addpath('C:\Users\jsherfey\Desktop\My World\Code\modelers\biosim');
%addpath(genpath('/project/crc-nak/sherfey/code/research/external/tallie'));

clear global; clear all; close all;
net = [];

% Define E-cell population
cell = [];
cell.label = 'E';         % --> cell name
cell.multiplicity = 40;   % --> number of cells
cell.mechanisms = {'itonic','ileak','iK','iNa','noise'};
cell.parameters = {'stim',10,'E_l',-54.4,'g_l',.3,'Cm',1,'gNa',120,'EKf',-77,'V_noise',3,'IC_noise',1};
cell.dynamics = 'V''=(current)';
net.cells(1) = cell;

% Define I-cell population
cell = [];
cell.label = 'I';         % --> cell name
cell.multiplicity = 10;   % --> number of cells
cell.mechanisms = {'itonic','ileak','iK','iNa','noise'};
cell.parameters = {'stim',-10,'E_l',-54.4,'g_l',.3,'Cm',1,'gNa',120,'EKf',-77,'V_noise',3,'IC_noise',1};
cell.dynamics = 'V''=(current)';
net.cells(2) = cell;

% Connect cells
from=1; to=2; % E->I
net.connections(from,to).label = [net.cells(from).label '-' net.cells(to).label];
net.connections(from,to).mechanisms = 'iSYN';
net.connections(from,to).parameters = {'g_SYN',.1,'tauRx',.4,'tauDx',2,'E_SYN',0,'IC_noise',1};

from=2; to=1; % I->E
net.connections(from,to).label = [net.cells(from).label '-' net.cells(to).label];
net.connections(from,to).mechanisms = 'iSYN';
net.connections(from,to).parameters = {'g_SYN',.2,'tauRx',1,'tauDx',15,'E_SYN',-80,'IC_noise',1};

netmodeler(net); % rt_biosim(net); data=biosim(net); biosim_plots(data);

%specpath='/space/mdeh3/9/halgdev/projects/jsherfey/code/modeler/database';
%rt_biosim(net);      data=biosim(net);      biosim_plots(data);
%rt_biosim(specpath); data=biosim(specpath); biosim_plots(data);
