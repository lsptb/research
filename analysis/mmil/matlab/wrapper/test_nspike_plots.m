% test_nspike_plots.m
%
% Purpose: 
% Demonstrates the use of nspike_plots() to generate online plots of average waveforms.
%
% Required input: 
%   epoch_data      -- may contain multiple conditions
%
% Optional input:
%   channels        -- vector of channel indices
%   cfg             -- configuration structure with preprocessing options
%
% Formats:
%   nspike_plots(epoch_data, cfg, channels);
%   nspike_plots(epoch_data, channels);
%
% Created by Jason Sherfey on 29-Oct-2008

datafile = '/home/jsherfey/matlab/wip/nspike/epoch_data.mat';
load(datafile);

cfg = [];
cfg.detrend   = 'yes';
cfg.blc       = 'yes';
cfg.blcwindow = [-0.4 -0.1];

channels = [1:90];

nspike_plots(epoch_data,cfg,channels);


