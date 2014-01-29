%
% Hardware definition file for nstream data acquisition system
% Configured for mobile 256 channel NSpike system plus one dataglove
% for Thomas Thesen's lab at the NYU Comprehensive Epilepsy Center
%
% Questions: adam@cns.nyu.edu
%

% path to nstream binaries
path(path, '/home/nspike/svn_code/nsuite/nstream/pkg/mex');

%
%
% masterdsp configuration
nstream.hardware.nspike.masterdsp(1).ip = '10.1.2.10';
nstream.hardware.nspike.masterdsp(2).ip = '10.1.2.20';

%
%
% auxdsp configuration
nstream.hardware.nspike.auxdsps(1).ip  = '10.1.2.11';
nstream.hardware.nspike.auxdsps(2).ip  = '10.1.2.12';
nstream.hardware.nspike.auxdsps(3).ip  = '10.1.2.13';
nstream.hardware.nspike.auxdsps(4).ip  = '10.1.2.14';
nstream.hardware.nspike.auxdsps(5).ip  = '10.1.2.15';
nstream.hardware.nspike.auxdsps(6).ip  = '10.1.2.16';
nstream.hardware.nspike.auxdsps(7).ip  = '10.1.2.17';
nstream.hardware.nspike.auxdsps(8).ip  = '10.1.2.18';
nstream.hardware.nspike.auxdsps(9).ip  = '10.1.2.21';
nstream.hardware.nspike.auxdsps(10).ip = '10.1.2.22';
nstream.hardware.nspike.auxdsps(11).ip = '10.1.2.23';
nstream.hardware.nspike.auxdsps(12).ip = '10.1.2.24';
nstream.hardware.nspike.auxdsps(13).ip = '10.1.2.25';
nstream.hardware.nspike.auxdsps(14).ip = '10.1.2.26';
nstream.hardware.nspike.auxdsps(15).ip = '10.1.2.27';
nstream.hardware.nspike.auxdsps(16).ip = '10.1.2.28';

%
%
% dataglove configuration
nstream.hardware.dataglove(1).serial_port = '/dev/ttyS0';
