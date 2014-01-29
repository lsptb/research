% erpave() - perform epoch selection plus simple artifact rejection, 
%            baseline-zeroing, and/or block averaging
%
% Usage:  >> [dataout,accept,reject,accidx,rejidx]  ...
%             = erpave(data,framelist,epochframes,epochoffset);
%          or = erpave(data,framelist,epochframes,epochoffset,chanlist, ...
%                                 baseframes,rejectlimits,epochsonly,verbose);
% Inputs:
%   data         = raw data matrix (chans,framestotal)
%   framelist    = list of epoch-centered indices into data matrix (provide the
%                  center frame index for each epoch.
%   epochframes  = number of time points in each epoch
%   epochoffset  = time points from framelist index to start of epoch (often < 0) 
%   chanlist     = list of channel numbers in output data (default|0 -> all)
%   baseframes   = vector: the mean of these baseline frames 
%                  is removed from each channel   {default|0 = none}
%   rejectlimits = out-of-bounds rejection limits | (0        -> no rejection, 
%                  epochs with out-of-bounds vals | [n]       -> +/-n bounds
%                  are rejected from averaging.   | [min max] -> low/high bounds)
%   epochsonly   = if non-0, dataout is selected epochs [epoch epoch ...] 
%                  instead of average {default|0 = dataout is average only}
%   verbose      = turn verbose on/off (1/0) {default on}
%
% Outputs:
%   dataout      = averaged data [length(chanlist),framestot] or accepted epochs 
%   accept       = list of accepted framelist values (nsums = length(accept))
%   reject       = list of rejected framelist values 
%   accidx       = list of accepted framelist indices (nsums = length(accept))
%   reject       = list of rejected framelist indices 
%
% Authors: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1997 
%
% Example: dataout = erpave(EEG.data, [100 600 1100], 300, -50);
%          extract 3 epochs in the data at frames 100, 600, 1100. The
%          length of each epoch is 300 points and their range is
%          [50 350] for the first one, [550 850] and [1050 1350] for the
%          others. THe 3 epochs are averaged before begin returned


% Copyright (C) 1-14-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: erpave.m,v $
% Revision 1.2  2002/12/05 22:57:38  arno
% adding example
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 1-23-97 fixed baseline-zero option -sm
% 2-04-97 improved arg testing and verbosity -sm
% 2-07-97 tested row vector args -sm
% 4-04-97 shortened name to erpave() -sm
% 12-20-98 inserted missing initializations, cleaned usage msg, added verbose -sm
% 11-24-99 added accidx, rejidx -sm
% 01-25-02 reformated help & license -ad 

function [dataout,accept,reject,accidx,rejidx] = erpave(data,framelist,epochframes,epochoffset,chanlist,baseframes,rejectlimits,epochsonly,verbose)

if nargin < 4
    fprintf('erpave() - needs at least four arguments!\n');
    return
end

nchans =    size(data,1);
framestot = size(data,2);
accept = [];
reject = [];
rejidx = [];
accidx = [];

if nargin < 9
    verbose = 1;
end
if nargin < 8,
    epochsonly = 0; % default: dataout is average of selected epochs
end

if nargin < 7, 
    rejectlimits = [-1e25 1e25];
end
if rejectlimits == 0;
    rejectlimits = [-1e25 1e25];
elseif length(rejectlimits) == 1
    if rejectlimits > 0,
        rejectlimits = [-1*rejectlimits rejectlimits];
    else
        rejectlimits = [rejectlimits -1*rejectlimits];
    end
elseif rejectlimits(2) <= rejectlimits(1),
  tmp= rejectlimits(2); 
  rejectlimits(2) = rejectlimits(1); % make [min,max]
  rejectlimits(1) = tmp;
end

if nargin < 6
	baseframes = 0;		% default no baseline zeroing
end
if baseframes ~= 0,
  if size(baseframes,1)>1
    baseframes = baseframes'; % make it a row vector
  end
  if max(baseframes) > epochframes
fprintf('erpave() - max baseframes frame (%d) not in epoch.\n',...
                 max(baseframes));
      return
  end
  for b=baseframes
    if b<1 | b > framestot
fprintf('erpave() - baseframes index %d out of data.\n',b);
      return
    end
  end
end

if nargin < 5  
    chanlist = [1:nchans];
elseif chanlist == 0
    chanlist = [1:nchans];
end

if size(chanlist,1) > 1
	chanlist = chanlist';	% make it a row vector
end

if max(chanlist) > nchans,
 fprintf('erpave() - data has only %d channels! ');
 fprintf('Cannot ask for channel %d\n',nchans,max(chanlist));
 return
elseif min(chanlist) <1
 fprintf('erpave() - chanlist entries must be > 0! ');
end

if length(chanlist) > size(data,1)
    fprintf('erpave() - data has only %d channels. Cannot average on %d!\n',size(data,1),length(chanlist));
    return
end
    
if epochsonly~=0,
  basemeans = zeros(1,length(chanlist));
  dataout = zeros(length(chanlist),epochframes*length(framelist));
end
    
if length(framelist) < 1
    fprintf('erpave() - framelist is empty!\n');
    return
end

if size(framelist,1) > 1
	framelist = framelist'; % make it a row vector
end


%%%%%%%%%%%%%%%%%%%%%%%%%% print header decription %%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose ~= 0
    fprintf('\nerpave(): ');
	fprintf('Selecting and processing %d (%d,%d) epochs',...
                           length(framelist),length(chanlist),epochframes); 
  if rejectlimits(1) ~= -1e25
    fprintf(', using artifact reject limits');
    fprintf(' [%g %g]',rejectlimits(1),rejectlimits(2));
  end
    fprintf('\n');
	fprintf('             ');
  if epochsonly==0,
    fprintf(' Returning average of accepted epochs');
  else
    fprintf(' Returning accepted epochs');
  end
  if baseframes~=0,
	fprintf(' with baseline means removed');
  end
  fprintf('\n'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Artifact rejection:

   errorno=0;
   nsums = 0;
   accept=[];
   n=1;

   for f=framelist, %%%%%%%%%%%%%%%%%% foreach epoch %%%%%%%%%%%%%%%%%%%%%%%%%

    if f+epochoffset > framestot | f+epochoffset < 1,
if verbose ~= 0
 fprintf('  epoch on frame %d rejected (epoch start out of bounds)\n',f)
end
            reject = [reject f];
            rejidx = [rejidx n];
            errorno = 1;
    elseif f+epochoffset+epochframes-1>framestot|f+epochoffset+epochframes-1<1,
if verbose ~= 0
 fprintf('  epoch on frame %d rejected (epoch end out of bounds)\n',f)
end
            reject = [reject f];
            rejidx = [rejidx n];
            errorno = 1;
    else
        epoch = data(chanlist,f+epochoffset:f+epochoffset+epochframes-1);
    end

    if errorno==0 % foreach allowable epoch ...
      if rejectlimits(1) ~= -1e25 % if meaningful artifact reject limits set
        mx = max(max(epoch));
        mn = min(min(epoch));
        if mx>=rejectlimits(2)
if verbose ~= 0
  fprintf('  epoch on frame %d rejected (out of bounds value %g)\n',f,mx)
end
            reject = [reject f];
            rejidx = [rejidx n];
            errorno = 131;
        elseif mn<=rejectlimits(1),
if verbose ~= 0
  fprintf('  epoch on frame %d rejected (out of bounds value %g)\n',f,mn)
end
            reject = [reject f];
            rejidx = [rejidx n];
            errorno = 132;
        end
      end
      if errorno == 0,
        accept = [accept f];
        accidx = [accidx n];
        if epochsonly ~=0,   % if output is epochs only - no average
            dataout(:,nsums*epochframes+[1:epochframes]) = epoch; 
                                        % return epochs if epochsonly set
            if baseframes ~= 0,
               basemeans = basemeans + mean(epoch(:,baseframes)');
            end
        end
        nsums = nsums+1;
      end % epoch accept
    end % epoch test 
    errorno=0;

    fprintf('.')
    if rem(n,64)==0
      fprintf('%d\n',n)
    end
    n=n+1;
  end % framelist f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if rem(n,64)~=0
    fprintf('\n')
  end

  if epochsonly ~=0,   % if output is epochs only - no average
    dataout(:,length(accept)*epochframes+1:end) = []; % truncate reserved epochs
    if baseframes ~= 0,
      fprintf('Subtracting baseline means: ');
       basemeans = basemeans/length(accept);
       if verbose ~= 0
         fprintf('  Baseframes channel means subtracted: ');
         for bc=1:length(chanlist)
           fprintf('%g ',basemeans(bc));
         end; fprintf('\n');
       end
       basesub = basemeans'*ones(1,epochframes); 
       for e=1:length(accept), % subtract an epoch at a time to save memory
          dataout(:,(e-1)*epochframes+1:e*epochframes) ...
            = dataout(:,(e-1)*epochframes+1:e*epochframes) - basesub;
       end
    end
  end

  if exist('accept')==0
	accept = [];
  end
  if exist('reject')==0
	reject = [];
  end

  % Average selected data epochs:

  if epochsonly == 0,  % if averaging requested...
    if length(accept) > 0 % if epochs to average ....
if verbose ~= 0
  fprintf('  Total of %d epochs accepted, %d rejected.\n',...
                                    length(accept),length(reject));
end
      dataout = zeros(length(chanlist),epochframes);
      nsums = length(accept);
      for f=accept    
        epoch = data(chanlist,f+epochoffset:f+epochoffset+epochframes-1);
        dataout = dataout + 1/nsums*epoch;
      end
      if baseframes ~= 0,    % if baseline-zero
         basemeans = mean(dataout(:,baseframes)');
         if verbose ~= 0
           fprintf('  Baseframes channel means subtracted: ');
           for bc=1:length(chanlist)
             fprintf('%g ',basemeans(bc));
           end; fprintf('\n');
         end
         dataout = dataout - basemeans'*ones(1,epochframes); 
      end
    end % if accepted epochs
  end % if averaging requested

  if length(accept) == 0,
     fprintf('  All %d requested epochs rejected!\n',length(framelist));
     return
  end
