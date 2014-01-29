function stim=rtpulsestim(sz,points,t,t0,args)
% params:
  % points = multiple subscripts to nodes receiving input in sz-grid (returned by getpts)
  % args (key/value pairs listing params)
    % type (integer identifier)
      % 1: pulse onclick. params: dur, amp
      % 2: random dots. params: rho, diam
stim=zeros(sz);
if numel(points)==0 || numel(t0)==0
  stim = reshape(stim,[prod(sz) 1]);
  return
end
if isequal(size(points),[1 4]) % rectangle from getrect(); convert to (x,y)
  n=sz(1); p=sz(2); sp=points;
  sp(1) = max(floor(sp(1)),1); % xmin
  sp(2) = max(floor(sp(2)),1); % ymin
  sp(3) = min(ceil(sp(1)+sp(3)),p); % xmax
  sp(4) = min(ceil(sp(2)+sp(4)),n); % ymax
  js=sp(1):sp(3); nj=numel(js);
  is=sp(2):sp(4); ni=numel(is);
  npts=nj*ni; cnt=0;
  points=zeros(npts,2);
  for i=1:ni
    for j=1:nj
      cnt=cnt+1;
      points(cnt,2)=js(j);
      points(cnt,1)=is(i);
    end
  end
end
parms = mmil_args2parms(args,{'type',1,[],'dur',.1,[],'amp',1,[]},0);
target=sub2ind(sz,floor(points(:,1)),floor(points(:,2))); % linear indices

switch parms.type
  case 1 % pulse onclick
    t0(t0<(t-parms.dur))=[];
    nstim=length(t0); I=0;
    for k=1:nstim % kick up the stimulus for each event in the last dur
      I = I + parms.amp*((t-t0(k))<parms.dur); % could make this exp decay...
    end
    stim(target) = I;
end
stim = reshape(stim,[prod(sz) 1]); % return vectorized form

%{
stim.txt:
  args = [{'type',1,'dur',.1,'amp',1}] % brackets may be unecessary
    % note: 'args' will be turned into Elabel_Mlabel_args in buildmodel()
  term => rt_pulsestim([sqrt(N) sqrt(N)],points,t,tlast,args)
in rt_biosim():
  ...
  tlast=-inf; target=[];
  % then numeric integration...
    integrate(step)
    if stateevent, tlast(end+1)=tlastevent (ex. spike in target node or LFP)
    if onclick,    tlast(end+1)=tlastclick; [u,v] = getpts(fig); points=[u,v];
    if callback:optionselect
      arglabel=[Elabel '_' opttype '_args']; % ex) opttype='stim'
      args = eval(arglabel);
      [UPDATE ARGS W/ NEW OPTIONS; user modified val of opt (string)]
      args{find(cellfun(@(x)isequal(x,opt)))+1} = val;
      eval([arglabel '=args;']); % ex) for 'stim': opt='type','dur','amp'
note: opttype's in GUI should match mech file names and each given its own
subpanel in the set-param area of the GUI figure; opt's for a type are keys
with defaults set in the mech file var 'args' passed to auxfuncs (which
set the absolute defaults); user opt values are set w/ edit and radio controls.
%}