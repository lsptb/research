function [model,IC,functions,auxvars,sys,Sodes,Svars] = buildmodel2(spec,varargin)
% get input structure into the right form

if ~isfield(spec,'connections')
  if isfield(spec,'mechanisms')
    if isfield(spec,'N')
      for i=1:length(spec)
        spec(i).multiplicity = spec(i).N;
      end
      spec = rmfield(spec,'N');
    end
    tmp.entities = spec;
    if isfield(spec,'files')
      tmp.files = spec(1).files;
    end
    spec=tmp; clear tmp    
  end
  spec.connections.label='';
  spec.connections.mechanisms={};
  spec.connections.parameters={};
end
if isfield(spec,'cells') && ~isfield(spec,'entities')
  spec.entities = spec.cells;
  spec = rmfield(spec,'cells');
end
if isfield(spec,'connections') && any(size(spec.connections)<length(spec.entities))
  n=length(spec.entities);
  spec.connections(n,n)=spec.connections(1);
  spec.connections(n,n).label=[];
  spec.connections(n,n).mechanisms=[];
  spec.connections(n,n).parameters=[];
end
if isfield(spec,'files')
  spec.files = unique(spec.files);
end
parms = mmil_args2parms( varargin, ...
                         {  ...
                            'logfid',1,[],...
                            'override',[],[],...
                            'dt',.01,[],...
                            'nofunctions',1,[],...
                            'verbose',1,[],...
                         }, false);
% note: override = {label,field,value,[arg]; ...}
fileID = parms.logfid;

if ischar(spec)
  spec=loadspec(spec);
end
Elabels = {spec.entities.label}; % Entity labels
Clabels = {spec.connections.label};

% override parameters in specification
if ~isempty(parms.override)
  o = parms.override;
  [nover,ncols] = size(o);
  for k = 1:nover
    l = o{k,1}; f = o{k,2}; v = o{k,3}; 
    if ncols>3, a = o{k,4}; else a = []; end
    if ~ischar(l) || ~ischar(f), continue; end
    if ismember(l,Elabels), type='entities';
    elseif ismember(l,Clabels), type='connections';
    else continue; 
    end
    n = strmatch(l,{spec.(type).label},'exact');
    if isequal(f,'parameters')
      if isempty(spec.(type)(n).(f)), continue; end
      matched = find(cellfun(@(x)isequal(v,x),spec.(type)(n).(f)));
      if ~isempty(matched)
        spec.(type)(n).(f){matched+1} = a;
      else
        spec.(type)(n).(f){end+1} = v;
        spec.(type)(n).(f){end+1} = a;
      end
    else
      spec.(type)(n).(f) = v;
    end
  end
  clear o nover ncols l f v a n
end

I=find([spec.entities.multiplicity]~=0);
spec.entities = spec.entities(I);
spec.connections = spec.connections(I,I);
Elabels = {spec.entities.label}; % Entity labels
Clabels = {spec.connections.label};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load complete spec (including all mech structs)
% specpath='/space/mdeh3/9/halgdev/projects/jsherfey/model/adv/new_spec/pfc_nmda';
% spec=loadspec(specpath);

sys = spec;
NE = [spec.entities.multiplicity];
N = length(spec.entities);
if issubfield(spec,'simulation.sim_dt')
  dt=spec.simulation.sim_dt;
else
  dt=parms.dt;
end

% combine intrinsic and connection mechanisms per entity; load mech models
mechtype={}; % 0=connection, 1=intrinsic
mechsrc=[]; mechdst=[];
for i=1:N % loop over entities
  clear m1 m2 p1 p2 inp1 inp2
  m1 = spec.entities(i).mechanisms;
  [p1{1:length(m1)}]=deal(spec.entities(i).parameters);
  inp1=num2cell(i*ones(1,length(m1)));
  m2={spec.connections(:,i).mechanisms};
  p2={spec.connections(:,i).parameters};
  sel=find(~cellfun(@isempty,m2));
  inp2=num2cell(sel);
  m2=m2(sel);
  p2=p2(sel);
  mechtype{i}=ones(size(m1));
  mechsrc{i}=i*ones(size(m1)); mechdst{i}=i*ones(size(m1)); 
  if ~isempty(m2) % there are connections to this entity
    tmpm={}; tmpp={}; tmpinp={};
    for j=1:length(m2)
      if ~iscell(m2{j})
        tmp=splitstr(m2{j},' '); % multiple mechanisms are separated by a space
      else
        tmp=m2{j};
      end
      if size(tmp,2)>1, tmp=tmp'; end
      tmpm=cat(2,tmpm,tmp');
      tmpp=cat(2,tmpp,cat(2,repmat(p2(j),[1 length(tmp)])));
      tmpinp=cat(2,tmpinp,cat(2,repmat(inp2(j),[1 length(tmp)])));
      mechtype{i} = [mechtype{i} zeros(1,length(tmp))];
      mechsrc{i} = [mechsrc{i} sel(j)*ones(1,length(tmp))];
      mechdst{i} = [mechdst{i} i*ones(1,length(tmp))];
    end
    m2=tmpm;
    p2=tmpp;
    inp2=tmpinp;
    %mechtype{i} = [mechtype{i} zeros(size(m2))];
  end
  if size(m2,1)>1, m2=m2'; end
  mechs=cat(2,m1,m2); % combine intrinsic and connection mechanisms
  sys.entities(i).mechanisms = mechs;
  sys.entities(i).parameters = cat(2,p1,p2);
  sys.entities(i).inputs = cat(2,inp1,inp2);
  if ~iscell(sys.entities(i).dynamics)
    sys.entities(i).dynamics = splitstr(sys.entities(i).dynamics,' ');
  end
  M=length(mechs);
  % load mechanism models from text files
  for j=1:M
    ML=mechs{j};    
    if issubfield(sys.entities(i),'mechs.label')
      if ismember(ML,{sys.entities(i).mechs.label}) && length(sys.entities(i).mechs)>=j
        continue; % mechanism model already provided by user
      end
    end
    idx = regexp(spec.files,sprintf('(^%s|/%s|\\\\%s).txt',ML,ML,ML));
    idx = find(~cellfun(@isempty,idx));
    if length(idx) ~= 1
      fprintf('Looking in known mech list for %s\n',ML);
      [jnk,knownfiles] = get_mechlist;
      idx = regexp(knownfiles,sprintf('(^%s|/%s|\\\\%s).txt',ML,ML,ML));
      idx = find(~cellfun(@isempty,idx));
      if length(idx) ~= 1
        error('Failed to find a distinct specification file for mechanism: %s',ML);
      else
        file = knownfiles{idx};
        fprintf('Found mech file. Using: %s\n',file);
      end
    else
      file = spec.files{idx};
    end
    this = parse_mech_spec(file,[]); % includes mech parameters
    try
      this.label = ML;
      sys.entities(i).mechs(j) = this;
    catch
      sys.entities(i).mechs(j) = this;
    end
  end
  % remove old mech models
  if issubfield(sys.entities(i),'mechs.label')
    test = setdiff({sys.entities(i).mechs.label},sys.entities(i).mechanisms);
    if ~isempty(test)
      sys.entities(i).mechs(ismember({sys.entities(i).mechs.label},test))=[];
    end
  end
end

% Compute total number of state vars, params, funcs, expressions
nGvar=sum(cellfun(@length,{sys.entities.dynamics}));
nMvar=sum(cellfun(@sum,arrayfun(@(x)sum(cellfun(@(y)size(y,1),{x.mechs.odes})),sys.entities,'unif',0)));
nGparm=0; nMparm=0;
for k=1:N
  nMparm=nMparm+sum(arrayfun(@(x)length(fieldnames(x.params)),sys.entities(k).mechs));
  nGparm=nGparm+sum(cellfun(@(x)length(x),sys.entities(k).parameters)/2);
end
nfunc=sum(cellfun(@sum,arrayfun(@(x)sum(cellfun(@(y)size(y,1),{x.mechs.functions})),sys.entities,'unif',0)));
nexpr=sum(cellfun(@sum,arrayfun(@(x)sum(cellfun(@(y)size(y,1),{x.mechs.auxvars})),sys.entities,'unif',0)));
nsubst=sum(cellfun(@sum,arrayfun(@(x)sum(cellfun(@(y)size(y,1),{x.mechs.substitute})),sys.entities,'unif',0)));
nmech=sum(cellfun(@length,{sys.entities.mechs}));
nvar=nGvar+nMvar;
nparm=nGparm+nMparm;
%[nGvar nMvar nvar nGparm nMparm nparm nfunc nexpr nmech]

% Preallocate indicator/mapping vectors
Svars = cell(nvar,3); % {label, prefix_label, indices, IC}
Sodes = cell(nvar,1); % F, var'=F (differential equation)
Spop = zeros(nvar,1);
Smech = zeros(nvar,1);
Stype = zeros(nvar,1);
Pdata = cell(nparm,2); % {label, value}
Ppop = zeros(nparm,1);
Pmech = zeros(nparm,1);
Ptype = zeros(nparm,1);
Cexpr = cell(nexpr,3); % {label, prefix_label, expression}
Cpop = zeros(nexpr,1);
Cmech = zeros(nexpr,1);
Hfunc = cell(nfunc,3); % {label, prefix_label, function}
Hpop = zeros(nfunc,1);
Hmech = zeros(nfunc,1);
Tsubst = cell(nsubst,2); % {entity ode label, mech term}
Tpop = zeros(nsubst,1);
Tmech = zeros(nsubst,1);
Minputs = cell(nmech,1); % entity #s of input mechs
Mlabels = cell(nmech,1);
Mid=zeros(nmech,1);
% Create list of state variables Svars (label, prefix_label, indices) & odes Sodes over all populations and mechanisms
%      Sg: pop_var. Smk: pop_mech_var. I=cnt+(1:NE). odes = Og & Fmk
%      + associated lists of entity ids Spop (entity indices), and vector Stype of type indicators (0 if in Sg; m if in Sm)
% Create list of parameters Pdata (label, value) over all populations and mechanisms
%      + associated lists of Ppop, Pmech, and Ptype (0 if in Pg; m if in Pm); note: do not explicitly store Pg & Pm separately
% Create list of expressions: Cexpr (LHS: label, prefix_label. RHS: expression) over all mechanisms|populations
%      Cmk: pop_mech_label
%      + associated lists of Cpop (entity index) and Cmech (mech index)
% Create list of functions: Hfunc (LHS: label, prefix_label, inputs) over all mechanisms|populations
%      Hmk: pop_mech_label
%      + associated lists of Hpop (entity index) and Hmech (mech index)

scnt=0; % state var label index
stateindx=0; % state var vector indices
pcnt=0; ccnt=0; hcnt=0; tcnt=0; mcnt=0;
EL = {sys.entities.label};
for i=1:N
  E=sys.entities(i); n=length(E.dynamics); I=scnt+(1:n);
  tmp=regexp(E.dynamics,'\w+''','match');
  tmp=cellfun(@(x)x{1}(1:end-1),tmp,'unif',0);
  [Svars{I,1}]=deal(tmp{:});                % varlabel
  for j=1:length(I)
    Svars{I(j),2}=[EL{i} '_' tmp{j}];       % E_varlabel
    Svars{I(j),3}=stateindx+(1:NE(i));      % state vector indices
    Svars{I(j),4}= ['[' num2str(zeros(1,NE(i))) ']'];           % initial conditions
    stateindx=stateindx+NE(i);
  end
  tmp=regexp(E.dynamics,'=(.)*','match');
  tmp=cellfun(@(x)x{1}(2:end),tmp,'unif',0);
  [Sodes{I}]=deal(tmp{:});                  % ODEs
  Spop(I)=i; 
  Smech(I)=0;
  Stype(I)=0;                               % var scope (0=global)
  scnt=scnt+n;
  for m=1:length(E.mechanisms)
    mcnt=mcnt+1;
    Mid(mcnt)=mcnt;
    M=E.mechs(m);
    Minputs{mcnt}=E.inputs{m};
    Mlabels{mcnt}=E.mechanisms{m};
    if Minputs{mcnt}~=i
      prefix = [EL{Minputs{mcnt}} '_' EL{i} '_' Mlabels{mcnt}];
    else
      prefix = [EL{i} '_' Mlabels{mcnt}];
    end
    if ~isempty(E.parameters{m})  % USER PARAMETERS
      n=length(E.parameters{m})/2; I=pcnt+(1:n);
      tmp=E.parameters{m}(1:2:end); [Pdata{I,1}]=deal(tmp{:});
      tmp=E.parameters{m}(2:2:end); [Pdata{I,2}]=deal(tmp{:});
      %[Pdata{I,1}]=deal(E.parameters{m}{1:2:end});  % user param key
      %[Pdata{I,2}]=deal(E.parameters{m}{2:2:end});  % user param value
      Ppop(I)=i;
      Pmech(I)=mcnt;
      Ptype(I)=0;
      pcnt=pcnt+n;
    end
    if ~isempty(M.statevars)      % MECH STATE VARS
      tmp=M.statevars;
      n=length(tmp); I=scnt+(1:n);
      [Svars{I,1}]=deal(tmp{:});                            % varlabel
      for j=1:length(I)
%         if Minputs{mcnt}~=i
%           Svars{I(j),2}=[prefix '_' tmp{j}]; % Pre_E_M_varlabel
%         else
          Svars{I(j),2}=[prefix '_' tmp{j}]; % E_M_varlabel
%         end
        Svars{I(j),3}=stateindx+(1:NE(i));                  % state vector indices
        Svars{I(j),4}=M.ic{j};                              % initial conditions
        stateindx=stateindx+NE(i);
      end      
      [Sodes{I}]=deal(M.odes{:});                           % ODEs
      Spop(I)=i; 
      Smech(I)=mcnt;
      Stype(I)=1;                                           % var scope (0=global)
      scnt=scnt+n;
    end
    if ~isempty(M.params)         % DEFAULT MECH PARAMETERS
      key=fieldnames(M.params); val=struct2cell(M.params);
      n=length(key); I=pcnt+(1:n);
      [Pdata{I,1}]=deal(key{:});  % default param key
      [Pdata{I,2}]=deal(val{:});  % default param value
      Ppop(I)=i;
      Pmech(I)=mcnt;
      Ptype(I)=1;
      pcnt=pcnt+n;
    end
    if ~isempty(M.auxvars)        % MECH EXPRESSIONS (label,prefix_label,expression)
      LHS=M.auxvars(:,1); RHS=M.auxvars(:,2);
      n=length(LHS); I=ccnt+(1:n);
      [Cexpr{I,1}]=deal(LHS{:});              % exprlabel
      [Cexpr{I,3}]=deal(RHS{:});              % expression
      Cpop(I)=i;
      Cmech(I)=mcnt;
      for j=1:length(I)
        Cexpr{I(j),2}=[prefix '_' LHS{j}];    % E_M_exprlabel
      end      
      ccnt=ccnt+n;
    end
    if ~isempty(M.functions)      % MECH FUNCTIONS  (label,prefix_label,function)
      LHS=M.functions(:,1); RHS=M.functions(:,2);
      n=length(LHS); I=hcnt+(1:n);
      [Hfunc{I,1}]=deal(LHS{:});              % funclabel
      [Hfunc{I,3}]=deal(RHS{:});              % function
      Hpop(I)=i;
      Hmech(I)=mcnt;
      for j=1:length(I)
        Hfunc{I(j),2}=[prefix '_' LHS{j}];    % E_M_exprlabel
      end      
      hcnt=hcnt+n;
    end
    if ~isempty(M.substitute)     % SUBSTITUTIONS {entity ode label, mech term}
      LHS=M.substitute(:,1); RHS=M.substitute(:,2);
      n=length(LHS); I=tcnt+(1:n);
      [Tsubst{I,1}]=deal(LHS{:});              % subst label
      [Tsubst{I,2}]=deal(RHS{:});              % term to insert
      Tpop(I)=i;
      Tmech(I)=mcnt;
      tcnt=tcnt+n;
    end
  end
end

% [Mid,Mlabels,Minputs] % mech ids, labels, input populations
% Sodes: (differential equation = F: var'=F )
% Svars: {label, prefix_label, indices, IC}
% Pdata: {label, value}
% [Stype  Spop Smech] % variables
% [Ptype  Ppop Pmech] % parameters
% Cexpr, [Cpop Cmech] % expressions:        {label, prefix_label, expression}
% Hfunc, [Hpop Hmech] % functions:          {label, prefix_label, function}
% Tsubst,[Tpop Tmech] % term substitutions: {label, mech_term}
Hfunc0=Hfunc; Cexpr0=Cexpr; Sodes0=Sodes; Tsubst0=Tsubst; Svars0=Svars;
for m=1:nmech
  f=Hfunc(Hmech==m,3); 
  e=Cexpr(Cmech==m,3); 
  o=Sodes(Smech==m); 
  t=Tsubst(Tmech==m,2);
  ic=Svars(Smech==m,4);
  % global user params: (expressions, functions, odes, terms)
  old=Pdata(Pmech==m & Ptype==0,1); new=Pdata(Pmech==m & Ptype==0,2);
  [f,e,o,t,ic]=substitute(old,new,f,e,o,t,ic);
  % default mech params | same pop & mech: (expressions, functions, odes, terms)
  old=Pdata(Pmech==m & Ptype==1,1); new=Pdata(Pmech==m & Ptype==1,2);
  [f,e,o,t,ic]=substitute(old,new,f,e,o,t,ic);
  % reserved keywords: (expressions, functions, odes, terms)
  old={'Npre','N[1]','Npost','N[0]','Npop','dt'};
  k0=unique(Ppop(Pmech==m)); n0=NE(k0);
  if m>0, k1=Minputs{m}(1); else k1=k0; end
  n1=NE(k1);
  new={n1,n1,n0,n0,n0,dt};
  [f,e,o,t,ic]=substitute(old,new,f,e,o,t,ic);
  % prefixes: expressions: (expressions, functions, odes, terms)
  old=Cexpr(Cmech==m,1); new=Cexpr(Cmech==m,2);
  [f,e,o,t]=substitute(old,new,f,e,o,t);
  % prefixes: functions: (functions, odes, terms)
  old=Hfunc(Hmech==m,1); new=Hfunc(Hmech==m,2);
  [f,o,t]=substitute(old,new,f,o,t);
  % prefixes: vars (Sg(E) => Sg(~E) => Sm(E) => S~m(E)): (odes, terms)
  E=unique(Ppop(Pmech==m));
  old=Svars(Spop==E & (Smech==m | Stype==0),1); 
  new=Svars(Spop==E & (Smech==m | Stype==0),2); 
  old2={}; new2={}; old3={}; new3={};
  for k=1:length(old), old2{k}=[old{k} '[0]']; new2{k}=new{k}; end
  for k=1:length(old)
    old3{k}=[old{k} '[1]']; 
    new3{k}= [EL{k1} new{k}(find(new{k}=='_',1,'first'):end)]; 
  end
  [o,t]=substitute(old3,new3,o,t);
  [o,t]=substitute(old2,new2,o,t);
  [o,t]=substitute(old,new,o,t);
  for k=1:length(old), old2{k}=[old{k} 'post']; new2{k}=new{k}; end
  for k=1:length(old)
    old3{k}=[old{k} 'pre']; 
    new3{k}= [EL{k1} new{k}(find(new{k}=='_',1,'first'):end)]; 
  end
  [o,t]=substitute(old3,new3,o,t);
  [o,t]=substitute(old2,new2,o,t);  
  Hfunc(Hmech==m,3)=f; 
  Cexpr(Cmech==m,3)=e; 
  Sodes(Smech==m)=o; 
  Tsubst(Tmech==m,2)=t;
  Svars(Smech==m,4)=ic;
end
E=unique(Spop(Stype==0));
for e=1:length(E)
  idx=(Stype==0 & Spop==E(e));
  o=Sodes(idx);
  old=Svars(idx,1); new=Svars(idx,2); 
  o=substitute(old,new,o);
  old=Pdata(Ptype==0,1); new=Pdata(Ptype==0,2);
  o=substitute(old,new,o);
  old={'Npost','N[0]','Npop','dt'}; n0=NE(e);
  new={n0,n0,n0,dt};
  o=substitute(old,new,o);
  Sodes(idx)=o;
end

if parms.nofunctions
  % Substitute functions into functions
  keep_going=1; cnt=0;
  while keep_going
    keep_going=0; cnt=cnt+1;
    for t=1:size(Hfunc,1)
      target=Hfunc{t,3};
      % match function labels between functions and their dependencies
      pat = @(s)sprintf('(\\W+%s)|(^%s)\\(',s,s); % function label pattern
      funcinds=find(~cellfun(@isempty,cellfun(@(s)regexp(target,pat(s)),Hfunc(:,2),'unif',0)));
      for f=1:length(funcinds)
        keep_going=1;
        ind = funcinds(f);
        submatch = regexp(target,[Hfunc{ind,2} '\([\w\s,]*\)'],'match');      
        %submatch = regexp(target,[Hfunc{ind,2} '\(.*\)'],'match');
        subvars = regexp(strrep(submatch,Hfunc{ind,2},''),'\([a-zA-Z]\w*\)','match');
        subvars = [subvars{:}];
        subvars = cellfun(@(s)strrep(s,'(',''),subvars,'unif',0);
        subvars = cellfun(@(s)strrep(s,')',''),subvars,'unif',0);
        str = Hfunc{ind,3};
        vars = regexp(str,'@\([\s\w,]*\)','match');
        expr = strtrim(strrep(str,vars,''));
        vars = regexp(vars,'\w+','match');
        vars = [vars{:}];
        if length(vars)~=length(subvars)
          error('coding error in buildmodel2 when inserting full expressions into terms Tmk');
        end
        expr = substitute(vars,subvars,expr);
        target = strrep(target,submatch{1},['(' expr{1} ')']); 
      end
      Hfunc{t,3}=target;
    end
  end  
  
  % Substitute functions into ODEs
  for t=1:size(Sodes,1)
    target=Sodes{t,1};
    % match function labels between functions and their dependencies
    pat = @(s)sprintf('(\\W+%s)|(^%s)\\(',s,s); % function label pattern
    funcinds=find(~cellfun(@isempty,cellfun(@(s)regexp(target,pat(s)),Hfunc(:,2),'unif',0)));
    for f=1:length(funcinds)
      ind = funcinds(f);
      submatch = regexp(target,[Hfunc{ind,2} '\([\w\s,]*\)'],'match');      
      %submatch = regexp(target,[Hfunc{ind,2} '\(.*\)'],'match');
      subvars = regexp(strrep(submatch,Hfunc{ind,2},''),'\([a-zA-Z]\w*\)','match');
      subvars = [subvars{:}];
      subvars = cellfun(@(s)strrep(s,'(',''),subvars,'unif',0);
      subvars = cellfun(@(s)strrep(s,')',''),subvars,'unif',0);
      str = Hfunc{ind,3};
      vars = regexp(str,'@\([\s\w,]*\)','match');
      expr = strtrim(strrep(str,vars,''));
      vars = regexp(vars,'\w+','match');
      vars = [vars{:}];
      if length(vars)~=length(subvars)
        error('coding error in buildmodel2 when inserting full expressions into terms Tmk');
      end
      expr = substitute(vars,subvars,expr);
      target = strrep(target,submatch{1},['(' expr{1} ')']); 
    end
    Sodes{t,1}=target;
  end  
  
  % Substitute functions into enitity-dynamic-substitution terms
  for t=1:size(Tsubst,1)
    % match function labels to substitution terms
    %pat = @(s)sprintf('(\\W+%s)|(^%s\\W+)',s,s);
    pat = @(s)sprintf('(\\W+%s)|(^%s)\\(',s,s); % function label pattern
    funcinds=find(~cellfun(@isempty,cellfun(@(s)regexp(Tsubst{t,2},pat(s)),Hfunc(:,2),'unif',0)));
    for f=1:length(funcinds)
      ind = funcinds(f);
      submatch = regexp(Tsubst{t,2},[Hfunc{ind,2} '\(.*\)'],'match');
      subvars = regexp(strrep(submatch,Hfunc{ind,2},''),'\w+','match');
      subvars = [subvars{:}];
      str = Hfunc{ind,3};
      vars = regexp(str,'@\([\s\w,]*\)','match');
      expr = strtrim(strrep(str,vars,''));
      vars = regexp(vars,'\w+','match');
      vars = [vars{:}];
      if length(vars)~=length(subvars)
        error('coding error in buildmodel2 when inserting full expressions into terms Tmk');
      end
      expr = substitute(vars,subvars,expr);
      Tsubst{t,2} = strrep(Tsubst{t,2},submatch{1},['(' expr{1} ')']); 
    end
  end
end

% Insert terms Tmk into entity dynamics
for t=1:size(Tsubst,1)
  idx = (Stype==0 & Spop==Tpop(t));
  o = Sodes(idx);
  %old=Tsubst{t,1}; new=['(' Tsubst{t,2} ')+' old];
  old=Tsubst{t,1}; new=['((' Tsubst{t,2} ')+' old ')'];
  Sodes(idx)=substitute(old,new,o);
end

% Evaluate ICs and determine state vector indices
stateindx=0;
for i=1:nvar
  s=Svars{i,1};
  ic=Svars{i,4};
  Pg=Pdata(Ptype==0 & Ppop==Spop(i),:);
  if any(find(cellfun(@(x)isequal(x,[s '_IC']),Pg(:,1))))
    ind2=find(cellfun(@(x)isequal(x,[s '_IC']),Pg(:,1)));
    icval=Pg{ind2(1),2};
    if numel(icval)==1
      ic = sprintf('[%s]',num2str(ones(1,NE(Spop(i)))*icval));    
    elseif numel(icval)==NE(Spop(i))
      ic = sprintf('[%s]',num2str(icval));
    elseif ischar(icval)
      ic = eval(icval);
    end
  end    
  if ischar(ic)
    ic=eval(ic);
  elseif numel(ic)==1
    ic=repmat(ic,[NE(Spop(i)) 1]);
  end
  if size(ic,1)<size(ic,2), ic=ic'; end
  if any(find(cellfun(@(x)isequal(x,[s '_IC_noise']),Pg(:,1))))
    ind2=find(cellfun(@(x)isequal(x,[s '_IC_noise']),Pg(:,1)));
    ic=ic+Pg{ind2(1),2}.*rand(size(ic));
  end  
  Svars{i,4}=ic;
  Svars{i,3}=stateindx+(1:length(ic));
  stateindx=stateindx+length(ic);
end
IC=cat(1,Svars{:,4});

% Substitute state vector indices
old=Svars(:,2);
new=cell(size(old));
for k=1:length(new)
  new{k}=sprintf('X(%g:%g)',Svars{k,3}(1),Svars{k,3}(end));
end
SodesVec=substitute(old,new,Sodes);
old=unique(Tsubst(:,1)); clear new
[new{1:length(old)}]=deal('0');
SodesVec=substitute(old,new,SodesVec);%,Sodes

% Prepare model string for evaluation
model = SodesVec;
for k=1:nvar, if model{k}(end) ~= ';', model{k} = [model{k} ';']; end; end
model = ['@(t,X) [' [model{:}] '];'];

% How to use:
if 0
  for k = 1:size(Cexpr,1)
    eval( sprintf('%s = %s;',Cexpr{k,2},Cexpr{k,3}) );
  end
  for k = 1:size(Hfunc,1)
    eval( sprintf('%s = %s;',Hfunc{k,2},Hfunc{k,3}) );
  end
  model = eval(model);
  Y=model(0,IC);
end

functions=Hfunc(:,[2 3 1]);
auxvars=Cexpr(:,[2 3 1]);

% Collect equations in system struct
% For spec.entities(i) & spec.connections(i,j): auxvars=Cexpr(:,2:3); functions=Hfunc(:,2:3)
% Set spec.variables.entity=Spop; spec.variables.labels=Svars(:,2)
% For each entity i: 
% 	spec.entities(i).var_index=cat(1,Svars{Spop==i,3})
% 	spec.entities(i).var_list=cat(1,Svars{Spop==i,2})
% 	spec.entities(i).orig_var_list=cat(1,Svars{Spop==i,1})
% 	also update: odes, ode_labels
sys.variables.entity=[];
for i=1:N
  sys.entities(i).auxvars = auxvars(Cpop==i,:);
  sys.entities(i).functions = functions(Hpop==i,:);
  sys.entities(i).odes = Sodes(Spop==i);
  sys.entities(i).ode_labels = Svars(Spop==i,2);
  sys.entities(i).var_index=cat(2,Svars{Spop==i,3});
  nhere=length(sys.entities(i).var_index);
  sys.entities(i).var_list={};%cell(1,length(sys.entities(i).var_index));
  sys.entities(i).orig_var_list={};%cell(1,length(sys.entities(i).var_index));
  ind=find(Spop==i);
  for j=1:length(ind)
    k=ind(j); nthis=length(Svars{ismember(Svars(:,2),Svars(k,2)),3});
    tmp=repmat(Svars(k,2),[1 nthis]);
    sys.entities(i).var_list      = {sys.entities(i).var_list{:} tmp{:}}; %repmat( Svars(Spop==i,2);
    tmp=repmat(Svars(k,1),[1 nthis]);
    sys.entities(i).orig_var_list = {sys.entities(i).orig_var_list{:} tmp{:}}; %sys.entities(i).orig_var_list=Svars(Spop==i,1);    
  end
  sys.variables.entity=[sys.variables.entity i*ones(1,nhere)];
  if any(mechtype{i}==0) % connection mechanisms
    sys.entities(i).connection_mechanisms = sys.entities(i).mechanisms(mechtype{i}==0);
    sys.entities(i).connection_parameters = sys.entities(i).parameters(mechtype{i}==0);
    sys.entities(i).connection_mechs = sys.entities(i).mechs(mechtype{i}==0);
    connis=find(mechtype{i}==0);
    for j=1:length(connis)
      ii=mechsrc{i}(connis(j));
      jj=mechdst{i}(connis(j));
      sys.connections(ii,jj).parameters = sys.entities(i).connection_parameters{j};
      if isempty(sys.connections(ii,jj).parameters) % pull default params from mech structure
        p = sys.entities(i).connection_mechs(j).params;
        keys = fieldnames(p); 
        vals = struct2cell(p); 
        tmp=[keys(:) vals(:)]';
        sys.connections(ii,jj).parameters = [tmp(:)]';
      end
      if ~isfield(sys.connections,'mechs') || isempty(sys.connections(ii,jj).mechs)
        sys.connections(ii,jj).mechs = sys.entities(i).connection_mechs(j);
      else
        sys.connections(ii,jj).mechs(end+1) = sys.entities(i).connection_mechs(j);
      end
    end
  end
  if any(mechtype{i}==1) % intrinsic mechanisms
    sys.entities(i).intrinsic_parameters = sys.entities(i).parameters(mechtype{i}==1);
    sys.entities(i).mechanisms = sys.entities(i).mechanisms(mechtype{i}==1);
    sys.entities(i).parameters = sys.entities(i).parameters{find(mechtype{i}==1,1,'first')};
    sys.entities(i).mechs = sys.entities(i).mechs(mechtype{i}==1);
  end
end
tmp={sys.entities.var_list};
sys.variables.labels = [tmp{:}];
sys.model.functions = functions;
sys.model.auxvars = auxvars;
sys.model.ode = model;
sys.model.IC = IC;

if parms.verbose
  % Print model info
  sgind=find(Stype==0);
  fprintf(fileID,'\nModel Description\n----------------------------------\n\n');
  fprintf(fileID,'Specification files:\n');
  for f = 1:length(spec.files)
    fprintf(fileID,'%s\n',spec.files{f});
  end
  fprintf(fileID,'\nPopulations:');
  for i = 1:N
    fprintf(fileID,'\n%-6.6s (n=%g):\n',EL{i},NE(i));
    ind=find(Spop==i & Stype==0);
    for j=1:length(ind)
      fprintf(fileID,'\tdynamics: %s'' = %s\n',Svars{ind(j),1},Sodes0{ind(j)});
    end
    M=spec.entities(i).mechanisms;
    P=spec.entities(i).parameters;
    if ~isempty(M)
      fprintf(fileID,'\tmechanisms: %s',M{1});
      for j = 2:length(M)
        fprintf(fileID,', %s',M{j});
      end
    end
    if length(P)>=2
      fprintf(fileID,'\n\tparameters: %s=%g',P{1},P{2});
      for j = 2:length(P)/2
        fprintf(fileID,', %s=%g',P{2*j-1},P{2*j});
      end
    end
  end
  fprintf(fileID,'\n\nConnections:\n');
  fprintf(fileID,'%-10.6s\t',' '); 
  for i = 1:N, fprintf(fileID,'%-10.6s\t',EL{i}); end; fprintf(fileID,'\n');
  for i = 1:N
    fprintf(fileID,'%-10.6s\t',EL{i});
    for j = 1:N
      if isempty(spec.connections(i,j).mechanisms)
        fprintf(fileID,'%-10.6s\t','-');
      else
        try
          fprintf(fileID,'%-10.10s\t',[spec.connections(i,j).mechanisms{:}]);
        catch
          fprintf(fileID,'%-10.10s\t',spec.connections(i,j).mechanisms);
        end
      end
    end
    fprintf(fileID,'\n');
  end
  fprintf(fileID,'\nConnection parameters:');
  CL={spec.connections.label};
  CP={spec.connections.parameters};
  for i = 1:length(CL)
    if isempty(CP{i}), continue; end
    fprintf(fileID,'\n%s: ',strrep(CL{i},'-','->'));
    if length(CP{i})>=2
      fprintf(fileID,'\n\tparameters: %s=%g',CP{i}{1},CP{i}{2});
      for j = 2:length(CP{i})/2
        fprintf(fileID,', %s=%g',CP{i}{2*j-1},CP{i}{2*j});
      end
    end
  end
  fprintf(fileID,'\n\nModel Equations:\n----------------\n');
  fprintf(fileID,'ODEs:\n');
  for i = 1:nvar
    fprintf(fileID,'\t%-20s = %s\n',[Svars{i,2} ''''],Sodes{i});
  end
  fprintf(fileID,'\nInitial Conditions:\n');
  for i = 1:nvar
    for j=1:length(Svars{i,4})
      if j==1
        fprintf(fileID,'\t%-20s=[',[Svars{i,2} '(0)']);
      end
      fprintf(fileID,'%.4g, ',Svars{i,4}(j));
    end
    fprintf(fileID,']\n');
  end
  fprintf(fileID,'\n');
  fprintf(fileID,'\nMatlab-formatted (copy and paste to repeat simulation):\n%%-----------------------------------------------------------\n');
  fprintf(fileID,'%% Auxiliary variables:\n');
  for i = 1:size(Cexpr,1)
    fprintf(fileID,'\t%-20s = %s;\n',Cexpr{i,2},Cexpr{i,3});
  end
  fprintf(fileID,'\n%% Anonymous functions:\n');
  for i = 1:size(Hfunc,1)
    %fprintf(fileID,'\t%-20s = %-40s\t%% (%-9s)\n',Hfunc{i,2},[Hfunc{i,3} ';'],functions{i,3});
    fprintf(fileID,'\t%-20s = %-40s\n',Hfunc{i,2},[Hfunc{i,3} ';']);
  end
  fprintf(fileID,'\n%% ODE Handle, ICs, integration, and plotting:\nF = %s\n',model);
  fprintf(fileID,'IC = [%s];\n',num2str(IC'));
  fprintf(fileID,'\n[t,y]=ode23(F,[0 100],IC);   %% numerical integration\nfigure; plot(t,y);           %% plot all variables/functions\n');
  if nvar>=1,fprintf(fileID,'try legend(''%s''',strrep(Svars{1,2},'_','\_')); end
  if nvar>=2,for k=2:nvar,fprintf(fileID,',''%s''',strrep(Svars{k,2},'_','\_')); end; end
  if nvar>=1,fprintf(fileID,'); end\n'); end
  fprintf(fileID,'%%-----------------------------------------------------------\n\n');
end
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
% default mech params | same pop & mech
  % expressions, functions, odes, terms
% reserved keywords
  % expressions, functions, odes, terms
% prefixes: expressions
  % expressions, functions, odes, terms
% prefixes: functions
  % functions, odes, terms
% prefixes: vars (Sg(E) => Sg(~E) => Sm(E) => S~m(E))
  % odes, terms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES & EQUATIONS

% Substitution procedure:
% During substitution, to avoid overwriting by substituting labels that are substrings of other labels: use substitution map by (1) finding matches, (2) storing the string to insert, (3) a unique identifier id, (4) and temporarily substituting the match by @(id). Then, after finding all do substitutions.

% Hfunc(:,3),Cexpr(:,3),Tsubst(:,2),Sodes
% [Stype  Spop Smech] % variables
% [Ptype  Ppop Pmech] % parameters
% Cexpr, [Cpop Cmech] % expressions:        {label, prefix_label, expression}
% Hfunc, [Hpop Hmech] % functions:          {label, prefix_label, function}
% Tsubst,[Tpop Tmech] % term substitutions: {label, mech_term}

% Substitute state vars into mech expression/functions/odes/substitution terms:
% 1. any E:Sg in E expression/function/ode/term without brackets is replaced by E_label (E_s)
%      (global var: within-entity)
% 2. any var[#] in E expression/function/ode/term is replaced by inputPop(#)_var if var in inputPop(#):Sg
%      note: not an option to substitute inputPop(#):Smk into different pop. (i.e. mech vars between-entity)
%      (global var: between-entity)
% 3. in E:m, any E:Sm (same mechanism) in expression/function/ode/term is replaced by E_m_label
%      (mech var: within-entity & within-mech)
% 4. in E:m, any var label (exposed in different mechanism E:Smk) in expression/function/ode/term is replaced by E_mk_label
%      (mech var: within-entity & between-mech)

% Update mech function/expression labels (within-E): substitute Cmk-label & Hmk-label by prefix_label in Cmk(RHS) & Hmk(RHS)
%      E:Cmk LHS prefix_label into Cmk(RHS), Hmk(RHS), Fmk, Tmk(?)
%      E:Hmk LHS prefix_label into Hmk(RHS), Fmk, Tmk(?)

% Insert terms Tmk into entity dynamics
% Prepare model string for evaluation
% Prepare expressions/functions for evaluation: insert prefix_label in Cmk(LHS) and Hmk(LHS)

function varargout = substitute(old,new,varargin)
  % INPUTS:
  % old = cell array of keys (strings)
  % new = cell array of strings or numeric values
  % cellmats: {cellmat, cellmat, ...}
  %   where cellmat: {x11,x12,x13,...; x21,x22,x23,...; ...}, x*=string
  varargout = varargin;
  if isempty(old) || isempty(new) || isempty(varargin), return; end
  cellmats = varargin;
  if ~iscell(old), old = {old}; end
  if ~iscell(new), new = {new}; end
  if size(old,1)>1, old = old'; end
  if iscellstr(old)
    [l,I] = sort(cellfun(@length,old),2,'descend');
    old = old(I);
    new = new(I);
  end
  % loop over cellmat elements: replace old=>new
  for i = 1:numel(cellmats)
    cellmat = cellmats{i};
    if ~iscell(cellmat), cellmat={cellmat}; end
    cells = cellmat(:);
    for j = 1:numel(cells)
      if ~ischar(cells{j}), continue; end      
      for k = 1:numel(old)
        this = cells{j};
        % handle reserved words first
        if isequal(old{k},'Npre') && isempty(regexp(cells{j},'[^A-Za-z]+Npre')), continue; end
        if isequal(old{k},'Npost') && isempty(regexp(cells{j},'[^A-Za-z]+Npost')), continue; end
        key = old{k};
        l = length(key);
        key = strrep(key,'[','\[');
        key = strrep(key,']','\]');
        if isnumeric(new{k})
          val=sprintf('(%g)',new{k});
        elseif ischar(new{k})
          val=new{k};
        end
        inds=regexp(this,['\<' key '\>']);
        if isempty(inds), continue; end
        if inds(1)>1
          tmp=this(1:inds(1)-1);
        else
          tmp=[];
        end
        for c=1:length(inds)
          tmp=[tmp val];
          if c<length(inds) && length(this)>=(inds(c)+l)
            tmp=[tmp this(inds(c)+l:inds(c+1)-1)];
          end
        end
        if inds(c)+l<=length(this)
          tmp=[tmp this(inds(c)+l:end)];
        end
        cells{j}=tmp;
      end
    end
    cellmats{i} = reshape(cells,size(cellmat));
  end
  varargout = cellmats;

