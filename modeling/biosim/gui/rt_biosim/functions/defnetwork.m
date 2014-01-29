function A = defnetwork(n,type,args)
if nargin<3, args={}; end
parms = mmil_args2parms(args,{'span',inf,[],'prob',1,[],'SelfConn',1,[]},0);
nD = length(n);
if nD==2, Npre=n(1); Npost=n(2); else Npre=n; Npost=n; end
Nmax=max(n);

switch type
  case 'fanout'
    fanout = Nmax*parms.span;   
    if nD==1
      Xpre = linspace(1,Nmax,Npre)'*ones(1,Npost)
      Xpost = (linspace(1,Nmax,Npost)'*ones(1,Npre))'
      A = abs(Xpre-Xpost)<=fanout;
    elseif nD==2
      N=Nmax;
      % Positions
      c=reshape(1:N^2,[N N])';
      xp=mod(c-1,N)+1;
      yp=xp';
      pos = @(k)  [xp(k) yp(k)];
      % Weighted Distance
      Dfun = @(ind,P) sqrt(((xp(c(ind))-xp(c(:))-P(1))).^2/P(3)^2 + ((yp(c(ind))-yp(c(:))-P(2))).^2/P(4)^2)';%-1
      Daux = @(P) cellfun(@(x)Dfun(x,P),num2cell(1:N^2),'uniformoutput',false);
      Dpq  = @(P) reshape(cell2mat(Daux(P)),[N^2 N^2]); % dist b/w cells in pops p & q
      % Absolute Distances
      Dact = Dpq([0 0 1 1]);
      % Kernel
      A = Dact<=fanout;
    else
      A=ones(n.^2);
    end    
  case 'random'
    p = parms.prob;
    A=(random('unif',0,1,[Nmax^2,Nmax^2])<=p);
  case 'smallworld'
    
  otherwise
    A=ones([Nmax^2,Nmax^2]);
end
A = A - (1-parms.SelfConn)*diag(diag(A)); % N^2 x N^2