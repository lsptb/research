function K = defweights(n,type,args)
if nargin<3, args={}; end
parms = mmil_args2parms(args,{'pdf','unifrnd',[],'a',0,[],'b',1,[],'span',1,[]},0);
nD = length(n);
if nD==2, Npre=n(1); Npost=n(2); else Npre=n; Npost=n; end
Nmax=max(n);

switch type
  case 'gaussian'
    sigma=(parms.span*Nmax);
    if nD==1
      Xpre = linspace(1,Nmax,Npre)'*ones(1,Npost)
      Xpost = (linspace(1,Nmax,Npost)'*ones(1,Npre))'
      K = scale*exp(-(Xpre-Xpost).^2/sigma^2);      
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
      %Dact = Dpq([0 0 1 1]);
      % Kernel
      K = exp(-Dpq([0 0 sigma sigma]).^2); % N^2 x N^2
    else
      K=ones(n.^2);
    end
    %{
    % Non-gaussian kernel component (ramp increasing from 0 at north to 1 at south)
    W = reshape(repmat(linspace(0,1,N)',[1 N])',[1 N^2]);
    W = repmat(W,[N^2 1]);
    % Kernel with Gaussian
    % Variable Kernels (parameterse are a function of eccentricity)
    Decc = Dfun(N^2/2-N/2,[0 0 1 1]); % distance b/w each cell & the fovea (center point)
    kvar = @(siglo,sighi) siglo+(2*Decc/(N*sqrt(2)))*(sighi-siglo);
    Dvaux = @(P) cellfun(@(x,u1,u2,s1,s2)Dfun(x,[u1 u2 s1 s2]),num2cell(1:N^2),...
                  num2cell(kvar(P(1),P(2))),num2cell(kvar(P(3),P(4))),...
                  num2cell(kvar(P(5),P(6))),num2cell(kvar(P(7),P(8))),'uniformoutput',false);
    Dvpq = @(P) reshape(cell2mat(Dvaux(P)),[N^2 N^2]); % dist b/w cells in pops p & q
    Kv = @(P) (((P(11)>0)*W + (P(11)==0)).*exp(-Dvpq(P(1:8)).^2))*P(9); % cell^2 x cell^2
    %}
  case 'distributed'
    
  otherwise
    K=ones(Nmax^2,Nmax^2);
end
