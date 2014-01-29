function [M,nnf]=ts_calc_inverse_operator(G,SNR,C,R);
%function [M,nnf]=ts_calc_inverse_operator(G,[SNR],[C],[R]);
%
% ts_calc_inverse_operator: calculate dSPM inverse operator
%   using cortically constrained minimum norm solution
%
% Input:
%   G: forward gain matrix of sensor amplitudes for each dipole
%      (or dipole orientation)   (num_sensors x num_sources)
%   SNR: estimated signal to noise ratio {default = 10}
%   C: sensor covariance matrix {default = identity matrix}
%   R: source covariance matrix {default = scaled identity matrix}
%
% Output:
%   M: 'whitenend' inverse operator
%   nnf: vector of source scaling factors (noise normalization factor)
%
% Acknowledgements:
%   method from: Dale et al Neuron, Vol. 26, 55-67, April, 2000
%   notation based on: MNE manual by M. Hamalainen Chap 4
%   code modified from: mn_inverse_operator.m by M. Huang, Ph.D. 05/19/05
%
% created 08/02/05 by Don Hagler
% last modified 02/08/07 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
  help(mfilename);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=[];
nnf=[];

if(isempty(G))
  error('gain matrix G is empty');
end
[num_sensors,num_sources] = size(G);

if ~exist('SNR','var'), SNR=[]; end;
if ~exist('C','var'), C=[]; end;
if ~exist('R','var'), R=[]; end;

if isempty(SNR), SNR=10; end;
if isempty(C)
  C = eye(num_sensors);
end;
if isempty(R)
  R=sparse(num_sources,num_sources);
  for i=1:num_sources
    R(i,i)=1;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate C^-1/2 for prewhitening
[Uc,gammac]=eig(C);
rank_C=rank(C);
sqrt_inv_gammac=zeros(size(gammac));
for i=1:num_sensors
  if (num_sensors-i+1)<=rank_C
    sqrt_inv_gammac(i,i)=sqrt(1.0/gammac(i,i));
  end;
end;
Cinv_sqrt=sqrt_inv_gammac*Uc';

% prewhiten gain matrix
G=Cinv_sqrt*G;

% scale source covariance matrix so that trace(GRG')/trace(I)=1
tmp = trace(G*R*G');
if isnan(tmp)
  error('trace of GRG'' is NaN');
end;
R=num_sensors*R/tmp;

% Cholesky factorization for computational efficiency
Rc=chol(R);             

% calculate inverse operator with singular value decomposition
A=G*Rc;
[V,lamda,U]=svd(A',0);  % economy SVD preventing memory overflow
lamda=diag(lamda);
lamda_sq_SNR=1.0/SNR^2; % power SNR is needed for regularization
gamma=diag(lamda./(lamda.^2+lamda_sq_SNR));
gamma=sparse(gamma);
Vbar_gamma=Rc*V*gamma;
M=Vbar_gamma*U'*Cinv_sqrt; % this is the inverse operator

% noise sensitivity normalization
nnf=zeros(num_sources,1);
for k=1:num_sources
  nnf(k)=sqrt(sum(Vbar_gamma(k,:).^2));
end

return;

