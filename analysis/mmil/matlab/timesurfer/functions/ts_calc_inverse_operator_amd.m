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
%   M: inverse operator
%   nnf: vector of source scaling factors (noise normalization factor)
%
% Acknowledgements:
%   method from: Dale et al Neuron, Vol. 26, 55-67, April, 2000
%
% created:   08/02/2005 by Don Hagler
% last mod:  01/20/2008 by Anders Dale
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
  R=speye(num_sources,num_sources);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scale source covariance matrix so that trace(GRG')/trace(I)=1
tmp = trace(G*R*G');
if isnan(tmp)
  error('trace of GRG'' is NaN');
end;
R=num_sensors*R/tmp;

% calculate inverse operator
lamda_sq_SNR=1.0/SNR^2; % power SNR is needed for regularization
M = R*G'*inv(G*R*G' + lamda_sq_SNR*mean(diag(G*R*G'))/mean(diag(C))*C);

% noise sensitivity normalization
nnf=zeros(num_sources,1);
for k=1:num_sources
  nnf(k)=sqrt(M(k,:)*C*M(k,:)');
end

return;

