function X = mtfilt(tapers, sampling, f0)
%  MTFILT Generates a bandpass filter using the multitaper method
%
%  X = MTFILT(TAPERS, SAMPLING, F0)
%

%  Author:  Bijan Pesaran 04/15/98
%


pr_op = dp_proj(tapers, sampling, f0);
N = size(pr_op,1);
X = zeros(1,2.*N);
pr = pr_op*pr_op';
for t = 0:N-1 
	X(t+1:t+N) = X(t+1:N+t)+pr(:,N-t)'; 
end
X = X./N;

