function [Xn_k,D] = hdft_calc(x,M,L,bins)

N = length(x);

% How many hypothetical windows should exist before x(M)?
q_D = floor((M-1)/L);           % Number of hypothetical windows
X_start = mod(M-1,L)+1;         % Where should it starts

% Number of windows demanded
Q = floor((N-M)/(L))+1;
Xn_k = zeros(Q,M);

% Initializing bins vector
if nargin < 4
   bins = 0:1:M-1; 
end

% Computing a lookup table for twiddle factors
W_M_Lk = zeros(1,M);
for k = bins
    W_M_Lk(k+1) = exp(1j*(2*pi/M)*L*(k));
end

if max(abs(W_M_Lk)) > 1
    disp('One pole is out of unitary circle!');
end

%% HDFT UVT calculation
D = raw_uvt(x,N,M,L,Q,q_D,X_start,bins);
%D = uvt_radix2(x,N,M,L,Q,q_D,X_start,bins);

%% HDFT spectrogram

% The first window assumes that the firs X_n(k) is zero for all bins
Xn_k(1,1:M)=W_M_Lk.*D(1,0+1:M);

% Processing the remaining windows
for n_i = 2:1:(Q+q_D)        % n=0 is calculated before
    w_i = 0;
    for k = bins
        w_i = w_i + 1;
        Xn_k(n_i,k+1)=W_M_Lk(k+1)*(Xn_k(n_i-1,k+1)+D(n_i,k+1));
    end
end

