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

%% HDFT UVT calculation
D = raw_uvt(x,N,M,L,Q,q_D,X_start,bins);

%% HDFT spectrogram

% The first window assumes that the firs X_n(k) is zero for all bins
Xn_k(1,1:M)=exp(1j*(2*pi/M)*L*(0:M-1)).*D(1,0+1:M);

% Processing the remaining windows
for n_i = 2:1:(Q+q_D)        % n=0 is calculated before
    for k = bins
        Xn_k(n_i,k+1)=exp(1j*(2*pi/M)*L*k)*(Xn_k(n_i-1,k+1)+D(n_i,k+1));
    end
end

