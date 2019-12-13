clear;clc
%%
Fs = 2^12;t = 0:1/Fs:2;x = sawtooth(2*pi*512*t,0.75);M = 273;L = 87;

N = length(x);

% How many hypothetical windows should exist before x(M)?
q_D = floor((M-1)/L);           % Number of hypothetical windows
X_start = mod(M-1,L)+1;         % Where should it starts


% Number of windows demanded
Q = floor((N-M)/(L))+1;
Xn_k = zeros(Q+q_D,M,2);

% Uncoment if you want just a few bins
% bins = 13:27;

bins = 0:M-1;

% calculating all the hypothetical windows before x(M)
for n_i = 1:1:q_D
    n_x = X_start + (n_i-1)*L;      % Index n of x(n) window 
    m_i = (n_x>=M)*(n_x-L)+1;
    Xn_k(n_i,:,1)=fft([ zeros(1,M - n_x) x(m_i:n_x)]);
end

% Matlab's FFT spectrogram
for n_i = q_D+1:1:(Q+q_D)
    n_l = M+(n_i-q_D-1)*L;
    Xn_k(n_i,:,1)=fft(x(n_l-M+1:n_l));
end

if length(bins) < M
    P = zeros(Q+q_D,M);
    O = ones(Q+q_D,M);
    P(:,bins+1) = P(:,bins+1)+O(:,bins+1);
    Xn_k(:,:,1) = Xn_k(:,:,1).*P;
end



%% HDFT function

D = raw_uvt(x,N,M,L,Q,q_D,X_start,bins);
%D = uvt_radix2(x,N,M,L,Q,q_D,X_start,bins);
%D = uvt_dit_r2(x,N,M,L,Q,q_D,X_start,bins);

%% HDFT Analysis

% Estimating D values to compare with those calculated with the HDFT
iW_M_L = exp(-1j*(2*pi/M)*L*(0:M-1));
D_reff = zeros(Q+q_D,M);
D_reff(1,:) = iW_M_L.*Xn_k(1,:,1);

for n_i = 2:Q+q_D
    D_reff(n_i,:) = iW_M_L.*Xn_k(n_i,:,1)-Xn_k(n_i-1,:,1);
end

D_diff = zeros(Q+q_D,M);
D_diff = D - D_reff;

figure
mesh(abs(D_diff(1:Q+q_D,:)))
zlabel('\epsilon')
xlabel('k')
ylabel('n')

