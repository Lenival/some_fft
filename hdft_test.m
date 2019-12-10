clear;clc
%% Define reference signal x and your DFT
%x = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
%x =  [1 2 3 4 3 2 1 2 3 4 3  2  1  2  3  4  3  2  1  2  3  4  3  2  1];
%x = randn(1,2^18);
Fs = 2^10;t = 0:1/Fs:2;x = 1+sin(2*pi*200*t);
%x = [1 2 3 4 1 2 3 4 1 2 3 4];
%Fs = 10e3;t = 0:1/Fs:2;x = vco(sawtooth(2*pi*t,0.5),[0.1 0.4]*Fs,Fs);
%Fs = 2^16;t = 0:1/Fs:2;x = sawtooth(2*pi*512*t,0.75);
%Fs = 2^6;t = 0:1/Fs:2;x = sawtooth(2*pi*512*t,0.75);

N = length(x);
M = 64;            % Window length
L = 2;            % Size of window increment

% How many hypothetical windows should exist before x(M)?
q_D = floor((M-1)/L);           % Number of hypothetical windows
X_start = mod(M-1,L)+1;         % Where should it starts

% Number of windows demanded
Q = floor((N-M)/(L))+1;
Xn_k = zeros(Q+q_D,M,2);

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


%% HDFT function

[Xn_k(:,:,2),D] = hdft_calc(x,M,L);

%% HDFT Analysis

% Estimating D values to compare with those calculated with the HDFT
D_reff(1,:) = exp(-1j*(2*pi/M)*L*(0:M-1)).*Xn_k(1,:,1);

for n_i = 2:Q+q_D
    D_reff(n_i,:) = exp(-1j*(2*pi/M)*L*(0:M-1)).*Xn_k(n_i,:,1)-Xn_k(n_i-1,:,1);
end

D_diff = D - D_reff;
Xn_k_diff = Xn_k(:,:,1)-Xn_k(:,:,2);


%% Graph plot to compare results

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
surf(abs(Xn_k(1:Q+q_D,:,1)))
%title('SDFT com FFT do MATLAB')
zlabel('X_{n}(k)')
xlabel('n')
ylabel('k')
colormap(pink)    % change color map
shading interp    % interpolate colors across lines and faces

subplot(1,2,2)
% stem(abs(fft(x)/N))
% imagesc(abs(Xn_k(:,:,2)))
surf(abs(Xn_k(1:Q+q_D,:,2)))
zlabel('X_{n}(k)')
xlabel('n')
ylabel('k')
colormap(pink)    % change color map
%title('HDFT implementada')
shading interp    % interpolate colors across lines and faces

figure
mesh(abs(Xn_k_diff(1:Q+q_D,:)))
zlabel('\epsilon')
xlabel('n')
ylabel('k')