clear;clc
%% Define reference signal x and your DFT
% x = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
x =   [1 2 3 4 3 2 1 2 3 4 3  2  1  2  3  4  3  2  1  2  3  4  3  2  1];
%x = randn(1,64*128);
%x = 1:64*128;x(1:10*128) = sin(2*pi*x(1:10*128));
%x = [1 2 3 4 1 2 3 4 1 2 3 4];
%Fs = 10e3;t = 0:1/Fs:2;x = vco(sawtooth(2*pi*t,0.5),[0.1 0.4]*Fs,Fs);
%Fs = 2048;t = 0:1/Fs:2-1/Fs;x = sawtooth(2*pi*512*t,0.75);

N = length(x);
M = 8;
L = 4;

% Number of windows demanded
Q = floor((N-M)/L+1)+2;
Xn_k = zeros(Q,M,2);

% Initializing with zero pading from x(-M+1) until x(-1)
Xn_k(1,:,1) = fft([zeros(1,M-1) x(1)]);
Xn_k(2,:,1) = fft([zeros(1,M-L-1) x(1:1+L)]);

% Matlab's FFT spectrogram
for n = 3:1:Q
    n_l = (n-1)*L+1
    Xn_k(n,:,1)=fft(x(n_l-M+1:n_l));
end




%% HDFT spectrogram
% k = 0:1:M-1;                % Frequency bin index
% W_M = exp(-1j*2*pi*L*k);    % Twiddle factors
W_M_1 = exp(1j*2*pi/M);

% Calculating the UVT
d = zeros(1,N);

% Calculando d(n) e D_n[k] p/condições iniciais nulas x[-(M+L)] até x[0]
d(1:M+L) = x(1:M+L);



for n = M+L:L:N-1
    for i = L-1:-1:0
        %n-i
        d(n-i+1)= x(n-i+1)-x(n-i-M+1);
    end
end

D = zeros(N,M);
% Calculando D_n[k] considerando condições iniciais nulas x[-(M+L)] até x[0]
D(0+1,0+1:7+1)=d(1).*exp(1j*2*pi*((0-L+1).*(0:7))./M);
for n = L:L:M+L
    for k = 0:1:M-1         % q=0 is calculated before
    	for m = 0:1:L-1
            D(n+1,k+1)=D(n+1,k+1)+d(n-m+1)*exp(1j*2*pi*((m-L+1)*k)/M);      % 
        end
    end
end
for n = M+L:L:N-L
    for k = 0:1:M-1         % q=0 is calculated before
        for m = 0:1:L-1
            D(n+1,k+1)=D(n+1,k+1)+d(n-m+1)*exp(1j*2*pi*((m-L+1)*k)/M);      % 
        end
    end
end
% Do directly the first and second FFT
% Xn_k(0+1,:,2)=fft(x(1:M));
% Xn_k(0+2,:,2)=fft(x(L+1:M+L));
% Xn_k(0+3,:,2)=fft(x(M+1:2*M));

% Do the first window
Xn_k(1,1:8,2)=exp(1j*2*pi*L*0:7).*D(1,1:8);

% Processing the remaining windows
for n = 1:1:Q-1         % n=0 is calculated before
    for k = 0:1:M-1
        Xn_k(n+1,k+1,2)=exp(1j*2*pi*L*k)*(Xn_k(n,k+1,2)+D(n*L+1,k+1));
    end
end


%% Graph plot to compare results
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
surf(abs(Xn_k(:,:,1))/M)
colormap(pink)    % change color map
shading interp    % interpolate colors across lines and faces

subplot(1,2,2)
% stem(abs(fft(x)/N))
surf(abs(Xn_k(:,:,2))/M)
colormap(pink)    % change color map
shading interp    % interpolate colors across lines and faces