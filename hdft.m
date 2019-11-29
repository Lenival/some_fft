clear;clc
%% Define reference signal x and your DFT
% x = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
%x =   [1 2 3 4 3 2 1 2 3 4 3  2  1  2  3  4  3  2  1  2  3  4  3  2  1];
%x = randn(1,64*128);
%x = 1:64*128;x(1:10*128) = sin(2*pi*x(1:10*128));
%x = [1 2 3 4 1 2 3 4 1 2 3 4];
%Fs = 10e3;t = 0:1/Fs:2;x = vco(sawtooth(2*pi*t,0.5),[0.1 0.4]*Fs,Fs);
Fs = 2048;t = 0:1/Fs:2-1/Fs;x = sawtooth(2*pi*512*t,0.75);

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
for n_i = 3:1:Q
    n_l = (n_i-1)*L+1;
    Xn_k(n_i,:,1)=fft(x(n_l-M+1:n_l));
end




%% HDFT UVT calculation
% k = 0:1:M-1;                % Frequency bin index
% W_M = exp(-1j*2*pi*L*k);    % Twiddle factors
W_M_1 = exp(1j*2*pi/M);

% Calculating d(n) for zero initial conditions from x[-(M-1)] até x[1]
d = zeros(1,N);
d(1:M) = x(1:M);
d(M+1:N)=x(M+1:N)-x(1:N-M);
% 
% for n = M+1:1:N
%     for i = L-1:-1:0
%         %n-i
%         d(n-i+1)= x(n-i+1)-x(n-i-M+1);
%     end
% end

D = zeros(Q,M);
% Calculando D_n[k] considerando condições iniciais nulas x[-(M+L)] até x[0]
D(0+1,0+1:M)=d(1).*exp(1j*2*pi*((0-L+1).*(0:M-1))./M);

for n_i = 2:1:Q
    n = (n_i-1)*L;
    for k = 0:1:M-1         % q=0 is calculated before
        for m = 0:1:L-1
            D(n_i,k+1)=D(n_i,k+1)+d(n-m+1)*exp(1j*2*pi*((m-L+1)*k)/M);      % 
        end
    end
end

%% HDFT spectrogram

% for n_i = M+L:L:N-L
%     for k = 0:1:M-1         % q=0 is calculated before
%         for m = 0:1:L-1
%             D(n_i+1,k+1)=D(n_i+1,k+1)+d(n_i-m+1)*exp(1j*2*pi*((m-L+1)*k)/M);      % 
%         end
%     end
% end
% Do directly the first and second FFT
% Xn_k(0+1,:,2)=fft(x(1:M));
% Xn_k(0+2,:,2)=fft(x(L+1:M+L));
% Xn_k(0+3,:,2)=fft(x(M+1:2*M));

% Do the first window
Xn_k(1,1:M,2)=exp(1j*(2*pi/M)*L*(0:M-1)).*D(1,0+1:M);

% Processing the remaining windows
for n_i = 2:1:Q        % n=0 is calculated before
    for k = 0:1:M-1
        %n_l = (n-1)*L+1;
        Xn_k(n_i,k+1,2)=exp(1j*(2*pi/M)*L*k)*(Xn_k(n_i-1,k+1,2)+D(n_i,k+1));
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