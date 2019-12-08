%clear;clc
%% Define reference signal x and your DFT
%x = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
%x =  [1 2 3 4 3 2 1 2 3 4 3  2  1  2  3  4  3  2  1  2  3  4  3  2  1];
%x = randn(1,512);
Fs = 2^10;t = 0:1/Fs:2;x = 1+sin(2*pi*20*t);
%x = 1:64*128;x(1:10*128) = sin(2*pi*x(1:10*128));
%x = [1 2 3 4 1 2 3 4 1 2 3 4];
%Fs = 10e3;t = 0:1/Fs:2;x = vco(sawtooth(2*pi*t,0.5),[0.1 0.4]*Fs,Fs);
%Fs = 2^16;t = 0:1/Fs:2;x = sawtooth(2*pi*512*t,0.75);
%Fs = 2^6;t = 0:1/Fs:2;x = sawtooth(2*pi*512*t,0.75);

N = length(x);
M = 256;            % Window length
L = 64;            % Size of window increment
S = M-L;          % Superposition length
q_D = floor((M-1)/L);             % Number of UVT before the start window
X_start = mod(M-1,L)+1;       % Start sample of UVT calculation


% Number of windows demanded
Q = floor((N-M)/(L))+1;
Xn_k = zeros(Q+q_D,M,2);

% Initializing with zero pading from x(-M+1) until x(-1)
%Xn_k(1,:,1) = fft([zeros(1,M-1) x(1)]);
%Xn_k(2,:,1) = fft([zeros(1,L-1) x(1:1+M_h)]);

% Padding zeroes x(n)

disp('-------------FFT--------------')
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


D_reff = zeros(Q+q_D,M);
D_reff(1,:) = exp(-1j*(2*pi/M)*L*(0:M-1)).*Xn_k(1,:,1);

for n_i = 2:Q+q_D
    D_reff(n_i,:) = exp(-1j*(2*pi/M)*L*(0:M-1)).*Xn_k(n_i,:,1)-Xn_k(n_i-1,:,1);
end


%% HDFT UVT calculation
disp('-------------HDFT--------------')
% k = 0:1:M-1;                % Frequency bin index
% W_M = exp(-1j*2*pi*L*k);    % Twiddle factors
W_M_1 = exp(1j*2*pi/M);
%exponentials = [];
% Calculating the first sample to 

% Calculating d(n) for zero initial conditions from x[-(M-1)] até x[1]
d = zeros(1,N);
d(1:M) = x(1:M);
d(M+1:N)=x(M+1:N)-x(1:N-M);

D = zeros(Q+q_D,M);
% Calculating D_n[k] to x[-(M+L)] até x[0] with null initial conditions 
% D(0+1,0+1:M)=d(1:M_h).*exp(1j*2*pi*((0-L+1).*(0:M-1))./M);

% sum(d(1:5).*W_M_1.^(-(3:7)*0))
for n_i = 1:1:q_D
    n_x = X_start + (n_i-1)*L;      % Index n of x(n) window 
    m_i = (n_x>=M)*(n_x-L)+1;
    %m_f = n;
    disp([m_i n_x])
    for k = 0:1:M-1         % q=0 is calculated before
        for m = m_i:1:n_x     % Matlab index related to n-m HDFT index
            D(n_i,k+1)=D(n_i,k+1)+d(m)*exp(1j*2*pi*(((n_x-m)-L+1)*k)/M);   
            %exponentials = [exponentials exp(1j*2*pi*(((n_x-m)-L+1)*k)/M)];% 
        end
    end
end

for n_i = q_D+1:1:(Q+q_D)
    %n_x = X_start + (n_i-1)*M_h;      % Index n of x(n) window
    n_x = M+(n_i-q_D-1)*L;
    m_i = (n_x-L)+1;
    %m_f = n;
    disp([m_i n_x])
    for k = 0:1:M-1         % q=0 is calculated before
        for m = m_i:1:n_x     % Matlab index related to n-m HDFT index
            D(n_i,k+1)=D(n_i,k+1)+d(m)*exp(1j*2*pi*(((n_x-m)-L+1)*k)/M); 
            %exponentials = [exponentials exp(1j*2*pi*(((n_x-m)-L+1)*k)/M)];% 
        end
    end
end


%% HDFT spectrogram

% Doing the first window
%Xn_k(1,1:M,2)=exp(1j*(2*pi/M)*L*(0:M-1)).*D(1,0+1:M);
Xn_k(q_D,1:M,2)= Xn_k(q_D,1:M,1);%Xn_k(2,1:M,2)= Xn_k(2,1:M,1);
%exponentials = [exponentials exp(1j*(2*pi/M)*L*(0:M-1))];
%Xn_k(q_D+1,1:M,2)=fft(x(1:1:M));


% Processing the remaining windows
for n_i = q_D+1:1:(Q+q_D)        % n=0 is calculated before
    for k = 0:1:M-1
        %n_l = (n-1)*L+1;
        Xn_k(n_i,k+1,2)=exp(1j*(2*pi/M)*L*k)*(Xn_k(n_i-1,k+1,2)+D(n_i,k+1));
        %exponentials = [exponentials exp(1j*(2*pi/M)*L*k)];
    end
end

D_diff = D - D_reff;
Diff_Xn_k = Xn_k(:,:,1)-Xn_k(:,:,2);
%% Graph plot to compare results
disp('-------------plots--------------')

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
surf(abs(Xn_k(q_D:Q+q_D,:,1)))
%title('SDFT com FFT do MATLAB')
zlabel('X_{n}(k)')
xlabel('n')
ylabel('k')
colormap(pink)    % change color map
shading interp    % interpolate colors across lines and faces

subplot(1,2,2)
% stem(abs(fft(x)/N))
% imagesc(abs(Xn_k(:,:,2)))
surf(abs(Xn_k(q_D:Q+q_D,:,2)))
zlabel('X_{n}(k)')
xlabel('n')
ylabel('k')
colormap(pink)    % change color map
%title('HDFT implementada')
shading interp    % interpolate colors across lines and faces

figure
mesh(abs(Diff_Xn_k(q_D:Q+q_D,:)))
zlabel('\epsilon')
xlabel('n')
ylabel('k')