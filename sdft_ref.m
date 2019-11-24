%x = [1 2 3 4 3 2 1 2 3 4 3 2 1 2 3 4 3 2 1 2 3 4 3 2 ];
%x = randn(1,64*128);
%x = 1:64*128;x(1:10*128) = sin(2*pi*x(1:10*128));
%x = [1 2 3 4 1 2 3 4 1 2 3 4];
%Fs = 10e3;t = 0:1/Fs:2;x = vco(sawtooth(2*pi*t,0.5),[0.1 0.4]*Fs,Fs);
Fs = 10000;t = 0:1/Fs:2;x = sawtooth(2*pi*512*t,0.75);

N = length(x);
M = 2048;
L = 1024;

% Quantidade de janelas calculadas
Q = floor((N-M)/L+1);
Xn_k = zeros(Q,M);


for q = 1:1:Q
    n_l = (q-1)*L+1;
    Xn_k(q,:)=fft(x(n_l:n_l+M-1));
end

%HeatMap(abs(Xn_k))
surf(abs(Xn_k)/M)
colormap(pink)    % change color map
shading interp    % interpolate colors across lines and faces

figure;
plot(abs(fft(x)/N))