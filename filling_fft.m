% Definição de pulso de 10 amostras de amplitude 1

% Suporte de 10 amostras em vetor de 10 amostras
x = ones(10,1);
X = fftshift(fft(x));

% Pulso num vetor de 20 amostras onde 10 são zeros
x_0 = [zeros(5,1);ones(10,1);zeros(5,1)];
X_0 = fftshift(fft(x_0));

% Pulso num vetor de 16 amostras onde 6 são complementos cíclicos
x_1 = [ones(3,1);ones(10,1);ones(3,1)];
X_1 = fftshift(fft(x_1));

% Pulso num vetor de 16 amostras onde 6 são zeros
x_2 = [zeros(3,1);ones(10,1);zeros(3,1)];
X_2 = fftshift(fft(x_2));


stem(abs(X_2),'c')
hold on
stem(abs(X_0))
stem(abs(X),'k')
stem(abs(X_1),'m')
