% Vers�o do Singleton (1969) utilizando la�o for

%x = [1 5 3 -1 -2 -3 1 -9 1 1 1 1 1 1 1 1 0 0 0 0 2 2 2 2 3 3 3 3 4 4 4 4];
t = -2:4/128:1.99999;x = 2*sin(6*pi*t)+8*sin(15*pi*t)+4*sin(10*pi/4*t);
%x = [1 1 1 1 1 1 1 1];


N = length(x);                  % Comprimento da amostra
N_2 = N/2;                      % Metada do comprimento da amostra
camadas = log2(N);              % N�mero de camadas 
W = exp(1j*2*pi*[0:N_2-1]/N);   % W_N com valores de e^j2pik/N

% Frequ�ncias digitais normalizadas para plot 
f = [0:1:(N-1)]/N;

X_0 = fft(x);
subplot(2,1,1)
stem(f,abs(X_0)/N,'k')
title('DFT calculada com FFT do MATLAB')
xlabel('\omega (x\pi rads/amostra)')
hold on


% Vetor de �ndices em nota��o bin�ria reversa
index = bi2de(de2bi(0:N-1,log2(N),'left-msb'))+1;

X_5 = x(index);                 % Entrada para algor�timo do Singleton
X_temp = zeros(1,N);            % Armazena resultados parciais das camada

% M�scara de bits para c�lculo do �ndice de W
w_index_mask = N-1;                     % Set os log2(N) primeiros bits
w_index_mask = w_index_mask*2^camadas;  % Desloca de log2(N)-1 bits

for c= 1:1:camadas
   for e = 1:2:N
        % C�lculos dos �ndices do elemento e dos W_N's
        i_x_t = bitshift(e,-1)+1;
        i_w = bitshift(bitand(e,w_index_mask),-1)+1;
        
        % C�lculo das amostras pares e �mpares
        X_temp(i_x_t)       = X_5(e)+W(i_w).*X_5(e+1);
        X_temp(i_x_t+N_2)   = X_5(e)-W(i_w).*X_5(e+1);
   end
   % Apontando os valores de sa�da para serem novas entradas
   X_5 = X_temp;
   % Desloca a m�scara para direita
   w_index_mask = bitshift(w_index_mask,-1);
end

subplot(2,1,2)
stem(f,abs(X_5)/N,'r')
title('DFT calculada com implementa��o pr�pria')
xlabel('\omega (x\pi rads/amostra)')

hold off
