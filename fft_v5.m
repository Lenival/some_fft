% Versão do Singleton (1969) utilizando laço for

x = [1 5 3 -1 -2 -3 1 -9 1 1 1 1 1 1 1 1 0 0 0 0 2 2 2 2 3 3 3 3 4 4 4 4];
%x = [1 1 1 1 1 1 1 1];


X_0 = fft(x);
subplot(2,1,1)
stem(abs(X_0))
hold on



N = length(x);
N_2 = N/2;
camadas = log2(N);

W = exp(1j*2*pi*[0:N_2-1]/N);

% Vetor de índices
index = bi2de(de2bi(0:N-1,log2(N),'left-msb'))+1;

X_5 = x(index);

% Vetor para armazenar resultados de cada camada
X_temp = zeros(1,N);

% Máscara de bits para cálculo do índice de W
w_index_mask = N-1; % Set os log2(N) primeiros bits
w_index_mask = w_index_mask*2^camadas; % Desloca de log2(N)-1 bits

%index = reshape(index,[2,N/2]);
X_temp = zeros(1,N);
aux = 0;
for c= 1:1:camadas
   for e = 1:2:N
        [bitshift(e,-1) bitshift(e,-1)+N_2 e e+1 bitshift(bitand(e,w_index_mask),-1) w_index_mask]
        i_x_t = bitshift(e,-1)+1;
        i_w = bitshift(bitand(e,w_index_mask),-1)+1;
        
        X_temp(i_x_t)       = X_5(e)+W(i_w).*X_5(e+1);
        X_temp(i_x_t+N_2)   = X_5(e)-W(i_w).*X_5(e+1);
%         X_temp(bitshift(e,-1)+1)	= X_4(e)+W(bitshift(bitand(e,w_index_mask),-1)+1).*X_4(e+1);
%         X_temp(bitshift(e,-1)+N_2+1)= X_4(e)-W(bitshift(bitand(e,w_index_mask),-1)+1).*X_4(e+1);
   end
   % Apontando os valores de saída para serem novas entradas
   X_5 = X_temp;
   % Desloca a máscara para direita
   w_index_mask = bitshift(w_index_mask,-1);
end

subplot(2,1,2)
stem(abs(X_5))

hold off
