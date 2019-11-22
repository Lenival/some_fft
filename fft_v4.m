% Versão do Singleton (1969) utilizando laço for

x = [1 5 3 -1 -2 -3 1 -9];
%x = [1 1 1 1 1 1 1 1];


X_0 = fft(x);
subplot(3,1,1)
stem(abs(X_0))
hold on



N = length(x);
N_2 = N/2;
camadas = log2(N);

W = exp(1j*2*pi*[0:N_2-1]/N);

% Vetor de índices
index = bi2de(de2bi(0:N-1,log2(N),'left-msb'))+1;

X_4 = x(index);
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
        
        X_temp(i_x_t)       = X_4(e)+W(i_w).*X_4(e+1);
        X_temp(i_x_t+N_2)   = X_4(e)-W(i_w).*X_4(e+1);
%         X_temp(bitshift(e,-1)+1)	= X_4(e)+W(bitshift(bitand(e,w_index_mask),-1)+1).*X_4(e+1);
%         X_temp(bitshift(e,-1)+N_2+1)= X_4(e)-W(bitshift(bitand(e,w_index_mask),-1)+1).*X_4(e+1);
   end
   % Apontando os valores de saída para serem novas entradas
   X_4 = X_temp;
   % Desloca a máscara para direita
   w_index_mask = bitshift(w_index_mask,-1);
end

subplot(3,1,2)
stem(abs(X_4))

X_temp = zeros(1,N);
% Camada 1
% 0 4 ; 2 6 ; 1 5 ; 3 7
% 0 4 ; 0 4 ; 0 4 ; 0 4 --> 0 0 ; 0 0 ; 0 0 ; 0 0 

X_temp(0+1) = X_5(0+1)+W(0+1).*X_5(1+1);
X_temp(4+1) = X_5(0+1)-W(0+1).*X_5(1+1);

X_temp(1+1) = X_5(2+1)+W(0+1).*X_5(3+1);
X_temp(5+1) = X_5(2+1)-W(0+1).*X_5(3+1);

X_temp(2+1) = X_5(4+1)+W(0+1).*X_5(5+1);
X_temp(6+1) = X_5(4+1)-W(0+1).*X_5(5+1);

X_temp(3+1) = X_5(6+1)+W(0+1).*X_5(7+1);
X_temp(7+1) = X_5(6+1)-W(0+1).*X_5(7+1);

X_5 = X_temp;

% Camada 2
% 0 2 ; 4 6 ; 1 3 ; 5 7
% 0 4 ; 2 6 ; 0 4 ; 2 6 -->  0 0 ; 2 2 ; 0 0 ; 2 2

X_temp(0+1) = X_5(0+1)+W(0+1).*X_5(1+1);
X_temp(4+1) = X_5(0+1)-W(0+1).*X_5(1+1);

X_temp(1+1) = X_5(2+1)+W(0+1).*X_5(3+1);
X_temp(5+1) = X_5(2+1)-W(0+1).*X_5(3+1);

X_temp(2+1) = X_5(4+1)+W(2+1).*X_5(5+1);
X_temp(6+1) = X_5(4+1)-W(2+1).*X_5(5+1);

X_temp(3+1) = X_5(6+1)+W(2+1).*X_5(7+1);
X_temp(7+1) = X_5(6+1)-W(2+1).*X_5(7+1);


X_5 = X_temp;

% Camada 3
% 0 1 ; 4 5 ; 2 3 ; 6 7
% 0 4 ; 1 5 ; 2 6 ; 3 7 -->  0 0 ; 1 1 ; 2 2 ; 3 3

X_temp(0+1) = X_5(0+1)+W(0+1).*X_5(1+1);
X_temp(4+1) = X_5(0+1)-W(0+1).*X_5(1+1);

X_temp(1+1) = X_5(2+1)+W(1+1).*X_5(3+1);
X_temp(5+1) = X_5(2+1)-W(1+1).*X_5(3+1);

X_temp(2+1) = X_5(4+1)+W(2+1).*X_5(5+1);
X_temp(6+1) = X_5(4+1)-W(2+1).*X_5(5+1);

X_temp(3+1) = X_5(6+1)+W(3+1).*X_5(7+1);
X_temp(7+1) = X_5(6+1)-W(3+1).*X_5(7+1);



X_5 = X_temp;

subplot(3,1,3)
stem(abs(X_5))


%
%for c = 1:1:log2(N)
%    c
%end
%
hold off
