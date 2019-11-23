% Versão do Singleton (1969) utilizando laço for

%x = [1 5 3 -1 -2 -3 1 -9 1 1 1 1 1 1 1 1 0 0 0 0 2 2 2 2 3 3 3 3 4 4 4 4];
%t = -2:4/512:1.99999;x = 2*sin(6*pi*t)+8*sin(15*pi*t)+4*sin(10*pi/4*t);
x = randn(1,128);
%x = [0 0 1 1 0 1 0 1];


N = length(x);                  % Comprimento da amostra
N_2 = N/2;                      % Metada do comprimento da amostra
camadas = log2(N);              % Número de camadas 
W = exp(-1j*2*pi*[0:N_2-1]/N);   % W_N com valores de e^j2pik/N

% Frequências digitais normalizadas para plot 
f = [0:1:(N-1)]/N;

X_0 = fft(x);
close
subplot(3,2,1)
stem(f,abs(X_0)/N,'k')
title('Amplitude da FFT do MATLAB')
xlabel('\omega (x\pi rads/amostra)')
hold on
subplot(3,2,2)
stem(f,phase(X_0)/N,'k')
title('Fase da FFT do MATLAB')
xlabel('\omega (x\pi rads/amostra)')


% Vetor de índices em notação binária reversa
index = bi2de(de2bi(0:N-1,log2(N),'left-msb'))+1;

X_5 = zeros(2,N);
X_5(1,:) = x(index);                 % Entrada para algorítimo do Singleton
%X_temp = zeros(1,N);            % Armazena resultados parciais das camada

% Máscara de bits para cálculo do índice de W
w_index_mask = N-1;                     % Set os log2(N) primeiros bits
w_index_mask = w_index_mask*2^camadas;  % Desloca de log2(N)-1 bits

all_index = zeros(camadas*N_2,12);

for c= 1:1:camadas
   for e = 1:2:N
        % Cálculos dos índices do elemento e dos W_N's
        i_x = bitshift(e,-1)+1;     % i_x-ésimo par de elementos
        i_x_in = ~bitand(c,1)+1;    % Linha da matriz que é entrada
        i_x_out = bitand(c,1)+1;    % Linha da matriz que é saída
        i_w = bitshift(bitand(e,w_index_mask),-1)+1;    % ìndice do W
        
        % Cálculo das amostras pares e ímpares
        X_5(i_x_out,i_x)       = X_5(i_x_in,e)+W(i_w).*X_5(i_x_in,e+1);
        X_5(i_x_out,i_x+N_2)   = X_5(i_x_in,e)-W(i_w).*X_5(i_x_in,e+1);
        
        %N_2^(c-1)+bitshift(e,-1)+1;
        all_index(N_2*(c-1)+bitshift(e,-1)+1,1:6) = [c e i_x i_x_in i_x_out i_w ];
   end
   % Apontando os valores de saída para serem novas entradas
   % X_5 = X_temp;
   % Desloca a máscara para direita
   w_index_mask = bitshift(w_index_mask,-1);
end

subplot(3,2,3)
stem(f,abs(X_5(i_x_out,:))/N,'b')
title('Amplitude da implementação com dois loops')
xlabel('\omega (x\pi rads/amostra)')

subplot(3,2,4)
stem(f,phase(X_5(i_x_out,:))/N,'b')
title('Fase da implementação com dois loops')
xlabel('\omega (x\pi rads/amostra)')


X_5 = zeros(2,N);
X_5(1,:) = x(index);                 % Entrada para algorítimo do Singleton

% Máscara de bits para cálculo do índice de W
l_mask = N_2-1;                         % Set os log2(N) primeiros bits
h_mask = bitshift(l_mask,camadas-1);    % Desloca de log2(N)-1 bits
i_x_mask = N_2-1;                       % Máscara dos índices das entradas
for g = 0:1:(camadas*N_2-1)
        g_l = bitand(g,l_mask);            % Parte alta de g
        g_h = bitshift(g,-(camadas-1));            % Parte baixa de g
        e_o = bitset(bitshift(g_l,1),1);            % Elemento ímpar
        e_e = e_o+1;          % Elemento par
        i_x = g_l+1;     % i_x-ésimo par de elementos
        i_x_in = (bitand(g,N_2)>0)+1;    % Linha da matriz que é entrada
        i_x_out = (~bitand(g,N_2)>0)+1;    % Linha da matriz que é saída
        % Índice do W: (g&(N/2-1))&((N/2-1)<<(g>>log(N/2)))
        % A parte alta de g indica quais bits da parte baixa usar
        i_w = (bitand(g_l,bitshift(h_mask,-bitshift(g,-(camadas-1))))+1); 
        % Cálculo das amostras pares e ímpares
        X_5(i_x_out,i_x)       = X_5(i_x_in,e_o)+W(i_w).*X_5(i_x_in,e_e);
        X_5(i_x_out,i_x+N_2)   = X_5(i_x_in,e_o)-W(i_w).*X_5(i_x_in,e_e);
        
        all_index(g+1,7:12) = [g_h+1 e_o i_x i_x_in i_x_out i_w ];
end

subplot(3,2,5)
stem(f,abs(X_5(i_x_out,:))/N,'r')
title('Amplitude da DFT implementada com 1 loop')
xlabel('\omega (x\pi rads/amostra)')
subplot(3,2,6)
stem(f,phase(X_5(i_x_out,:))/N,'r')
title('Fase da DFT implementada com 1 loop')
xlabel('\omega (x\pi rads/amostra)')

hold off