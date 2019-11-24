x = [1 2 3 4 3 2 1 2 3 4 3 2 1 2 3 4 3 2 1 2 3 4 3 2 ];

%x = [1 2 3 4 1 2 3 4 1 2 3 4];

N = length(x);
M = 4;
L = 2;

% Quantidade de janelas calculadas
Q = (N-M)/L+1;
Xn_k = zeros(Q,M);


for q = 1:1:Q
    n_l = (q-1)*L+1;
    Xn_k(q,:)=fft(x(n_l:n_l+M-1));
end