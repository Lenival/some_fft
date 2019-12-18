function D = uvt_dit_r2(x,N,M,L,Q,q_D,X_start)

% How many bins will be calc
bins = 0:1:M-1;

% How many zeroes before start
n_pad = L-(X_start);

% Update vector
d = zeros(1,N+n_pad);
d(n_pad+1:M+n_pad) = x(1:M);
d(M+n_pad+1:N+n_pad)=x(M+1:N)-x(1:N-M);

% Twiddle factors for decimation of UVT calc
W_M_mk = exp(-1j*(2*pi/M)*(0:M/2-1));

v = log2(L);

max_index = (M/2-1);

initial_index = 0:1:max_index;
W_i = zeros(v,max_index);
disp('------------------------------------------')
disp('               --------                   ')

for l = 1:1:v

    for i = initial_index
        W_i(v+1-l,i+1) = bitshift(bitshift(i,-(l-1)),(l-1));
    end
end

% Index examples for M=16 and L=4
%De1 0 4; 1 5;  2 6;  3 7; 8 12; 9 13; 10 14; 11 15
%We1 +0-; +2-;  +4-;  +6-;  +0-;  +2-;   +4-;   +6-     
%De2 0 8; 1 9; 2 10; 3 11; 4 12; 5 13;  6 14;  7 15
%We2 +0-; +1-;  +2-;  +3-;  +4-;  +5-;   +6-;   +7-


D_i = zeros(v,max_index);
for l = 1:1:v
    D_index = 0;
    max_decimation = M/(2^l);
    for i = 1:1:max_decimation
        for c = 0:1:(2^l-1)
            D_index = D_index + 1;
            D_i(v+1-l,D_index) = max_decimation*c+(i-1);
        end
    end
end


% Updating Vector Transform (UVT)
D = zeros(Q+q_D,M);

bit_rindex = bi2de(de2bi(L-1:-1:0,log2(L),'left-msb'))+1;

for n_i = 1:1:Q+q_D
    for m = 1:1:L
        for l = 1:1:(M/L)
            %disp([m l (m-1)*(M/L)+l (n_i-1)*L+m])
%             D(n_i,(m-1)*(M/L)+l) = d((n_i-1)*L+m);
            D(n_i,(m-1)*(M/L)+l) = d((n_i-1)*L+bit_rindex(m));
        end
    end
end

% for c = 1:1:v
%     for e = 1:2:M
%         D() = ;
%     end
% end


% UVT calc
for n_i = 1:1:(Q+q_D)
    n_x = M+n_pad+(n_i-q_D-1)*L;        % Window index
    
    for c = 1:1:v
        for e = 1:1:M/2
%             disp([n_i c e W_i(c,e) D_i(c,2*e-1) D_i(c,2*e) ])
%             disp([n_x D(n_i,D_i(c,2*e-1)+1) D(n_i,D_i(c,2*e)+1)])
            
            
            D_temp = W_M_mk(W_i(c,e)+1)*D(n_i,D_i(c,2*e-1)+1);
            D(n_i,D_i(c,2*e-1)+1) = D(n_i,D_i(c,2*e)+1)+D_temp;
            D(n_i,D_i(c,2*e)+1)   = D(n_i,D_i(c,2*e)+1)-D_temp;
        end
    end
%     for k = bins
%         for m = 0:1:L-1     % Matlab index related to n-m HDFT index
%             D(n_i,k+1)= D(n_i,k+1)+d(n_x-m)*W_M_mk(m+1,k+1);
%         end
%     end
end

