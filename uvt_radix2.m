function D = uvt_radix2(x,N,M,L,Q,q_D,X_start,bins)

% How many zeroes before start
n_pad = L-(X_start);

% Update vector
d = zeros(1,N+n_pad);
d(n_pad+1:M+n_pad) = x(1:M);
d(M+n_pad+1:N+n_pad)=x(M+1:N)-x(1:N-M);

% Updating Vector Transform (UVT)
D = zeros(Q+q_D,M);

% Twiddle factors for UVT
W_M_mk = zeros(L,M);
for k = bins
    for m = 0:1:L-1
        W_M_mk(m+1,k+1) = exp(1j*2*pi*((m-L+1)*k)/M);
    end
end

% UVT calc
for n_i = 1:1:(Q+q_D)
    n_x = M+n_pad+(n_i-q_D-1)*L;        % Window index
    for k = bins
        for m = 0:1:L-1     % Matlab index related to n-m HDFT index
            D(n_i,k+1)= D(n_i,k+1)+d(n_x-m)*W_M_mk(m+1,k+1);
        end
    end
end

