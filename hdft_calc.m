function Xn_k = hdft_calc(x,M,L)

N = length(x);


% Number of windows demanded
Q = floor((N-M)/(L))+1;
Xn_k = zeros(Q,M);

d = zeros(1,N);
d(1:M) = x(1:M);
d(M+1:N)=x(M+1:N)-x(1:N-M);

D = zeros(Q,M);

for n_i = 1:1:Q
    %n_x = X_start + (n_i-1)*M_h;      % Index n of x(n) window
    n_x = M+(n_i-1)*L;
    m_i = (n_x-L)+1;
    %m_f = n;
    for k = 0:1:M-1         % q=0 is calculated before
        for m = m_i:1:n_x     % Matlab index related to n-m HDFT index
            D(n_i,k+1)=D(n_i,k+1)+d(m)*exp(1j*2*pi*(((n_x-m)-L+1)*k)/M); 
        end
    end
end


%% HDFT spectrogram
Xn_k(1,1:M)=fft(x(1:1:M));


% Processing the remaining windows
for n_i = 2:1:Q        % n=0 is calculated before
    for k = 0:1:M-1
        %n_l = (n-1)*L+1;
        Xn_k(n_i,k+1)=exp(1j*(2*pi/M)*L*k)*(Xn_k(n_i-1,k+1)+D(n_i,k+1));
    end
end

Xn_k = fftshift(Xn_k);
