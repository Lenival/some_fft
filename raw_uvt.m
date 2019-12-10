function D = raw_uvt(x,N,M,L,Q,q_D,X_start,bins)

d = zeros(1,N);
d(1:M) = x(1:M);
d(M+1:N)=x(M+1:N)-x(1:N-M);

D = zeros(Q+q_D,M);

for n_i = 1:1:q_D
    n_x = X_start + (n_i-1)*L;      % Index n of x(n) window 
    m_i = (n_x-L);                  % First index to calculate D_n(k)
    m_i = (m_i>0)*m_i+1;            % Start in 1 for not positive index
    
    for k = bins         % q=0 is calculated before
        for m = m_i:1:n_x     % Matlab index related to n-m HDFT index
            D(n_i,k+1)=D(n_i,k+1)+d(m)*exp(1j*2*pi*(((n_x-m)-L+1)*k)/M);
        end
    end
end

for n_i = q_D+1:1:(Q+q_D)
    
    n_x = M+(n_i-q_D-1)*L;
    m_i = (n_x-L)+1;
    
    for k = bins         % q=0 is calculated before
        for m = m_i:1:n_x     % Matlab index related to n-m HDFT index
            D(n_i,k+1)=D(n_i,k+1)+d(m)*exp(1j*2*pi*(((n_x-m)-L+1)*k)/M); 
        end
    end
end
