x = [1 5 3 -1 -2 -3 1 -9];


X_0 = fft(x);
subplot(2,1,1)
stem(abs(X_0))
hold on



N = length(x);
W = exp(j*2*pi*[0:(N/2)-1]/N);


X_2 =  x;



% Camada 1
% 0 4 ; 2 6 ; 1 5 ; 3 7
% 0 4 ; 0 4 ; 0 4 ; 0 4 --> 0 0 ; 0 0 ; 0 0 ; 0 0 
temp = X_2(0+1);
X_2(0+1) = temp+W(0+1).*X_2(4+1);
X_2(4+1) = temp-W(0+1).*X_2(4+1);
temp = X_2(2+1);
X_2(2+1) = temp+W(0+1).*X_2(6+1);
X_2(6+1) = temp-W(0+1).*X_2(6+1);
temp = X_2(1+1);
X_2(1+1) = temp+W(0+1).*X_2(5+1);
X_2(5+1) = temp-W(0+1).*X_2(5+1);
temp = X_2(3+1);
X_2(3+1) = temp+W(0+1).*X_2(7+1);
X_2(7+1) = temp-W(0+1).*X_2(7+1);

% Camada 2
% 0 2 ; 4 6 ; 1 3 ; 5 7
% 0 4 ; 2 6 ; 0 4 ; 2 6 -->  0 0 ; 2 2 ; 0 0 ; 2 2
temp = X_2(1);
X_2(1) = temp+W(0+1).*X_2(3);
X_2(3) = temp-W(0+1).*X_2(3);
temp = X_2(5);
X_2(5) = temp+W(2+1).*X_2(7);
X_2(7) = temp-W(2+1).*X_2(7);
temp = X_2(2);
X_2(2) = temp+W(0+1).*X_2(4);
X_2(4) = temp-W(0+1).*X_2(4);
temp = X_2(6);
X_2(6) = temp+W(2+1).*X_2(8);
X_2(8) = temp-W(2+1).*X_2(8);

% Camada 3
% 0 1 ; 4 5 ; 2 3 ; 6 7
% 0 4 ; 1 5 ; 2 6 ; 3 7 -->  0 0 ; 1 1 ; 2 2 ; 3 3
temp = X_2(1);
X_2(1) = temp+W(0+1).*X_2(2);
X_2(2) = temp-W(0+1).*X_2(2);
temp = X_2(5);
X_2(5) = temp+W(1+1).*X_2(6);
X_2(6) = temp-W(1+1).*X_2(6);
temp = X_2(3);
X_2(3) = temp+W(2+1).*X_2(4);
X_2(4) = temp-W(2+1).*X_2(4);
temp = X_2(7);
X_2(7) = temp+W(3+1).*X_2(8);
X_2(8) = temp-W(3+1).*X_2(8);

shuffleButerfly = bi2de(de2bi(0:N-1,log2(N),'left-msb'))+1;
X_2 = X_2(shuffleButerfly);

subplot(2,1,2)
stem(abs(X_2))


%
%for c = 1:1:log2(N)
%    c
%end
%
hold off
