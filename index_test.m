%% HDFT start point and number of hypothetical windows
I_x = zeros(20,20);
Ms = 1:1:20;
s = 1:1:20;
for m = Ms
    for l = s
        I_x(m,l) = (m>l)*mod(m-1,m-l);
    end
end

Q_D = zeros(20,20);
for m = Ms
    for l = s
        Q_D(m,l) = (m>l)*floor((m-1)/(m-l));
    end
end


%% Twiddle factors and UVT decimation index construction
%clear;
M_test = 12; L_test = 4;
layers = log2(L_test);
max_index = (M_test/2-1);

initial_index = 0:1:max_index;
W_index_table = zeros(layers,max_index);
disp('------------------------------------------')

disp('------------------------------------------')
disp('               --------                   ')
disp('               --------                   ')
disp('')

for l = 1:1:layers

    for i = initial_index
disp('------------------------------------------')
        disp([l i ])
        W_index_table(layers+1-l,i+1) = bitshift(bitshift(i,-(l-1)),(l-1));
    end
end

% Index examples for M=16 and L=4
%De1 0 4; 1 5;  2 6;  3 7; 8 12; 9 13; 10 14; 11 15
%We1 +0-; +2-;  +4-;  +6-;  +0-;  +2-;   +4-;   +6-     
%De2 0 8; 1 9; 2 10; 3 11; 4 12; 5 13;  6 14;  7 15
%We2 +0-; +1-;  +2-;  +3-;  +4-;  +5-;   +6-;   +7-

disp('------------------------------------------')
disp('------------------------------------------')
disp('               --------                   ')
disp('               --------                   ')
disp('')

D_index_table = zeros(layers,max_index);
for l = 1:1:layers
    D_index = 0;
    max_decimation = M_test/(2^l);
    for i = 1:1:max_decimation
disp('------------------------------------------')
        disp([l i ])
        for c = 0:1:(2^l-1)
            D_index = D_index + 1;
            D_index_table(layers+1-l,D_index) = max_decimation*c+(i-1);
        end
    end
end











