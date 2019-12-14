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


%% Twiddle factors index construction
%clear;
M_test = 16; L_test = 4;
layers = log2(L_test);
% max_index = (2^layers-1);
max_index = (M_test/2-1);
initial_index = 0:1:max_index;
index_table = zeros(layers,max_index);

for l = 1:1:layers
    for i = initial_index
        index_table(layers+1-l,i+1) = bitshift(bitshift(i,-(l-1)),(l-1));
    end
end

