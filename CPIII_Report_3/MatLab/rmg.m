function [rl,cp] = rmg(input)

N = length(input);

rl = rand(N);
cp = rl + 1i*rand(N);
    
v_rl = orth(rl);
v_cp = orth(cp);

rl = zeros(N);
cp = zeros(N);
for i=1:N
    rl = rl + v_rl(:,i)*input(i)*(v_rl(:,i)');
    cp = cp + v_cp(:,i)*input(i)*(v_cp(:,i)');
end
rl = 0.5*(rl+rl');
cp = 0.5*(cp+cp');
end