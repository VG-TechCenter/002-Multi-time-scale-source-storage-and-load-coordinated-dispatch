function Mbd = getMbdTR(A,bus)
%%输入所在节点向量，bus参数
nd = max(size(A));
nb = size(bus,1);

Mbd = zeros(nb,nd);

for ii = 1:nb
    for jj = 1:nd
        if ii == A(jj)
            Mbd(ii,jj) = 1;
        end
    end
end
