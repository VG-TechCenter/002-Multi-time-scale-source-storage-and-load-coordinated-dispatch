function Mbd = getMbdTR(A,bus)
%%�������ڽڵ�������bus����
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
