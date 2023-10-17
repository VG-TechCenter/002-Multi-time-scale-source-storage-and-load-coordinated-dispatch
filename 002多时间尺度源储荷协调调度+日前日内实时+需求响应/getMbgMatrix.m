function Mbg = getMbgMatrix(gen,bus)
%% ‰»Îgen£¨bus≤Œ ˝
ng = size(gen,1);
nb = size(bus,1);

Mbg = zeros(nb,ng);

for ii = 1:nb
    for jj = 1:ng
        if ii == gen(jj,1)
            Mbg(ii,jj) = 1;
        end
    end
end