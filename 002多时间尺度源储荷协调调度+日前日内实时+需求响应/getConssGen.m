function cons = getConssGen(x_P_sg, Psmax, Psmin, Horizon)
cons = [];
% 1. 机组上下限约束
cons = [cons, repmat(Psmin,1,Horizon) <=x_P_sg <= repmat(Psmax,1,Horizon)];

