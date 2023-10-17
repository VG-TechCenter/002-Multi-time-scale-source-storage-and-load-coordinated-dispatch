function cons = getConsEES(x_P_ch, x_P_dis, x_u_ch, x_u_dis, EESmax, EESmin, capmax, Horizon,theta)
cons = [];
soc0=0.5;
yitac=0.95;
yitad=0.95;
socmax=0.9;
socmin=0.2;
soc = sdpvar(1, Horizon);
% 1. 充放电功率约束，充电记为正数，放电记为负数
cons = [cons, x_u_ch.*repmat(EESmin, 1, Horizon) <= x_P_ch <= x_u_ch.*repmat(EESmax, 1, Horizon)];
cons = [cons, - x_u_dis.*repmat(EESmax, 1, Horizon) <=  x_P_dis <= - x_u_dis.*repmat(EESmin, 1, Horizon)];
% 2. 充放电次数约束
% cons = [cons, sum(x_u_ch, 2) <= num_chmax];
% cons = [cons, sum(x_u_dis, 2) <= num_dismax];
% 3. 不同时充放电约束
cons = [cons, 0 <= x_u_ch + x_u_dis <= 1];
% 4. SOC约束
cons = [cons,socmin <= soc0*(1-theta)+x_P_ch(:,1)*yitac/capmax+x_P_dis(:,1)/yitad/capmax <= socmax];
%cons = [cons,soc(1)== soc0+x_P_ch(:,1)];
soc(1)= soc0+x_P_ch(:,1);
for ii = 2:Horizon
    %cons = [cons, 0 <= sum(x_P_ch(:,1:ii) + x_P_dis(:,1:ii),2) <= capmax - capmin];
    cons = [cons,soc(ii)==soc(ii-1)*(1-theta)+x_P_ch(:,ii)*yitac/capmax+x_P_dis(:,ii)/yitad/capmax];
    cons = [cons,socmin <= soc(ii)<= socmax];
end
cons = [cons, sum(x_P_ch(:,1:Horizon) + x_P_dis(:,1:Horizon),2) == 0];
