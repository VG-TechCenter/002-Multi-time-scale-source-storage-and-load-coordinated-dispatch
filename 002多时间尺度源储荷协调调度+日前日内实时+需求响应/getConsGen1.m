function cons = getConsGen1(Pgone,Pgtwo,mp,T,ugone,ugtwo)
%% ��ȡ����Լ��
Pgmin=mp(:,3);
Pgmax=mp(:,2);
rud=mp(:,7);
On_min=2;
Off_min=2;
OnOff_history = [0 0];
cons = [];
% 1. ����������Լ��
cons = [cons, ugone.*repmat(Pgmin,1,T) <=Pgone<= ugone.*repmat(Pgmax,1,T)];
cons = [cons, ugtwo.*repmat(Pgmin,1,T) <=Pgtwo<= ugtwo.*repmat(Pgmax,1,T)];
% 2. ����Լ��
cons = [cons, abs([Pgone(:, 2:end),Pgone(:, 1)] - Pgone)<=repmat(rud,1,T)];
% 3. ����������Լ������С��ͣʱ��Լ��

for i=1:6
for t=1:T   
    cons=[cons,
%             ugone(i,t)*Pgmin(i) <= Pgone(t) <= ugone(i,t)*Pgmax(i), 
            consequtiveON([OnOff_history ugone(i,:)],On_min),
            consequtiveON(1-[OnOff_history ugone(i,:)],Off_min)
        ];
end
end
%����2Լ��
% 2. ����Լ��
cons = [cons, abs([Pgtwo(:, 2:end),Pgtwo(:, 1)] - Pgtwo)<=repmat(rud,1,T)];
% 3. ����������Լ������С��ͣʱ��Լ��
for i=1:6
for t=1:T   
    cons=[cons,
%             ugtwo(i,t)*Pgmin(i) <= Pgtwo(t) <= ugtwo(i,t)*Pgmax(i), 

            consequtiveON([OnOff_history ugtwo(i,:)],On_min),
            consequtiveON(1-[OnOff_history ugtwo(i,:)],Off_min)
        ];
end
end