function [upwone,upwtwo,Perssone,Persstwo,Pidrbone,Pidrbtwo]=main2(k,pload,ugone,ugtwo,Psone,Pstwo,Ppdrone,Ppdrtwo,Pidraone,Pidratwo)
pw10=[61 60.5 70 63 64 102 131 165 167 182 153 179 182 121 99 293 272 245 210 182 145 105 96 85];%风电正调峰
pw20=[92 150 125 162 148 132 141 96 43 19 72 96 49 132 141 180 165 192 148 196 143 81 21 45];%风电反调峰
T=4*4;%时间
caseName = case30;
p30=sum(caseName.bus(:,3));%当前节点功率之和
pbl=pload./p30;
pload1=repmat(pbl,30,1).*repmat(caseName.bus(:,3),1,T);
% P_D=0.85.*pload1;%固定负荷
% for i=4*(k-1)+1:4*k
%     for j=4*(i-4*(k-1)-1)+1:4*(i-4*(k-1))
%         pw12(:,j)=pw1(:,i);
%         pw22(:,j)=pw2(:,i);
%         pload1(:,j)=pload1(:,i);
%     end
% end
for i=1:24
    for j=4*(i-1)+1:4*i
    ugone1(:,j)=ugone(:,i);
    ugtwo1(:,j)=ugtwo(:,i);
    Psone1(:,j)=Psone(:,i);
    Pstwo1(:,j)=Pstwo(:,i);
    Ppdrone1(:,j)=Ppdrone(:,i);
    Ppdrtwo1(:,j)=Ppdrtwo(:,i);
    Pidraone1(:,j)=Pidraone(:,i);
    Pidratwo1(:,j)=Pidratwo(:,i);
    pw1(:,j)=pw10(:,i);
    pw2(:,j)=pw20(:,i);
    end
end
% ugone1=ugone;ugtwo1=ugtwo;Psone1=Psone;Pstwo1=Pstwo;Ppdrone1=Ppdrone;Ppdrtwo1=Ppdrtwo;Pidraone1=Pidraone;Pidratwo1=Pidratwo;

ugone2=zeros(6,T);ugtwo2=zeros(6,T);
Psone2=zeros(1,T);
Ppdrone2=zeros(30,T);
Pidraone2=zeros(30,T);
Pstwo2=zeros(1,T);
Ppdrtwo2=zeros(30,T);
Pidratwo2=zeros(30,T);
% for i=4*(k-1)+1:4*k
%     for j=4*(i-4*(k-1)-1)+1:4*(i-4*(k-1))
        ugone2(:,1:T)=ugone1(:,16*(k-1)+1:16*k);
        Psone2(:,1:T)=Psone1(:,16*(k-1)+1:16*k);
        Ppdrone2(:,1:T)=Ppdrone1(:,16*(k-1)+1:16*k);
        Pidraone2(:,1:T)=Pidraone1(:,16*(k-1)+1:16*k);
        ugtwo2(:,1:T)=ugtwo1(:,16*(k-1)+1:16*k);
        Pstwo2(:,1:T)=Pstwo1(:,16*(k-1)+1:16*k);
        Ppdrtwo2(:,1:T)=Ppdrtwo1(:,16*(k-1)+1:16*k);
        Pidratwo2(:,1:T)=Pidratwo1(:,16*(k-1)+1:16*k);
        pw12(:,1:T)=pw1(:,16*(k-1)+1:16*k);
        pw22(:,1:T)=pw2(:,16*(k-1)+1:16*k);
%     end
% end
%参数
Ns=2;%场景数为2
ps=[0.3 0.7];%每种场景概率
Nerss=1;%储能电站个数
Nbt=1;%电化学储能电站个数
PIbt=2.3;%电化学储能单位时间折旧成本
NG=6;%常规机组数量
NDG=1;%DG数量
kc=2.5;%弃风惩罚
socmin=0.1;
socmax=0.9;
Pech=-50;%充电功率限制
Pedis=50;%放电功率限制
Se=200;%储能容量限制
Se0=100;
Pwmin=-50;%蓄水电站最小功率
Pwmax=50;%蓄水电站最大功率
Vs0=100;
Vpmin=0;%蓄水容量最小值
Vpmax=200;%蓄水容量最大值
Ca=100;Cb=125;Cc=150;Cd=150;%IDR补偿成本系数
mp=[1 150 50 0.0375 20 372.5 72 2;%常规机组参数
    2 60  20  0.175 17.5 352.3 48 2;
    3 60  15  0.625 10 316.5  30 2;
    4 50  10  0.0834 32.5 329.2 30 2;
    5 40  10  0.25  30  276.4 18 2;
    6 45  12  0.25 30 232.2  24 2];
%pload=[390 382 354 340 328 409 460 460 472 511 465 458 442 456 457 471 475 503 512 469 445 438 410 362];%负荷

idg=[1 2 5 8 11 13];
iw=2;
ie=2;
is=8;
%变量
%场景1变量
thetaone=sdpvar(30,T,'full');
Pgone=sdpvar(NG,T,'full');
Pwone=sdpvar(1,T,'full');
upwone=binvar(1,T,'full');
Perssone=sdpvar(1,T,'full');
uerssone=binvar(1,T,'full');
Pdgone=sdpvar(1,T,'full');
Pidrbone=sdpvar(30,T,'full');
Pidrcone=sdpvar(30,T,'full');
plossone=sdpvar(30,T,'full');
%场景2变量
thetatwo=sdpvar(30,T,'full');
Pgtwo=sdpvar(NG,T,'full');
Pwtwo=sdpvar(1,T,'full');
upwtwo=binvar(1,T,'full');
Pstwo=sdpvar(1,T,'full');
Persstwo=sdpvar(1,T,'full');
uersstwo=binvar(1,T,'full');
Pdgtwo=sdpvar(1,T,'full');
Pidrbtwo=sdpvar(30,T,'full');
Pidrctwo=sdpvar(30,T,'full');
plosstwo=sdpvar(30,T,'full');
%约束
con=[];
con_gen=getConsGen1(Pgone,Pgtwo,mp,T,ugone2,ugtwo2);%传统机组约束--出力，爬坡，时间
con=[con,con_gen];
con=[con,0<=Pwone<=upwone.*pw12,0<=Pwtwo<=upwtwo.*pw22];%风电约束
% % 线路传输约束
%   [cons_pf, pf] = getConsPF(caseName, x_theta, Horizon);
%   cons = [cons, cons_pf];
%   cons = [cons, getConsAgl(x_theta)];
%潮流约束
baseMVA = 100;
    bus = caseName.bus;
    brch = caseName.branch;
    f = brch(:,1);
    t = brch(:,2);
    x = brch(:,4);
    nbus = size(bus,1);
    nbrch = size(brch,1);
    Cft = zeros(nbrch,nbus);
    for ii=1:nbrch
        Cft(ii,f(ii)) = 1;
        Cft(ii, t(ii)) = -1;
    end
    Bf = zeros(nbrch,nbus);
    for ii=1:nbrch
        Bf(ii,f(ii)) = 1./x(ii);
        Bf(ii,t(ii)) = -1./x(ii);
    end
    Bbus = Cft.' * Bf;
    %% 生成Pbus
    Mbg = getMbgMatrix(idg',bus);
    Mbe = getMbdMatrix(ie,bus);
    Mbw = getMbdMatrix(iw,bus);
    Mbs = getMbdMatrix(is,bus);
    Pbusone = [];
    for ii = 1:T
        Pbusone = [Pbusone, Mbg * Pgone(:,ii) + Mbs * Psone(1,ii) + Mbe * Perssone(1,ii) + ...
            Mbw * Pwone(1,ii) - pload1(:,ii) - Ppdrone2(:,ii) - Pidraone2(:,ii) - Pidrbone(:,ii) - Pidrcone(:,ii) - plossone(:,ii)];
    end
    %% 潮流平衡约束，弧度制，有名值
     for ii = 1:T
         con = [con, Pbusone(:,ii) == Bbus * thetaone(:,ii).*baseMVA];
     end
%场景2
    Pbustwo = [];
    for ii = 1:T
        Pbustwo = [Pbustwo, Mbg * Pgtwo(:,ii) + Mbs * Pstwo(1,ii) + Mbe * Persstwo(1,ii) + ...
            Mbw * Pwtwo(1,ii) - pload1(:,ii) - Ppdrtwo2(:,ii) - Pidratwo2(:,ii) - Pidrbtwo(:,ii)- Pidrctwo(:,ii) - plosstwo(:,ii)];
    end
    %% 潮流平衡约束，弧度制，有名值
     for ii = 1:T
         con = [con, Pbustwo(:,ii) == Bbus * thetatwo(:,ii).*baseMVA];
     end
    %潮流约束
for j=1:T
d_thetaone(j,:) = thetaone(f,j) - thetaone(t,j);
pfone = d_thetaone(j,:)' ./ x .* baseMVA; 
con = [con, -500 <=  pfone <=500];
end
%场景2约束
for j=1:T
d_thetatwo(j,:) = thetatwo(f,j) - thetatwo(t,j);
pftwo = d_thetatwo(j,:)' ./ x .* baseMVA; 
con = [con, -500 <=  pftwo <=500];
end
con = [con,-pi <= thetaone <= pi,-pi <= thetatwo <= pi];
%DR约束

con=[con,-0.05.*pload1<=Pidraone2+Pidrbone+Pidrcone<=0.05.*pload1,-0.05.*pload1<=Pidrbone<=0.05.*pload1,-0.05.*pload1<=Pidrcone<=0.05.*pload1,sum(Pidraone2+Pidrbone+Pidrcone,2)==0];
con=[con,-0.05.*pload1<=Pidratwo2+Pidrbtwo+Pidrctwo<=0.05.*pload1,-0.05.*pload1<=Pidrbtwo<=0.05.*pload1,-0.05.*pload1<=Pidrctwo<=0.05.*pload1,sum(Pidratwo2+Pidrbtwo+Pidrctwo,2)==0];

%电储能约束
con=[con,Pech.*uerssone<=Perssone<=Pedis.*uerssone];
for i=1:T
    con=[con,Se*socmin<=Se0+sum(Perssone(:,1:i),2)<=Se*socmax];
end
%场景2
con=[con,Pech.*uersstwo<=Persstwo<=Pedis.*uersstwo];
for i=1:T
    con=[con,Se*socmin<=Se0+sum(Persstwo(:,1:i),2)<=Se*socmax];
end
%抽水蓄能
% con=[con,Pwmin<=Psone<=Pwmax];
% for i=1:T
%     con=[con,Vpmin<=Vs0+sum(Psone(:,1:i),2)<=Vpmax];
% end
%场景2
% con=[con,Pwmin<=Pstwo<=Pwmax];
% for i=1:T
%     con=[con,Vpmin<=Vs0+sum(Pstwo(:,1:i),2)<=Vpmax];
% end
%目标函数
%常规机组运行成本a2*p+bp+c分段线性化
tn=T;
Pmax=[mp(1,2),mp(2,2),mp(3,2),mp(4,2),mp(5,2),mp(6,2)];
Pmin=[mp(1,3),mp(2,3),mp(3,3),mp(4,3),mp(5,3),mp(6,3)];
gn=5;%分段函数线性化，下同
x_pf=sdpvar(6, tn,'full');
gw1=sdpvar(gn+1,tn,'full');
gw2=sdpvar(gn+1,tn,'full');
gw3=sdpvar(gn+1,tn,'full');
gw4=sdpvar(gn+1,tn,'full');
gw5=sdpvar(gn+1,tn,'full');
gw6=sdpvar(gn+1,tn,'full');
gz1=binvar(gn, tn,'full');gz2=binvar(gn, tn,'full');gz3=binvar(gn, tn,'full');gz4=binvar(gn, tn,'full');gz5=binvar(gn, tn,'full');gz6=binvar(gn, tn,'full');
gl1=(Pmax-Pmin)./gn;
gl2=zeros(6,gn+1);
for i=1:6
gl2(i,:)=Pmin(i):gl1(i):Pmax(i);
end
con = [con, x_pf(1,:)==gl2(1,:).^2*gw1];
con = [con, x_pf(2,:)==gl2(2,:).^2*gw2];
con = [con, x_pf(3,:)==gl2(3,:).^2*gw3];
con = [con, x_pf(4,:)==gl2(4,:).^2*gw4];
con = [con, x_pf(5,:)==gl2(5,:).^2*gw5];
con = [con, x_pf(6,:)==gl2(6,:).^2*gw6];
con = [con, gw1(1,:)<=gz1(1,:)];
for i=2:gn
    con = [con, gw1(i,:)<=gz1(i-1,:)+gz1(i,:)];
end
con = [con, gw1(gn+1,:)<=gz1(gn,:)];
con = [con, sum(gz1)==ones(1,tn)];
con = [con, gw2(1,:)<=gz2(1,:)];
for i=2:gn
    con = [con, gw2(i,:)<=gz2(i-1,:)+gz2(i,:)];
end
con = [con, gw2(gn+1,:)<=gz2(gn,:)];
con = [con, sum(gz2)==ones(1,tn)];
con = [con, gw3(1,:)<=gz3(1,:)];
for i=2:gn
    con = [con, gw3(i,:)<=gz3(i-1,:)+gz3(i,:)];
end
con = [con, gw3(gn+1,:)<=gz3(gn,:)];
con = [con, sum(gz3)==ones(1,tn)];
con = [con, gw4(1,:)<=gz4(1,:)];
for i=2:gn
    con = [con, gw4(i,:)<=gz4(i-1,:)+gz4(i,:)];
end
con = [con, gw4(gn+1,:)<=gz4(gn,:)];
con = [con, sum(gz4)==ones(1,tn)];
con = [con, gw5(1,:)<=gz5(1,:)];
for i=2:gn
    con = [con, gw5(i,:)<=gz5(i-1,:)+gz5(i,:)];
end
con = [con, gw5(gn+1,:)<=gz5(gn,:)];
con = [con, sum(gz5)==ones(1,tn)];
con = [con, gw6(1,:)<=gz6(1,:)];
for i=2:gn
    con = [con, gw6(i,:)<=gz6(i-1,:)+gz6(i,:)];
end
con = [con, gw6(gn+1,:)<=gz6(gn,:)];
con = [con, sum(gz6)==ones(1,tn)];
con = [con, Pgone(1,:)==gl2(1,:)*gw1];
con = [con, Pgone(2,:)==gl2(2,:)*gw2];
con = [con, Pgone(3,:)==gl2(3,:)*gw3];
con = [con, Pgone(4,:)==gl2(4,:)*gw4];
con = [con, Pgone(5,:)==gl2(5,:)*gw5];
con = [con, Pgone(6,:)==gl2(6,:)*gw6];
con = [con, gw1>=0];
con = [con, gw2>=0];
con = [con, gw3>=0];
con = [con, gw4>=0];
con = [con, gw5>=0];
con = [con, gw6>=0];
%场景2线性化
x_pfz=sdpvar(6, tn,'full');
gw1z=sdpvar(gn+1,tn,'full');
gw2z=sdpvar(gn+1,tn,'full');
gw3z=sdpvar(gn+1,tn,'full');
gw4z=sdpvar(gn+1,tn,'full');
gw5z=sdpvar(gn+1,tn,'full');
gw6z=sdpvar(gn+1,tn,'full');
gz1z=binvar(gn, tn,'full');gz2z=binvar(gn, tn,'full');gz3z=binvar(gn, tn,'full');gz4z=binvar(gn, tn,'full');gz5z=binvar(gn, tn,'full');gz6z=binvar(gn, tn,'full');
gl1z=(Pmax-Pmin)./gn;
gl2z=zeros(6,gn+1);
for i=1:6
gl2z(i,:)=Pmin(i):gl1z(i):Pmax(i);
end
con = [con, x_pfz(1,:)==gl2z(1,:).^2*gw1z];
con = [con, x_pfz(2,:)==gl2z(2,:).^2*gw2z];
con = [con, x_pfz(3,:)==gl2z(3,:).^2*gw3z];
con = [con, x_pfz(4,:)==gl2z(4,:).^2*gw4z];
con = [con, x_pfz(5,:)==gl2z(5,:).^2*gw5z];
con = [con, x_pfz(6,:)==gl2z(6,:).^2*gw6z];
con = [con, gw1z(1,:)<=gz1z(1,:)];
for i=2:gn
    con = [con, gw1z(i,:)<=gz1z(i-1,:)+gz1z(i,:)];
end
con = [con, gw1z(gn+1,:)<=gz1z(gn,:)];
con = [con, sum(gz1z)==ones(1,tn)];
con = [con, gw2z(1,:)<=gz2z(1,:)];
for i=2:gn
    con = [con, gw2z(i,:)<=gz2z(i-1,:)+gz2z(i,:)];
end
con = [con, gw2z(gn+1,:)<=gz2z(gn,:)];
con = [con, sum(gz2z)==ones(1,tn)];
con = [con, gw3z(1,:)<=gz3z(1,:)];
for i=2:gn
    con = [con, gw3z(i,:)<=gz3z(i-1,:)+gz3z(i,:)];
end
con = [con, gw3z(gn+1,:)<=gz3z(gn,:)];
con = [con, sum(gz3z)==ones(1,tn)];
con = [con, gw4z(1,:)<=gz4z(1,:)];
for i=2:gn
    con = [con, gw4z(i,:)<=gz4z(i-1,:)+gz4z(i,:)];
end
con = [con, gw4z(gn+1,:)<=gz4z(gn,:)];
con = [con, sum(gz4z)==ones(1,tn)];
con = [con, gw5z(1,:)<=gz5z(1,:)];
for i=2:gn
    con = [con, gw5z(i,:)<=gz5z(i-1,:)+gz5z(i,:)];
end
con = [con, gw5z(gn+1,:)<=gz5z(gn,:)];
con = [con, sum(gz5z)==ones(1,tn)];
con = [con, gw6z(1,:)<=gz6z(1,:)];
for i=2:gn
    con = [con, gw6z(i,:)<=gz6z(i-1,:)+gz6z(i,:)];
end
con = [con, gw6z(gn+1,:)<=gz6z(gn,:)];
con = [con, sum(gz6z)==ones(1,tn)];
con = [con, Pgtwo(1,:)==gl2z(1,:)*gw1z];
con = [con, Pgtwo(2,:)==gl2z(2,:)*gw2z];
con = [con, Pgtwo(3,:)==gl2z(3,:)*gw3z];
con = [con, Pgtwo(4,:)==gl2z(4,:)*gw4z];
con = [con, Pgtwo(5,:)==gl2z(5,:)*gw5z];
con = [con, Pgtwo(6,:)==gl2z(6,:)*gw6z];
con = [con, gw1z>=0];
con = [con, gw2z>=0];
con = [con, gw3z>=0];
con = [con, gw4z>=0];
con = [con, gw5z>=0];
con = [con, gw6z>=0];    
for gt=1:T
    for gi=1:6
        fpg1(gi,gt)=mp(gi,4).*x_pf(gi,gt)+mp(gi,5).*Pgone(gi,gt)+mp(gi,6);
        fpg2(gi,gt)=mp(gi,4).*x_pfz(gi,gt)+mp(gi,5).*Pgtwo(gi,gt)+mp(gi,6);
    end
end
 yyone=binvar(6,T-1,'full');
 yytwo=binvar(6,T-1,'full');
 for i=1:T-1
     for j=1:6
     con=[con,yyone(j,i)<=ugone2(j,i+1)];
     con=[con,yyone(j,i)<=1-ugone2(j,i)];
     con=[con,yyone(j,i)>=ugone2(j,i+1)-ugone2(j,i)];
     con=[con,yytwo(j,i)<=ugtwo2(j,i+1)];
     con=[con,yytwo(j,i)<=1-ugtwo2(j,i)];
     con=[con,yytwo(j,i)>=ugtwo2(j,i+1)-ugtwo2(j,i)];
     end
 end
  yygone=binvar(1,T-1,'full');
  yygtwo=binvar(1,T-1,'full');
 for i=1:T-1
     con=[con,yygone(1,i)<=upwone(1,i+1)];
     con=[con,yygone(1,i)<=1-upwone(1,i)];
     con=[con,yygone(1,i)>=upwone(1,i+1)-upwone(1,i)];
     con=[con,yygtwo(1,i)<=upwtwo(1,i+1)];
     con=[con,yygtwo(1,i)<=1-upwtwo(1,i)];
     con=[con,yygtwo(1,i)>=upwtwo(1,i+1)-upwtwo(1,i)];
 end
 con=[con,0<=plossone<=pload1,0<=plosstwo<=pload1];
fg=sum(sum(ps(1).*fpg1+ps(2).*fpg2))+0.9*sum(sum(ps(1).*yyone+ps(2).*yytwo));
fess=sum(ps(1)*(0.53*abs(Perssone)+0.12*abs(Perssone))+ps(2)*(0.53*abs(Persstwo)+0.12*abs(Persstwo)))+sum(0.012.*(uerssone+uersstwo));
fdg=ps(1)*sum(0.02.*Pwone)+ps(2)*(sum(0.02.*Pwtwo)+0.9*sum(ps(1).*yygone+ps(2).*yygtwo))+ps(1)*5.2*sum(pw12-Pwone)+ps(2)*5.2*sum(pw22-Pwtwo);
fload=ps(1)*(sum(sum(125.*abs(Pidrbone))+sum(150.*abs(Pidrcone)))+15.9.*sum(sum(plossone)))+ps(2)*(sum(sum(100.*abs(Pidratwo2)))+sum(150.*abs(Pidrctwo))+15.9.*sum(sum(plosstwo)));
f1=fg+fess+fdg+fload;
ops=sdpsettings('solver','gurobi');
result=optimize(con,f1,ops);
upwone=value(upwone);
upwtwo=value(upwtwo);
Perssone=value(Perssone);
Persstwo=value(Persstwo);
Pidrbone=value(Pidrbone);
Pidrbtwo=value(Pidrbtwo);



