% function [Pwone,Pwtwo]=mainsj()
%参数
T=288;%时间
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
pload=[390 382 354 340 328 409 460 460 472 511 465 458 442 456 457 471 475 503 512 469 445 438 410 362];%负荷
pw10=[61 60.5 70 63 64 102 131 165 167 182 153 179 182 121 99 293 272 245 210 182 145 105 96 85];%风电正调峰
pw20=[92 150 125 162 148 132 141 96 43 19 72 96 49 132 141 180 165 192 148 196 143 81 21 45];%风电反调峰
pload=[389.428017265482,390.430961345508,392.218173996095,392.165162193146,389.929894586542,389.948757029671,389.464992918487,389.134422632792,390.758448567705,394.084495522225,394.349496596856,393.595158025653,378.290646192536,380.276680216358,378.357797954302,382.196103662072,382.059123514806,382.225732627690,378.724206962968,378.830593678088,379.291252493939,378.396284010726,378.021491905031,378.148406421688,356.153642293038,355.480666026004,354.957447665065,352.530984106758,352.628580680505,352.152274667114,355.171420727130,354.438944295227,353.912180564402,354.483968302795,353.961102250217,355.581613656752,343.505887554493,340.711353112437,341.009983515865,338.631651291761,337.537403645821,339.605534068396,340.768710246172,342.748624181301,340.768013213536,339.839475055451,338.256706262778,339.027128583511,329.123973717114,331.095453620829,331.468623949010,329.691223118458,330.796543435869,329.504246193230,326.502227646863,327.487290331055,326.720463091661,328.202897816432,328.084266136988,328.588453114973,409.355898448878,409.097643440222,409.911716063805,410.326066600440,409.023524890055,409.044077532690,410.480495436214,410.066661584182,411.161626138360,406.333062720895,407.701527921002,408.935968427345,460.323350600984,458.133263494380,457.136825780952,457.701459379432,457.680594915214,457.595324430484,465.177258178918,462.217718885221,464.888955825110,461.169583942045,460.011173636885,459.486887298867,461.298892893389,458.316720492664,458.264118899936,462.047384721840,463.100962120089,461.358096029197,461.344652722225,458.816478860923,462.499156207277,461.543886480304,462.184477422192,459.217174699639,473.875006374725,473.588794278000,470.619992440415,474.205415888422,473.260980531030,474.306890705834,471.734858105312,471.811416167328,470.704759943417,475.494596690575,474.375154734321,478.680066248766,508.905862552304,510.323092356121,506.415097825698,509.106100409702,508.465089142008,507.363624685553,506.707081897033,505.782182808081,508.288504290761,505.165506467012,504.332337503043,507.894034458113,464.907164864853,465.111422346378,464.771916157323,464.143831322554,464.012436250578,465.988401671050,462.905408345730,465.210775216702,465.092674913046,469.097252804794,467.476941562167,468.084743126815,461.241909626817,458.914230190048,457.151262516651,462.180191848034,461.095453092476,459.823967048720,464.260919562787,460.690707738104,462.669209209547,462.815709778907,461.809007680764,460.676029020511,437.191555300444,438.930994465328,437.889854161138,439.798367093852,438.238789506838,440.306439739175,444.774353526960,444.050161072628,442.301087227799,440.414583444226,437.380518856295,440.270771460560,454.971441217340,455.654846188356,454.349734044403,456.916058489263,454.819468421398,454.200913299570,460.318616575712,458.592885957561,460.031297867470,455.332463644008,456.442386539738,453.681246338166,457.898497822139,457.396817333882,457.147320060373,460.002712278806,458.031762042065,459.104684333865,459.551773847355,459.343096563369,458.172206455855,456.034053185887,456.962180096669,455.938428120177,470.995043430990,471.848030950712,467.358457933872,467.869776611545,465.771722007706,469.115864783169,474.118826700975,469.767792139643,473.884117022591,466.500820747758,469.621169692110,470.089450339848,470.832989957298,472.486094121083,472.550319466806,478.023119803175,475.692973424768,476.074589141465,471.054398016874,470.391301063946,470.699038311675,475.005291374512,476.538806043345,473.762283081872,504.438004518502,507.334656651529,507.138909946703,502.469856950460,504.033081014236,502.922218903753,504.552295839856,506.198222665779,504.084649862115,500.914324747707,502.845145882172,502.080531792765,513.423876834795,514.405613883020,515.165546377854,514.518880919151,514.355705539564,514.519073417351,511.469655592025,509.274548206422,512.532109319490,509.596385630423,506.342846727996,509.470136303062,466.695727778819,469.220204904506,467.201684623237,467.884611712069,466.968190518057,466.537000861055,465.953063645161,467.900694114593,468.165711087216,467.061721803006,463.873266821851,467.758282309784,449.754788646768,449.076823110607,445.815059799903,443.834988605786,445.284337407335,445.396202124389,447.214678813774,447.415484755106,449.227121762087,443.931586921103,444.580868537974,442.861261065150,438.935632457896,441.809051214612,441.532851462674,438.994980379810,435.381992194464,435.706923164748,438.935786032316,440.657472290600,439.350707726404,441.977326580290,438.469529809805,438.227735083517,409.228977585856,410.109133696159,408.225488376299,406.362658634161,407.014823542067,406.489701901571,412.423206284298,414.057109825589,412.912997452076,412.035900365934,413.232775333649,411.841810974827,360.756201872275,360.505256795744,362.546530533889,364.371819114072,363.187835174402,364.672745888114,360.245318402440,361.158096073750,363.100920847942,365.935952558994,365.192839073185,364.845782949341];
for i=1:24
    pw1(:,12*(i-1)+1:12*i)=pw10(:,i);
    pw2(:,12*(i-1)+1:12*i)=pw20(:,i);
end
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
Psone=sdpvar(1,T,'full');
Perssone=sdpvar(1,T,'full');
uerssone=binvar(1,T,'full');
Pdgone=sdpvar(1,T,'full');
Pidrbone=sdpvar(30,T,'full');
plossone=sdpvar(30,T,'full');
ugone=binvar(NG,T,'full');    
Ppdrone=sdpvar(30,T,'full');
Pidraone=sdpvar(30,T,'full');
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
plosstwo=sdpvar(30,T,'full');
ugtwo=binvar(NG,T,'full');   
Ppdrtwo=sdpvar(30,T,'full');
Pidratwo=sdpvar(30,T,'full');
%约束
caseName = case30;
p30=sum(caseName.bus(:,3));%当前节点功率之和
pbl=pload./p30;
pload1=repmat(pbl,30,1).*repmat(caseName.bus(:,3),1,T);
P_D=0.85.*pload1;%固定负荷

con=[];
con_gen=getConsGen1(Pgone,Pgtwo,mp,T,ugone,ugtwo);%传统机组约束--出力，爬坡，时间
con=[con,con_gen];
con=[con,0<=Pwone<=upwone.*pw1,0<=Pwtwo<=upwtwo.*pw2];%风电约束
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
            Mbw * Pwone(1,ii) - pload1(:,ii) - Ppdrone(:,ii) - Pidraone(:,ii) - Pidrbone(:,ii)-plossone(:,ii)];
    end
    %% 潮流平衡约束，弧度制，有名值
     for ii = 1:T
         con = [con, Pbusone(:,ii) == Bbus * thetaone(:,ii).*baseMVA];
     end
%场景2
    Pbustwo = [];
    for ii = 1:T
        Pbustwo = [Pbustwo, Mbg * Pgtwo(:,ii) + Mbs * Pstwo(1,ii) + Mbe * Persstwo(1,ii) + ...
            Mbw * Pwtwo(1,ii) - pload1(:,ii) - Ppdrtwo(:,ii) - Pidratwo(:,ii) - Pidrbtwo(:,ii)-plosstwo(:,ii)];
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
con=[con,-0.1.*pload1<=Ppdrone<=0.1.*pload1,sum(Ppdrone,2)==0];
con=[con,-0.1.*pload1<=Ppdrtwo<=0.1.*pload1,sum(Ppdrtwo,2)==0];
con=[con,-0.05.*pload1<=Pidraone+Pidrbone<=0.05.*pload1,-0.05.*pload1<=Pidrbone<=0.05.*pload1,-0.05.*pload1<=Pidraone<=0.05.*pload1,sum(Pidraone+Pidrbone,2)==0];
con=[con,-0.05.*pload1<=Pidratwo+Pidrbtwo<=0.05.*pload1,-0.05.*pload1<=Pidratwo<=0.05.*pload1,-0.05.*pload1<=Pidrbtwo<=0.05.*pload1,sum(Pidratwo+Pidrbtwo,2)==0];
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
con=[con,Pwmin<=Psone<=Pwmax];
for i=1:T
    con=[con,Vpmin<=Vs0+sum(Psone(:,1:i),2)<=Vpmax];
end
%场景2
con=[con,Pwmin<=Pstwo<=Pwmax];
for i=1:T
    con=[con,Vpmin<=Vs0+sum(Pstwo(:,1:i),2)<=Vpmax];
end
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
for gt=1:24
    for gi=1:6
        fpg1(gi,gt)=mp(gi,4).*x_pf(gi,gt)+mp(gi,5).*Pgone(gi,gt)+mp(gi,6);
        fpg2(gi,gt)=mp(gi,4).*x_pfz(gi,gt)+mp(gi,5).*Pgtwo(gi,gt)+mp(gi,6);
    end
end
 yyone=binvar(6,23,'full');
 yytwo=binvar(6,23,'full');
 for i=1:23
     for j=1:6
     con=[con,yyone(j,i)<=ugone(j,i+1)];
     con=[con,yyone(j,i)<=1-ugone(j,i)];
     con=[con,yyone(j,i)>=ugone(j,i+1)-ugone(j,i)];
     con=[con,yytwo(j,i)<=ugtwo(j,i+1)];
     con=[con,yytwo(j,i)<=1-ugtwo(j,i)];
     con=[con,yytwo(j,i)>=ugtwo(j,i+1)-ugtwo(j,i)];
     end
 end
  yygone=binvar(1,23,'full');
  yygtwo=binvar(1,23,'full');
 for i=1:23
     con=[con,yygone(1,i)<=upwone(1,i+1)];
     con=[con,yygone(1,i)<=1-upwone(1,i)];
     con=[con,yygone(1,i)>=upwone(1,i+1)-upwone(1,i)];
     con=[con,yygtwo(1,i)<=upwtwo(1,i+1)];
     con=[con,yygtwo(1,i)<=1-upwtwo(1,i)];
     con=[con,yygtwo(1,i)>=upwtwo(1,i+1)-upwtwo(1,i)];
 end
 con=[con,0<=plossone<=pload1,0<=plosstwo<=pload1];
fg=sum(sum(ps(1).*fpg1+ps(2).*fpg2))+0.9*sum(sum(ps(1).*yyone+ps(2).*yytwo));
fess=sum(ps(1)*(0.53*abs(Perssone)+0.12*abs(Perssone))+ps(2)*(0.53*abs(Persstwo)+0.12*abs(Persstwo)))+sum(0.012.*(uersstwo+uersstwo));
fdg=ps(1)*sum(0.02.*Pwone)+ps(2)*(sum(0.02.*Pwtwo)+0.9*sum(ps(1).*yygone+ps(2).*yygtwo))+ps(1)*5.2*sum(pw1-Pwone)+ps(2)*5.2*sum(pw2-Pwtwo);
fload=ps(1)*(sum(100.*abs(Pidraone)+125.*abs(Pidrbone))+15.9.*sum(sum(plossone)))+ps(2)*(sum(100.*abs(Pidratwo)+125.*abs(Pidrbtwo))+15.9.*sum(sum(plosstwo)));
f1=fg+fess+fdg+fload;
ops=sdpsettings('solver','gurobi');
result=optimize(con,f1,ops);
%结果
thetaone=value(thetaone);
Pgone=value(Pgone);
Pwone=value(Pwone);
upwone=value(upwone);
Psone=value(Psone);
Perssone=value(Perssone);
uerssone=value(uerssone);
Pdgone=value(Pdgone);
Pidrbone=value(Pidrbone);
plossone=value(plossone);
ugone=value(ugone);    
Ppdrone=value(Ppdrone);
Pidraone=value(Pidraone);
%场景2变量
thetatwo=value(thetatwo);
Pgtwo=value(Pgtwo);
Pwtwo=value(Pwtwo);
upwtwo=value(upwtwo);
Pstwo=value(Pstwo);
Persstwo=value(Persstwo);
uersstwo=value(uersstwo);
Pdgtwo=value(Pdgtwo);
Pidrbtwo=value(Pidrbtwo);
plosstwo=value(plosstwo);
ugtwo=value(ugtwo);   
Ppdrtwo=value(Ppdrtwo);
Pidratwo=value(Pidratwo);