function [Pgtwo23,Pidrdtwo23]=main3(upwtwo23,Persstwo23,Pidrbtwo23,pw23,ugtwo23,Pstwo23,Pidratwo23,Ppdrtwo23,pload)
% k=2;
T=1;%ʱ��
caseName = case30;
p30=sum(caseName.bus(:,3));%��ǰ�ڵ㹦��֮��
pbl=pload./p30;
pload23=repmat(pbl,30,1).*repmat(caseName.bus(:,3),1,T);
% P_D=0.85.*pload1;%�̶�����
% for i=4*(k-1)+1:4*k
%     for j=4*(i-4*(k-1)-1)+1:4*(i-4*(k-1))
%         pw12(:,j)=pw1(:,i);
%         pw22(:,j)=pw2(:,i);
%         pload1(:,j)=pload1(:,i);
%     end
% end
% for i=1:24
%     for j=4*(i-1)+1:4*i
%     ugone1(:,j)=ugone(:,i);
%     ugtwo1(:,j)=ugtwo(:,i);
%     Psone1(:,j)=Psone(:,i);
%     Pstwo1(:,j)=Pstwo(:,i);
%     Ppdrone1(:,j)=Ppdrone(:,i);
%     Ppdrtwo1(:,j)=Ppdrtwo(:,i);
%     Pidraone1(:,j)=Pidraone(:,i);
%     Pidratwo1(:,j)=Pidratwo(:,i);
%     pw1(:,j)=pw10(:,i);
%     pw2(:,j)=pw20(:,i);
%     end
% end
% ugone1=ugone;ugtwo1=ugtwo;Psone1=Psone;Pstwo1=Pstwo;Ppdrone1=Ppdrone;Ppdrtwo1=Ppdrtwo;Pidraone1=Pidraone;Pidratwo1=Pidratwo;

% ugone2=zeros(6,T);ugtwo2=zeros(6,T);
% Psone2=zeros(1,T);
% Ppdrone2=zeros(30,T);
% Pidraone2=zeros(30,T);
% Pstwo2=zeros(1,T);
% Ppdrtwo2=zeros(30,T);
% Pidratwo2=zeros(30,T);
% % for i=4*(k-1)+1:4*k
% %     for j=4*(i-4*(k-1)-1)+1:4*(i-4*(k-1))
%         ugone2(:,1:T)=ugone1(:,16*(k-1)+1:16*k);
%         Psone2(:,1:T)=Psone1(:,16*(k-1)+1:16*k);
%         Ppdrone2(:,1:T)=Ppdrone1(:,16*(k-1)+1:16*k);
%         Pidraone2(:,1:T)=Pidraone1(:,16*(k-1)+1:16*k);
%         ugtwo2(:,1:T)=ugtwo1(:,16*(k-1)+1:16*k);
%         Pstwo2(:,1:T)=Pstwo1(:,16*(k-1)+1:16*k);
%         Ppdrtwo2(:,1:T)=Ppdrtwo1(:,16*(k-1)+1:16*k);
%         Pidratwo2(:,1:T)=Pidratwo1(:,16*(k-1)+1:16*k);
%         pw12(:,1:T)=pw1(:,16*(k-1)+1:16*k);
%         pw22(:,1:T)=pw2(:,16*(k-1)+1:16*k);
%     end
% end
%����
Ns=2;%������Ϊ2
ps=[0.3 0.7];%ÿ�ֳ�������
Nerss=1;%���ܵ�վ����
Nbt=1;%�绯ѧ���ܵ�վ����
PIbt=2.3;%�绯ѧ���ܵ�λʱ���۾ɳɱ�
NG=6;%�����������
NDG=1;%DG����
kc=2.5;%����ͷ�
socmin=0.1;
socmax=0.9;
Pech=-50;%��繦������
Pedis=50;%�ŵ繦������
Se=200;%������������
Se0=100;
Pwmin=-50;%��ˮ��վ��С����
Pwmax=50;%��ˮ��վ�����
Vs0=100;
Vpmin=0;%��ˮ������Сֵ
Vpmax=200;%��ˮ�������ֵ
Ca=100;Cb=125;Cc=150;Cd=150;%IDR�����ɱ�ϵ��
mp=[1 150 50 0.0375 20 372.5 72 2;%����������
    2 60  20  0.175 17.5 352.3 48 2;
    3 60  15  0.625 10 316.5  30 2;
    4 50  10  0.0834 32.5 329.2 30 2;
    5 40  10  0.25  30  276.4 18 2;
    6 45  12  0.25 30 232.2  24 2];
%pload=[390 382 354 340 328 409 460 460 472 511 465 458 442 456 457 471 475 503 512 469 445 438 410 362];%����

idg=[1 2 5 8 11 13];
iw=2;
ie=2;
is=8;
%����
thetatwo23=sdpvar(30,T,'full');
Pgtwo23=sdpvar(NG,T,'full');
Pwtwo23=sdpvar(1,T,'full');
% upwone=binvar(1,T,'full');
% Perssone=sdpvar(1,T,'full');
% uerssone=binvar(1,T,'full');
% Pidrbone=sdpvar(30,T,'full');
Pidrctwo23=sdpvar(30,T,'full');
Pidrdtwo23=sdpvar(30,T,'full');
plosstwo23=sdpvar(30,T,'full');
rsp=sdpvar(1,T,'full');
rsb=sdpvar(1,T,'full');
rgp=sdpvar(NG,T,'full');
rgb=sdpvar(NG,T,'full');
%����2����
% thetatwo=sdpvar(30,T,'full');
% Pgtwo=sdpvar(NG,T,'full');
% Pwtwo=sdpvar(1,T,'full');
% % upwtwo=binvar(1,T,'full');
% % Pstwo=sdpvar(1,T,'full');
% % Persstwo=sdpvar(1,T,'full');
% % uersstwo=binvar(1,T,'full');
% % Pidrbtwo=sdpvar(30,T,'full');
% Pidrctwo=sdpvar(30,T,'full');
% plosstwo=sdpvar(30,T,'full');
%Լ��
con=[];
%% ��ȡ����Լ��
Pgmin=mp(:,3);
Pgmax=mp(:,2);
cons = [];
% 1. ����������Լ��
cons = [cons, ugtwo23.*repmat(Pgmin,1,T) <=Pgtwo23<= ugtwo23.*repmat(Pgmax,1,T)];
con=[con,0<=Pwtwo23<=upwtwo23.*pw23];%���Լ��
%����Լ��
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
    %% ����Pbus
    Mbg = getMbgMatrix(idg',bus);
    Mbe = getMbdMatrix(ie,bus);
    Mbw = getMbdMatrix(iw,bus);
    Mbs = getMbdMatrix(is,bus);
    Pbustwo23 = [];
    for ii = 1:T
        Pbustwo23 = [Pbustwo23, Mbg * Pgtwo23(:,ii) + Mbs * Pstwo23(1,ii) + Mbe * Persstwo23(1,ii) + ...
            Mbw * Pwtwo23(1,ii) - pload23(:,ii) - Ppdrtwo23(:,ii) - Pidratwo23(:,ii) - Pidrbtwo23(:,ii)- Pidrctwo23(:,ii) - Pidrdtwo23(:,ii) - plosstwo23(:,ii)];
    end
    %% ����ƽ��Լ���������ƣ�����ֵ
     for ii = 1:T
         con = [con, Pbustwo23(:,ii) == Bbus * thetatwo23(:,ii).*baseMVA];
     end
    %����Լ��
for j=1:T
d_thetatwo(j,:) = thetatwo23(f,j) - thetatwo23(t,j);
pftwo = d_thetatwo(j,:)' ./ x .* baseMVA; 
con = [con, -500 <=  pftwo <=500];
end
con = [con,-pi <= thetatwo23 <= pi];
%DRԼ��
con=[con,-0.05.*pload23<=Pidratwo23+Pidrbtwo23+Pidrctwo23<=0.05.*pload23,-0.05.*pload23<=Pidrctwo23<=0.05.*pload23,-0.03.*pload23<=Pidrdtwo23<=0.03.*pload23];
%��ת����Լ��
alfab=0.95;gama=0.8;
 con=[con,pload23(:,ii)+Ppdrtwo23(:,ii)+Pidratwo23(:,ii)+Pidrbtwo23(:,ii)+Pidrctwo23(:,ii)+Pidrdtwo23(:,ii)+plosstwo23(:,ii)-(sum(Pgtwo23(:,ii))+Pstwo23(1,ii)+Persstwo23(1,ii) + ...
            Pwtwo23(1,ii))+(2-2*alfab)*(rsb+sum(rgb))+(2*alfab-1)*gama*(rsb+sum(rgb))<=0];
 con=[con,pload23(:,ii)+Ppdrtwo23(:,ii)+Pidratwo23(:,ii)+Pidrbtwo23(:,ii)+Pidrctwo23(:,ii)+Pidrdtwo23(:,ii)+plosstwo23(:,ii)-(sum(Pgtwo23(:,ii))+Pstwo23(1,ii)+Persstwo23(1,ii) + ...
            Pwtwo23(1,ii))-(2-2*alfab)*(rsp+sum(rgp))-(2*alfab-1)*gama*(rsp+sum(rgp))<=0];
 con=[con,Pgtwo23+rgp<=Pgmax,Pgtwo23-rgb>=Pgmin,Pstwo23+rsp<=Pwmax,Pstwo23-rsb>=Pwmin,rgp>=0,rsp>=0,rgb>=0,rsb>=0];
%Ŀ�꺯��
%����������гɱ�a2*p+bp+c�ֶ����Ի�
tn=T;
Pmax=[mp(1,2),mp(2,2),mp(3,2),mp(4,2),mp(5,2),mp(6,2)];
Pmin=[mp(1,3),mp(2,3),mp(3,3),mp(4,3),mp(5,3),mp(6,3)];
gn=5;%�ֶκ������Ի�����ͬ
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
con = [con, Pgtwo23(1,:)==gl2z(1,:)*gw1z];
con = [con, Pgtwo23(2,:)==gl2z(2,:)*gw2z];
con = [con, Pgtwo23(3,:)==gl2z(3,:)*gw3z];
con = [con, Pgtwo23(4,:)==gl2z(4,:)*gw4z];
con = [con, Pgtwo23(5,:)==gl2z(5,:)*gw5z];
con = [con, Pgtwo23(6,:)==gl2z(6,:)*gw6z];
con = [con, gw1z>=0];
con = [con, gw2z>=0];
con = [con, gw3z>=0];
con = [con, gw4z>=0];
con = [con, gw5z>=0];
con = [con, gw6z>=0];    
for gt=1:T
    for gi=1:6
        fpg2(gi,gt)=mp(gi,4).*x_pfz(gi,gt)+mp(gi,5).*Pgtwo23(gi,gt)+mp(gi,6);
    end
end
 con=[con,0<=plosstwo23<=pload23];
fg=sum(sum(fpg2));
fess=sum((0.53*abs(Persstwo23)+0.12*abs(Persstwo23)));%+sum(0.012.*uersstwo23);
fdg=sum(0.02.*Pwtwo23)+5.2*sum(pw23-Pwtwo23);
fload=sum(sum(150.*abs(Pidrctwo23)+150.*abs(Pidrdtwo23)))+15.9.*sum(sum(plosstwo23));
f1=fg+fess+fdg+fload;
ops=sdpsettings('solver','gurobi');
result=optimize(con,f1,ops);




