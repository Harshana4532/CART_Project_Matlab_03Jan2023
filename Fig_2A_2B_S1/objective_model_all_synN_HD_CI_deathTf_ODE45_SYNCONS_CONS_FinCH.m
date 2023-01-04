function [yF]=objective_model_all_synN_HD_CI_deathTf_ODE45_SYNCONS_CONS_FinCH(x0)

global cutoff_time cal calR minS binS maxS 

global divOm_Ha divOs_Ha dat_scal_ConH 


divH1=1/cal; %scaling U distribution w.r.t. H;
divR1=1/calR;%scaling T distribution w.r.t. R;

%%%%%%%%%%%%%%%%%
% CONS high
data_HD=([3.20	44.79; 4.77	97.75; 5.39	88.13; 5.59	98.46; 6.29	98.10; 7.00	94.45]);
%data_HD=([3.20	40.79; 4.77	95.75; 5.39	88.13; 5.59	98.46; 6.29	98.10; 7.00	94.45]);

KD=17.6*(0.6*113*0.002);
koff_hd=9E-5;

alp=x0(3)^2;
sigma=x0(4)^2;
k=x0(1)^2;
A=x0(2)^2;

xb=minS:binS:maxS;

%%%%%%%%%%%%%%%%
dataUH1=[24	20000];

dataUH2=[500	3168.316832;
1000	3960.39604;
1500	3564.356436;
2000	2574.257426;
2500	1980.19802;
3000	1386.138614;
3500	1089.108911;
4000	792.0792079;
4500	693.0693069;
5000	495.049505;
5500	297.029703];

dataUH3=[2000	4552.352049;
4000	6676.783005;
6000	3945.371775;
8000	2427.921093;
10000	1213.960546;
12000	606.9802731;
14000	303.4901366;
16000	151.7450683;
18000	121.3960546];

dataUH4=[10000	4545.454545;
20000	5303.030303;
30000	3787.878788;
40000	2651.515152;
50000	1515.151515;
60000	946.969697;
70000	568.1818182;
80000	378.7878788;
90000	189.3939394;
100000	113.6363636];

dataUH5=[20000	3880.597015;
40000	6268.656716;
60000	4179.104478;
80000	2388.059701;
100000	1492.537313;
120000	895.5223881;
140000	447.761194;
160000	298.5074627;
180000	149.2537313];


dataUH6=[50000	3296.703297;
100000	6153.846154;
150000	4835.164835;
200000	2857.142857;
250000	1538.461538;
300000	659.3406593;
350000	329.6703297;
400000	219.7802198;
450000	109.8901099];

sig=max(sigma,(0.51-0.49)*rand()+0.49);% setting upper limit 
%%%%%%%%%%%%%%%%

for p=1:6
clear Uh H lenH Tr R lenR rho Lam sumT sumU y0 x sumUHL hLIFE_U...
    UHL hLIFE_Ui sumTDT dTIME_T sumTDTi dTIME_Ti rho Lam x dataTR1 Tr R LN

switch p
    case 1 
Uh(1,:)=dataUH1(:,2);
H=dataUH1(:,1);
lenH=length(H);
    case 2 
Uh(1,:)=dataUH2(:,2);
H=dataUH2(:,1);
lenH=length(H); 
    case 3 
Uh(1,:)=dataUH3(:,2);
H=dataUH3(:,1);
lenH=length(H);
    case 4 
Uh(1,:)=dataUH4(:,2);
H=dataUH4(:,1);
lenH=length(H);
    case 5 
Uh(1,:)=dataUH5(:,2);
H=dataUH5(:,1);
lenH=length(H);
    case 6 
Uh(1,:)=dataUH6(:,2);
H=dataUH6(:,1);
lenH=length(H);
end

H=H./divH1;
Uh(1,:)=Uh(1,:);

alpmu=alp;

%%%%%%
LN = (lognpdf(xb,alpmu,sig)./sum(lognpdf(xb,alpmu,sig))).*10000;
        
dataTR1(:,1)=xb';
dataTR1(:,2)=LN';
%%%%%%
Tr(1,:)=dataTR1(:,2);
R=dataTR1(:,1)./divR1;
lenR=length(R);
Tr(1,:)=Tr(1,:);

%%%%%%
% x [HR]complexes 
for i=1:lenR
    for j=1:lenH
    x(i,j)=0.5.*(R(i)+H(j)+KD).*(1-sqrt(1-(4.*R(i).*H(j))./((R(i)+H(j)+KD).^2)));
    Lam(i,j)=k.*(x(i,j));
    end
end
for i=1:lenR
    for j=1:lenH
    rho(i,j)=A.*(x(i,j));; 
    end
end


y0(1:lenH,1)=Uh(1,:);
y0((lenH+1):(lenH+lenR),1)=Tr(1,:);  
tspan = [0 cutoff_time];
opt = odeset('RelTol',1e-2,'AbsTol',1e-4);
[tR,yR] = ode45(@(tR,yR) odefcnCONS_LD(tR,yR,Lam,rho,lenH,lenR), tspan, y0);

YY=(1-(sum(yR(end,1:lenH))./sum(yR(1,1:lenH)))).*100;

simT=yR(end,lenH+1:end);

%%%%%%%
meanTR1(p)=sum(xb.*simT)./sum(simT);
stdTR1(p)=(sum(((xb-meanTR1(p)).^2).*(simT)./sum(simT)));
%%%%%%%

mYY(p)=YY;
end

yF=[(mYY)./dat_scal_ConH, (meanTR1)./divOm_Ha, ((stdTR1)./divOs_Ha)]';


