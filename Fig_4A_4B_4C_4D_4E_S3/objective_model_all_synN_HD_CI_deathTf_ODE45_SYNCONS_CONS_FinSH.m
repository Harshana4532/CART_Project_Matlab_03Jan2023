%%%%%  Objective function for synN High Affinity data fitting

function [yF]=objective_model_all_synN_HD_CI_deathTf_ODE45_SYNCONS_CONS_FinSH(x0)

global cutoff_time cal calR minS binS maxS  
global divOm_Ha_syn divOs_Ha_syn dat_scal_synH

divH1=1/cal;  %scaling U distribution w.r.t. H;
divR1=1/calR; %scaling T distribution w.r.t. R;

%%%%%%%%%%%%%%%%%
%synN High Affinity %Lysis
data_HD=([3.20	0.45; 4.77	30.56; 5.39	29.38; 5.59	69.88; 6.29	92.14; 7.00	79.08]);

% K_D disassociatoin rate
KD=17.6*(0.6*113*0.002);% 
koff_hd=9E-5;

% Parameters
k=x0(1)^2;
Kph=x0(2)^2;
N_hd=x0(3)^2;
A=x0(4)^2;
alp=x0(7)^2;
sigma=x0(8)^2;
hil=x0(9).^2;
np=x0(10).^2;

xb=minS:binS:maxS;  % bins for initial LN T cell distribution w.r.t. R

%%%%%%%%%%%%%%%% U(H) data

dataUH1=[25	2258.064516;
50	3548.387097;
75	4032.258065;
100	3709.677419;
125	2741.935484;
150	1774.193548;
175	967.7419355;
200	483.8709677;
225	241.9354839;
250	161.2903226;
275	80.64516129];

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

sig=max(sigma,(1.52-1.48)*rand()+1.48); % setting an upper limit
%%%%%%%%%%%%%%%%%%%

for p=1:6

clear Uh H lenH TrR lenR rho Lam sumT sumU y0 x sumUHL hLIFE_U...
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

alpmu=alp*((10^data_HD(p,1)).^np./(hil^np+(10^data_HD(p,1)).^np)); % mean for initial LN T distribution w.r.t. R

%%%%%% Initial LN T cell distribution w.r.t. R 
LN = (lognpdf(xb,alpmu,sig)./sum(lognpdf(xb,alpmu,sig))).*10000;
        
dataTR1(:,1)=xb';
dataTR1(:,2)=LN';
%%%%%%
Tr(1,:)=dataTR1(:,2);
R=dataTR1(:,1)./divR1;
lenR=length(R);
Tr(1,:)=Tr(1,:);


%%%%%%

% [HR]complexes 
for i=1:lenR
    for j=1:lenH
    bet=1+Kph/koff_hd; KDp=(KD+(Kph/(koff_hd/KD)))/bet; 
    x(i,j)=0.5.*(R(i)+H(j)+KDp).*(1-sqrt(1-(4.*R(i).*H(j))./((R(i)+H(j)+KDp).^2)));
    Lam(i,j)=k.*(x(i,j)).*(Kph./(Kph+koff_hd)).^N_hd;
    end
end
for i=1:lenR
    for j=1:lenH
    rho(i,j)=A.*(x(i,j)).*(Kph./(Kph+koff_hd)).^N_hd; 
    end
end


%%%%%%%%%%% ODE45 numerical solver
y0(1:lenH,1)=Uh(1,:);
y0((lenH+1):(lenH+lenR),1)=Tr(1,:);  
tspan = [0 cutoff_time];
opt = odeset('RelTol',1e-2,'AbsTol',1e-4);
[tR,yR] = ode45(@(tR,yR) odefcnCONS_LD(tR,yR,Lam,rho,lenH,lenR), tspan, y0);%, opt);
YY=(1-(sum(yR(end,1:lenH))./sum(yR(1,1:lenH)))).*100;
simT=yR(end,lenH+1:end);

%%%%%%%
if p>1
meanTR1(p-1)=sum(R.*simT')./sum(simT);
stdTR1(p-1)=(sum(((R-meanTR1(p-1)).^2).*(simT')./sum(simT)));
mYY(p-1)=YY;
end
 
end

yF=[(mYY)./dat_scal_synH, (meanTR1)./divOm_Ha_syn', (stdTR1)./divOs_Ha_syn']'; % for Weighted Cost function

