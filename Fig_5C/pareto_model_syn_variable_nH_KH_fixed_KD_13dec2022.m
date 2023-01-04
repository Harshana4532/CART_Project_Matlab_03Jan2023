function [YY_LD,YY_HD]=pareto_model_syn_variable_nH_KH_fixed_KD_13dec2022(k, Kph, N_hd, A, mu, sig, np, hil)
warning('off')

global  calR2 minS binS maxS KD YY_LD YY_HD koff cutoff_time

clear np_st1 np_st5 kon_s koff_s

divH1=1;       %1/cal; %scaling U distribution w.r.t. H;
divR1=1/calR2; %scaling T distribution w.r.t. R;

%%%%%%%%%%%%%%%%%
koff_hd=koff;
alp=mu;
xb=minS:binS:maxS;

clear sum_simT simT sum_simU simU np_st1 np_st5 np_st
clear Uh H lenH Tr R lenR rho Lam sumT sumU y0 x sumUHL hLIFE_U...
    UHL hLIFE_Ui sumTDT dTIME_T sumTDTi dTIME_Ti rho Lam x dataT Tr R LN RT lenRT

xp=20:100:2000; % For binning of U(t=0) distribution w.r.t. H

% The simulations of Log Normal U(t=0) distributions of Healthy (U2) vs Tumor (U1) cells
U1 = lognpdf(xp,6.9,0.3)./sum(lognpdf(xp,6.9,0.3));
U2 = lognpdf(xp,4.5,0.3)./sum(lognpdf(xp,4.5,0.3));
Uh(1,:)=(0.5.*U2 + 0.5.*U1); % Total U(t=0) distribution 
H=10.^(log(xp));
lenH=length(xp);
H=H./divH1;
Uh(1,:)=Uh(1,:).*20000; % Total No. of U cells =20K

% mean value of initial Log Normal T cells w.r.t. R
alpmu=alp*((10^(6.9)).^np./(hil^np +(10^6.9).^np));

% The simulation of initial Log Normal T cells w.r.t. R
LN = (lognpdf(xb,alpmu,sig)./sum(lognpdf(xb,alpmu,sig))).*10000; % Total No. of T cells =10K
dataTR1(:,1)=xb';
dataTR1(:,2)=LN';

%%%%%%
Tr(1,:)=dataTR1(:,2);
R=dataTR1(:,1)./divR1;
lenR=length(R);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solution of the model
y0(1:lenH,1)=Uh(1,:);
y0((lenH+1):(lenH+lenR),1)=Tr(1,:);  
tspan = [0 cutoff_time];
opt = odeset('RelTol',1e-2,'AbsTol',1e-4);
[tR,yR] = ode45(@(tR,yR) odefcnCONS_LD(tR,yR,Lam,rho,lenH,lenR), tspan, y0);%, opt);
lenHx=3; % Cutoff points of the distribution : healthy vs tumor

% %Lysis of Healthy cells
YY_LD=(1-(sum(yR(end,1:lenHx))./sum(yR(1,1:lenHx)))).*100;
xYY_HD=(1-(sum(yR(end,(1+lenHx):lenH))./sum(yR(1,(1+lenHx):lenH)))).*100;
% 1/%Lysis of Tumor cells
YY_HD=min(1/xYY_HD,10^12);

