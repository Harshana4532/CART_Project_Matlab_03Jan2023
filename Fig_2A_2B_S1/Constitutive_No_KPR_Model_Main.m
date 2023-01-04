%%% Constitutive with No KPR model parameter estimation program and manuscript figures

clear all
warning('off')


%global Kph N_hd 
global divOm_Ha divOs_Ha divOm_La divOs_La yF dat_scal_ConH dat_scal_ConL dkd
global cutoff_time simTime cal calR cal2 calR2 minS binS maxS
global yF 

%%%%%%%%%%%%%%%%%% Const
dkd=1;

cal2=58;                 %H data converion factor       
calR2=1.8; divRx=calR2;  %R data converion factor      

%%%%%%%%%%%%%%%%%% synN
cal=58;                   %H data converion factor 
calR=1.8;  divRy=calR;    %R data converion factor 

simTime=3; % Time length (No. of days) of model simulation
cutoff_time=simTime; % Same simulation time is used as the cutoff time for the analysis of the data
%%%%%%%

% initial value of parameters (Left: notations as per the model in the manuscript) 
k =1.065e-08;    % λc -coefficient of target cells lysis
Kph =0.0065;     % kp -Transition rate of states in the kinetic proof reading model
N_hd =7.3789;    % N  -No. of states formed in the kinetic proof reading model
A =3.7474e-09;   % ρc -coefficint of T cell proliferation
muC =7.2138;     % μc -mean of initial T cell distrobution w.r.t. R -Constitutive model
sigC =0.6;       % σc -stdv of initial T cell distrobution w.r.t. R -Constitutive model
x0=sqrt([k, A, muC, sigC]);% Initial value array
minS=250;binS=500;maxS=10000;  % starting bin level, bin sizes, and ending bin level of initial LN T cell distributions w.r.t R    
rng(914052671); % for repeatability
%%%%%%%

%%%%%%%%%%%% % Lysis data
% Const_High_Affinity
dat_cHD=([3.20	44.79; 4.77	97.75; 5.39	88.13; 5.59	98.46; 6.29	98.10; 7.00	94.45]);
% Const_Low_Affinity
dat_cLD=([3.20	0.34; 4.80	69.12; 5.40	67.67; 5.60	70.37; 6.30	80.10; 7.00	88.18]);

% Means, stdev, and kurtosys of T cell distributions w.r.t. R in Constitutive experiments on day 3
[meanTR1c,meanTR2c,stdTR1c,stdTR2c,kurtTR1c,kurtTR2c]=meanstd_CONS_Tr(divRx);% statistics from initial Constitutive T cell distribution data files
divOm_Ha=stdTR1c; divOs_Ha=kurtTR1c; divOm_La=stdTR2c;divOs_La=kurtTR2c;

% stdev of % lysis data on day 3 converted from standard error of the mean SEM from graphs: 2.stdev=SEM/sqrt(no.of samples) 
dat_scal_ConH=(sqrt(3).*([12.89, 4.89, 16.97, 5.38, 5.06, 5.06]))./2;                
dat_scal_ConL=(sqrt(3).*([14.11, 4.89, 8.81, 6.69, 4.57, 5.22]))./2;                  

% Means of target cells U distributions w.r.t. H in Constitutive High Affinity experiments on day3 
MNT1=([meanTR1c,meanTR1c,meanTR1c,meanTR1c,meanTR1c,meanTR1c])';
% Stdev of target cells U distributions w.r.t. H in Constitutive High Affinity experiments on day3
SGT1=([stdTR1c,stdTR1c,stdTR1c,stdTR1c,stdTR1c,stdTR1c])';
SGT1=SGT1.^2;  % variance of the same above
% Means of target cells U distributions w.r.t. H in Constitutive Low Affinity experiments on day 3
MNT2=([meanTR2c,meanTR2c,meanTR2c,meanTR2c,meanTR2c,meanTR2c])';
% Stdev of target cells U distributions w.r.t. H in Constitutive High Affinity experiments on day 3
SGT2=([stdTR2c,stdTR2c,stdTR2c,stdTR2c,stdTR2c,stdTR2c])';
SGT2=SGT2.^2; % variance of the same above

% Weighted %lysis of target cells, and mean and variance of T cell distributions w,r,t, R on day 3 -High Affinity, 
my1c=[dat_cHD(1:end,2)./dat_scal_ConH';[MNT1./divOm_Ha;SGT1./divOs_Ha]];
% Weighted %lysis of target cells, and mean and variance of T cell distributions w,r,t, R on day 3 -Low Affinity, 
my2c=[dat_cLD(1:end,2)./dat_scal_ConL';[MNT2./divOm_La;SGT2./divOs_La]];
xc=[my1c,my2c]; yc=([my1c,my2c]); % Data handling to insert to functions 

yD=([my1c;my2c])'; 
xD=([my1c;my2c])';

xD_dat_C=([dat_cHD(1:end,2);dat_cLD(1:end,2)])'./[dat_scal_ConH,dat_scal_ConL];
xD_dat_C_M=([MNT1./divOm_Ha;MNT2./divOm_La]);
xD_dat_C_S=([SGT1./divOs_Ha;SGT2./divOs_La]);

% Non linear fit function
[x00,residual,jacobian,covb]=nlinfit(xD,yD,@Objectivefn_Main_CI_deathTf_estconv_ODE45__CONS_Fin,x0);
% Confidence intervals display
C = nlparci(x00,residual,'covar',covb, 'alpha', 0.05).^2 %CI
RSS=sum(sum((residual.^2))); % Residual Sum of Squares display
residuals=residual.^2; % Residual squares vector
x00_st=x00; % Storing parameter values

%%%%%%%%%%%%%  Parameter value display
x00=x00_st;
k=x00(1)^2
A=x00(2)^2
muC=x00(3)^2
sigC=x00(4)^2

%%%%%%%%%%%%%% Preparing const function data for plotting after weighting:  
xD_mod_C=([yF(1:6),yF(19:24)])';
xD_mod_C_M=([yF(7:12),yF(25:30)])';
xD_mod_C_S=([yF(13:18),yF(31:36)])'; 

%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting Cost function without weighting 
%%% Figure S1
figure(1002); clf; 
xu=linspace(10^-2,10^8);
pp0=loglog(xu,xu,'-','LineWidth',2); hold on
xxD_dat_C=xD_dat_C'.*([dat_scal_ConH,dat_scal_ConL])';   
xxD_mod_C=xD_mod_C.*([dat_scal_ConH,dat_scal_ConL])'; 
xxD_dat_C_M=[xD_dat_C_M(1:6).*divOm_Ha;xD_dat_C_M(7:12).*divOm_La]; 
xxD_mod_C_M=[xD_mod_C_M(1:6).*divOm_Ha;xD_mod_C_M(7:12).*divOm_La]; 
xxD_dat_C_S=[xD_dat_C_S(1:6).*divOs_Ha;xD_dat_C_S(7:12).*divOs_La]; 
xxD_mod_C_S=[xD_mod_C_S(1:6).*divOs_Ha;xD_mod_C_S(7:12).*divOs_La];
pp1=loglog(xxD_dat_C,   xxD_mod_C,'o', 'MarkerSize',8, 'linewidth',1); 
pp2=loglog(xxD_dat_C_M, xxD_mod_C_M,'+','MarkerSize',8, 'linewidth',1); 
pp3=loglog(xxD_dat_C_S, xxD_mod_C_S,'^','MarkerSize',8, 'linewidth',1);
legend([pp0 pp1(1) pp2(1) pp3(1)],{'1:1','Const %Lysis','Const T dist means','Const T dist std.'},'Location','northwest')
xlabel('Computations from the data','Fontsize',12)
ylabel('Estimations from the model','Fontsize',12)
box on
title('Fig. S1')

%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %Lysis fit and other figures
%%% Figure 2A
figure(20); clf; hold on
[fig1, xmn_CH, RmnCH, HmnCH]=Plot_synN_HD_ODE45_SYNCONS_CONS_FinCH2(k, A, muC, sigC);
figure(20); hold on
[fig2, xmn_CL, RmnCL, HmnCL]=Plot_synN_LD_ODE45_SYNCONS_CONS_FinCL2(k, A, muC, sigC);
Complexes=([xmn_CH; xmn_CL])';
xlabel('H/cell (log10)', 'Fontsize',12)
ylabel('%Lysis', 'Fontsize',12)
box on

figure(20); hold on
legend([fig1,fig2],{'Const High Aff','Const Low Aff'}...
    ,'Location','northwest')
legend box off
xlabel('HER2/cell (log10)', 'Fontsize',14)
ylabel('%Lysis', 'Fontsize',14)
box on
title('Fig. 2A')


%%%%%%%%%%%%%% Upper CI -%Lysis
x00=sqrt(C(:,2));
k=x00(1)^2;
A=x00(2)^2;
muC=x00(3)^2;
sigC=x00(4)^2;


%%%%%%%%%%%%%% Lower CI -%Lysis
x00=sqrt(C(:,1));
k=x00(1)^2;
A=x00(2)^2;
muC=x00(3)^2;
sigC=x00(4)^2;


%%%%%%%%%% Figure 2B
%Number of complexes at mean R & mean HER2/cell
x00=x00_st;
k=x00(1)^2;
Kph=x00(2)^2;
muC=x00(3)^2;
sigC=x00(4)^2;

KD=17.6*(0.6*113*0.002);
koff_hd=9E-5;
Rmn1= 2534.3; 
clear xmnx Hmn1
vale=3;
for i=1:40;
Hmn1(i)=10^vale;
    xmnx(i)=0.5.*(Rmn1+Hmn1(i)+KD).*(1-sqrt(1-(4.*Rmn1.*Hmn1(i))./((Rmn1+Hmn1(i)+KD).^2)));
vale=vale+0.1;
end
figure(84); clf
plot(log10(Hmn1),xmnx,'o-'); hold on
KD=210*(0.6*113*0.002);
koff_hd=6.8E-4;
Rmn1= 2534.3;
clear xmnx Hmn1
vale=3;
for i=1:40;
Hmn1(i)=10^vale;
    xmnx(i)=0.5.*(Rmn1+Hmn1(i)+KD).*(1-sqrt(1-(4.*Rmn1.*Hmn1(i))./((Rmn1+Hmn1(i)+KD).^2)));
vale=vale+0.1;
end
hold on; figure(84)
plot(log10(Hmn1),xmnx,'o-'); hold on
title('Fig. 2B')
axis([-inf,inf,0,3000])
legend('High Affinity','Low Affinity')
xlabel('HER2 antigen density (molec/cell)')
ylabel('CAR-HER2')







