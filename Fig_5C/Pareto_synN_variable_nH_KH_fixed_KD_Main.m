
% Pareto frontier construction between %Lysis of Healthy cells vs 1/%Lysis of Tumor cells, 
% given the synN model compartment of the complex model, 
% with variable KH and nH of the Hill function 
% with fixed KD given by arbitary combinations of koff and kon

clear all
global YY_LD YY_HD cutoff_time

%%%%%%%%%%%%%%%%%% Const
cal2=58;                 %H data converion factor       
calR2=1.8; divRx=calR2;  %R data converion factor      
%%%%%%%%%%%%%%%%%% synN
cal=58;                  %H data converion factor 
calR=1.8; divRy=calR;    %R data converion factor 

% Time length of model simulations
simTime=5; % days
% Point in time used in the analysis 
cutoff_time=simTime; 

%%%%%%%
% Defining %Lysis of Healthy cells & 1/%Lysis of Tumor cells
YY_LD=1; YY_HD=1; 


%%%%%% Best estimated Parameters from the complex model
k =2.0787e-08;  % λc in the model
Kph =0.0072;    % kp in the model
N_hd =7.3735;   % N in the model
A =5.3725e-09;  % ρc in the model
muC =7.2431;    % μc in the model
sigC =0.5691;   % σc in the model
muS =6.4049;    % μs in the model
sigS =1.0999;   % σs in the model
hil =2.4179e+05;% KH in the model
np =3.9683;     % nH in the model

%%%%%%%%%%%%%%%%%%%%%%%%%
mu=muS;
sig=sigS;
%%%%%%%%%%%%%%%%%%%%%%%%%%

% For binning of initial T cell Log Normal distributions w.r.t. R 
minS=250;binS=500;maxS=8000;     

global simTime cal calR cal2 calR2 minS binS maxS 

global k Kph N_hd A mu sig KD koff

% Pareto options in the objective function 
options = optimoptions('gamultiobj','PopulationSize',100,...
          'ParetoFraction',0.7,'PlotFcn',@gaplotpareto);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The choice of kon AND koff values in the simulations 
kon=3.7711e-05;     %  1e-8; %3.7711e-05 ; 17.6*(0.6*113*0.002)
koff=9E-5;    %1.8e-3; %9E-5
%%%%%%%%%
%kon=1e-08;     %  1e-8; %3.7711e-05
%koff=5E-4;    %1.8e-3; %9E-5
%%%%%%%%%
%kon=5e-09;     %  1e-8; %3.7711e-05
%koff=5E-4;    %1.8e-3; %9E-5
%%%%%%%%%
%kon=5e-09;     %  1e-8; %3.7711e-05
%koff=7E-4;    %1.8e-3; %9E-500
%%%%%%%%%
             %kon=1e-09;     %  1e-8; %3.7711e-05
             %koff=7E-4;    %1.8e-3; %9E-5
%%%%%%%%%
%kon=1e-09;     %  1e-8; %3.7711e-05
%koff=15E-4;    %1.8e-3; %9E-5 

%%%%%%%%%%%%%%%%%%%%%%%%%%
KD=koff/kon; % KD value

% Lower and Upper bound for K_H (in Hill function) and n_H (in Hill function) 
lb = [1e03 2]; %[hillx np]
ub = [5e07 7];%[hillx np]
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pareto frontier construction function
[solution,ObjectiveValue] = gamultiobj(@mymulti1_syn_variable_nH_KH_fixed_KD_13dec2022,2,...
                          [],[],[],[],lb,ub,options);

% Plot
figure(1); clf; 
tb6=table(ObjectiveValue(:,1),ObjectiveValue(:,2),solution(:,1)./10^5,solution(:,2));
figure(10); clf; hold on
hold on
sc2=scatter(ObjectiveValue(:,1),ObjectiveValue(:,2),40.*solution(:,2),log10(solution(:,1)),'filled');
axis([0,100,0.01,0.035]);
colorbar
sc2.Marker='o';
grid on;
box on
xlabel('%Lysis Healthy')
ylabel('1/%Lysis Tumor')
title("Fig. 4D")


