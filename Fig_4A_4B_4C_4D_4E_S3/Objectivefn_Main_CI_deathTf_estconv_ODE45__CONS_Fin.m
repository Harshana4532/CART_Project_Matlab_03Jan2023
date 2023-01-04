% Combined objective function for symultaneous data fitting to above four secnarios  

function [yF]=Objectivefn_Main_CI_deathTf_estconv_ODE45__CONS_Fin(x0,xd)

global yF 

[yF1]=objective_model_all_synN_HD_CI_deathTf_ODE45_SYNCONS_CONS_FinCH(x0); %  Objective function for Constitutive High Affinity data fitting
[yF2]=objective_model_all_synN_LD_CI_deathTf_ODE45_SYNCONS_CONS_FinCL(x0); %  Objective function for Constitutive Low Affinity data fitting
[yF3]=objective_model_all_synN_HD_CI_deathTf_ODE45_SYNCONS_CONS_FinSH(x0); %  Objective function for synN High Affinity data fitting
[yF4]=objective_model_all_synN_LD_CI_deathTf_ODE45_SYNCONS_CONS_FinSL(x0); %  Objective function for synN Low Affinity data fitting

yF=(max(min(reshape([yF1(:,1);yF2(:,1);yF3(:,1);yF4(:,1)],1,[]),100000),0)); 

