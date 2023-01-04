function [f]= mymulti1_synN_variable_KD_only_13dec2022(x)

global k Kph N_hd A mu sig KD np hil YY_LD YY_HD

% Variables used in optimization function 
koff=x(1); 
kon=x(1)/KD;

% Model that generates the %Lysis of Healthy (YY_LD) vs % 1/%Lysis of tumor cells (YY_HD) 
[YY_LD,YY_HD]=pareto_model_synN_variable_KD_only_13dec2022(k, Kph, N_hd, A, mu, sig, kon, koff, np, hil );

% The two objective functions 
f(1) =YY_LD;
f(2) =YY_HD; 

end

