%% ODE for ODE45

function dydt = odefcnCONS_LD(tR,yR,Lam,rho,lenH,lenR)

dydt = zeros(lenR+lenH,1);
for kp=1:lenH
    for m=1:lenR
        sumU(m)=Lam(m,kp).*yR(lenH+m).*yR(kp);
    end
    dydt(kp)=-sum(sumU);
end
for kp=1:lenR  
    for m=1:lenH
       sumT(m)=rho(kp,m).*yR(lenH+kp).*yR(m) ;
    end
    dydt(lenH+kp)=sum(sumT);  
end
end

