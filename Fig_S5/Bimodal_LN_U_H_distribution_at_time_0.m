%%% Binned bimiodal target cell U(H) probability distribution at time zero

clear all

figure(1); clf

xp=20:10:2000; % For binning of U(t=0) distribution w.r.t. H

% The simulations of Log Normal U(t=0) distributions of Healthy (U2) vs Tumor (U1) cells
U1 = lognpdf(xp,6.9,0.3)./sum(lognpdf(xp,6.9,0.3));
U2 = lognpdf(xp,4.5,0.3)./sum(lognpdf(xp,4.5,0.3));
Uh(1,:)=(0.5.*U2 + 0.5.*U1); % Total U(t=0) distribution 
H=10.^(log(xp));
lenH=length(xp);
Uh(1,:)=Uh(1,:).*20000; % Total No. of U cells =20K

bar(log10(H),Uh(1,:),'BarWidth',2)
xlabel('HER2/cell', 'Fontsize',14)
ylabel('Initial frequency of target cells', 'Fontsize',14)
box on
title('Fig. S5')

