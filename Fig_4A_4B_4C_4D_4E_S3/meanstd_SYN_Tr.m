% For statistics extraction from initial synN T cell distribution data files

function [meanTR1,meanTR2,meanTR3,meanTR4,meanTR5,meanTR6,...
    stdTR1,stdTR2,stdTR3,stdTR4,stdTR5,stdTR6,...
    kurtTR1,kurtTR2,kurtTR3,kurtTR4,kurtTR5,kurtTR6]=meanstd_SYN_Tr(divRx)    

dataTR1=[25	1522.84264;
75	3654.822335;
125	2436.548223;
175	1421.319797;
225	609.1370558;
275	253.8071066;
325	101.5228426];
dataTR1(:,1)=dataTR1(:,1).*divRx;
meanTR1=sum(dataTR1(:,1).*dataTR1(:,2))./sum(dataTR1(:,2));
stdTR1=sqrt(sum(((dataTR1(:,1)-meanTR1).^2).*(dataTR1(:,2)))./sum(dataTR1(:,2)));
kurtTR1=sqrt(sum(((dataTR1(:,1)-meanTR1).^4).*dataTR1(:,2))./sum(dataTR1(:,2))); 



dataTR2=[250	5693.069307;
500	     1980.19802;
750	     866.3366337;
1000	495.049505;
1250	371.2871287;
1500	247.5247525;
1750	198.019802;
2000	123.7623762;
2250	24.75247525];
dataTR2(:,1)=dataTR2(:,1).*divRx;
meanTR2=sum(dataTR2(:,1).*dataTR2(:,2))./sum(dataTR2(:,2));
stdTR2=sqrt(sum(((dataTR2(:,1)-meanTR2).^2).*(dataTR2(:,2)))./sum(dataTR2(:,2)));
kurtTR2=sqrt(sum(((dataTR2(:,1)-meanTR2).^4).*dataTR2(:,2))./sum(dataTR2(:,2))); 

dataTR3=[500	3892.944039;
1000	2433.090024;
1500	1581.508516;
2000	851.5815085;
2500	486.6180049;
3000	364.9635036;
3500	243.3090024;
4000	121.6545012;
4500	24.33090024];
dataTR3(:,1)=dataTR3(:,1).*divRx;
meanTR3=sum(dataTR3(:,1).*dataTR3(:,2))./sum(dataTR3(:,2));
stdTR3=sqrt(sum(((dataTR3(:,1)-meanTR3).^2).*(dataTR3(:,2)))./sum(dataTR3(:,2)));
kurtTR3=sqrt(sum(((dataTR3(:,1)-meanTR3).^4).*dataTR3(:,2))./sum(dataTR3(:,2))); 

dataTR4=[500	4136.253041;
1000	2433.090024;
1500	1459.854015;
2000	729.9270073;
2500	486.6180049;
3000	364.9635036;
3500	243.3090024;
4000	121.6545012;
4500	24.33090024];
dataTR4(:,1)=dataTR4(:,1).*divRx;
meanTR4=sum(dataTR4(:,1).*dataTR4(:,2))./sum(dataTR4(:,2));
stdTR4=sqrt(sum(((dataTR4(:,1)-meanTR4).^2).*(dataTR4(:,2)))./sum(dataTR4(:,2)));
kurtTR4=sqrt(sum(((dataTR4(:,1)-meanTR4).^4).*dataTR4(:,2))./sum(dataTR4(:,2))); 

dataTR5=[500	4347.826087;
1000	2415.458937;
1500	1449.275362;
2000	724.6376812;
2500	483.0917874;
3000	241.5458937;
3500	193.236715;
4000	120.7729469;
4500	24.15458937];
dataTR5(:,1)=dataTR5(:,1).*divRx;
meanTR5=sum(dataTR5(:,1).*dataTR5(:,2))./sum(dataTR5(:,2));
stdTR5=sqrt(sum(((dataTR5(:,1)-meanTR5).^2).*(dataTR5(:,2)))./sum(dataTR5(:,2)));
kurtTR5=sqrt(sum(((dataTR5(:,1)-meanTR5).^4).*dataTR5(:,2))./sum(dataTR5(:,2))); 

dataTR6=[500	2094.594595;
1000	2027.027027;
1500	1756.756757;
2000	1216.216216;
2500	945.9459459;
3000	743.2432432;
3500	540.5405405;
4000	405.4054054;
4500	202.7027027;
5000	67.56756757];
dataTR6(:,1)=dataTR6(:,1).*divRx;
meanTR6=sum(dataTR6(:,1).*dataTR6(:,2))./sum(dataTR6(:,2));
stdTR6=sqrt(sum(((dataTR6(:,1)-meanTR6).^2).*(dataTR6(:,2)))./sum(dataTR6(:,2)));
kurtTR6=sqrt(sum(((dataTR6(:,1)-meanTR6).^4).*dataTR6(:,2))./sum(dataTR6(:,2))); 


end
