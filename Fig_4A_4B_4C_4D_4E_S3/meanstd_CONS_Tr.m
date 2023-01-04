% For statistics extraction from initial Constitutive T cell distribution data files

function [meanTR1,meanTR2,...
    stdTR1,stdTR2,kurtTR1,kurtTR2]=meanstd_CONS_Tr(divRx)

%High Aff
dataTR1=[250	6349.809882;
750	1863.117876;
1250	893.5361208;
1750	494.2965774;
2250	247.1482887;
2750	152.0912547;
3250	52.0912547];

dataTR1(:,1)=dataTR1(:,1).*divRx;

meanTR1=sum(dataTR1(:,1).*dataTR1(:,2))./sum(dataTR1(:,2));
stdTR1=sqrt(sum(((dataTR1(:,1)-meanTR1).^2).*(dataTR1(:,2)))./sum(dataTR1(:,2)));
kurtTR1=sqrt(sum(((dataTR1(:,1)-meanTR1).^4).*dataTR1(:,2))./sum(dataTR1(:,2))); 


%Low Aff
dataTR2=[250	3631.71355;
750	1227.621485;
1250	1048.593352;
1750	869.5652178;
2250	818.414323;
2750	767.2634283;
3250	383.6317133;
3750	306.9053703;
4250	204.6035808];


dataTR2(:,1)=dataTR2(:,1).*divRx;

meanTR2=sum(dataTR2(:,1).*dataTR2(:,2))./sum(dataTR2(:,2));
stdTR2=sqrt(sum(((dataTR2(:,1)-meanTR2).^2).*(dataTR2(:,2)))./sum(dataTR2(:,2)));
kurtTR2=sqrt(sum(((dataTR2(:,1)-meanTR2).^4).*dataTR2(:,2))./sum(dataTR2(:,2))); 


%%%%%%%%%%
