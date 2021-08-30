% Determine cell and BMP4 level
Brahma=0;% Set Brahma (1 for WT and 0 for KO cells)
BMP4L=0.25; % Set BMP4 level (0.25 for WT and 1 for High BMP4)
% Model
GOAL= NetfluxODE(Brahma,BMP4L);
%  Run single simulation 
y0 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
tspan = [0 10]; 
k1=2;
k2=1;
k3=1;
k4=1;
k5=2;
k6=1;
k7=1;
k8=1;
k9=2;
tt=find (GOAL(:,3)==6); % Transcription Factors activity at Day=6
XX=GOAL(tt(1),1);
YY=GOAL(tt(1),2);
if XX-YY> 0.2
    XX=1;
    YY=0;
end
if YY-XX> 0.2
    YY=1;
    XX=0;
end
rpar = [k1,k2,k3,k4,k5,k6,k7,k8,k9,XX,YY];
params = {rpar};
options = []; 
[t,y] = ode23(@ODEfun,tspan,y0,options, params); 
plot(t,y); 
legend('ESC','EB','MES','CP','CM','NS','NP', 'NCI1', 'OST', 'FIB');
hold on
function dydt=ODEfun(t,y, params) 
dydt = zeros(10,1);
rpar=params{:}; 
ESC=1;EB=2;MES=3;CP=4;CM=5;NS=6;NP=7;NCI1=8;OST=9;FIB=10;
XX=rpar(10);
YY=rpar(11);
J1=rpar(1)*y(ESC);
J2=rpar(2)*y(EB)*heaviside(t-2)*(t-2);
J3=rpar(3)*y(MES)*heaviside(t-4)*(t-4)*YY;
J4=rpar(4)*y(CP)*heaviside(t-6)*(t-6);
J5=rpar(5)*y(MES)*heaviside(t-4)*(t-4)*XX;
J6=rpar(5)*y(NS)*heaviside(t-6)*(t-6);
J7=rpar(7)*y(MES)*heaviside(t-4)*(t-4)*XX;
J8=rpar(8)*y(NCI1)*heaviside(t-6)*(t-6);
J9=rpar(9)*y(NCI1)*heaviside(t-6)*(t-6);
dydt(ESC) = -J1; 
dydt(EB) = J1-J2; 
dydt(MES) = J2-J5-J3-J7; 
dydt(CP) = J3-J4; 
dydt(CM) = J4; 
dydt(NS) = J5-J6;
dydt(NP) = J6;
dydt(NCI1) =J7-J8-J9;
dydt(OST) = J8;
dydt(FIB) = J9;
end