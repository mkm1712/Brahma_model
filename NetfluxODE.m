function GOAL= NetfluxODE(Brahma,BMP4L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Brahma signaling minimal model
%Determination of Potential Landscape
% created by Ali Khalilimeybodi 8/1/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% species parameters 

speciesNames = {'Act','CSMES','CSD6','BRM','BRMA','BRG1','POU','BMP4','BMPS','REST','FGF8','GATA4','BMPSI','CAR','NER',}; 
tau = [0.3, 0.3, 0.3, 0.01, 0.3, 0.3, 2, 0.01, 0.3, 0.3, 0.5, 0.5, 0.3, 0.5, 0.5, ]; % Time constants
ymax = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]; %Ymax
 
% reaction parameters 
w = [0, 0, 0, 1, 0, 1, 1, 1, 1, 0.5, 1, 1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]; % Reaction weights
n = [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, ]; %Hill nonlinear parameters
EC50 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, ]; %EC50
rpar = [w;n;EC50];

% Defining species for plot based on speciesNames order
var1=15;
var2=14;
var3=12;
var4=11;
% Initializing
XXX=[];
YYY=[];
VQ=[];
ymax(4)=Brahma;% Set Brahma (1 for WT and 0 for KO cells)

% Run model to steady state to obtain initial values of main simulation
y0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]; 
tspan = [0 10];
w(23)=0.5;
options = [];
params = {rpar,tau,ymax,speciesNames};
[t,y] = ode45(@ODEfun,tspan,y0,options,params); 

% Main Simulation t=0-2 days
y0 = y(end,:);
tspan = [0 2]; 
w(1)=1;
rpar = [w;n;EC50];
params = {rpar,tau,ymax,speciesNames};
options = []; 
[t,y] = ode45(@ODEfun,tspan,y0,options,params); 
delta=diff(t);
dXX=10*diff(y(:,var1))./delta;
dYY=10*diff(y(:,var2))./delta;
dVq=-1*(dXX.^2+dYY.^2).*delta;
Vq1= cumsum ([0;dVq]);
t1=t;
X1=y(:,var1);
Y1=y(:,var2);
C1=y(:,var3);
N1=y(:,var4);

% Plot variations of GATA4 and FGF8 in time
% plot (t, y(:,var3), '--k');
% hold on
% plot (t, y(:,var4), '--r');
% hold on

% Main Simulation t=2-4 days
w(3)=1;
w(5)=BMP4L;% Set BMP4 level (0.25 for WT and 1 for High BMP4)
y0 = y(end,:);
tspan = [0 2]; 
rpar = [w;n;EC50];
params = {rpar,tau,ymax,speciesNames};
options = []; 
[t,y] = ode45(@ODEfun,tspan,y0,options,params); 
delta=diff(t);
dXX=10*diff(y(:,var1))./delta;
dYY=10*diff(y(:,var2))./delta;
dVq=-1*(dXX.^2+dYY.^2).*delta;
Vq2=cumsum ([Vq1(end);dVq]);
t2=t;
X2=y(:,var1);
Y2=y(:,var2);
C2=y(:,var3);
N2=y(:,var4);

% Plot variations of GATA4 and FGF8 in time
% plot (t+2, y(:,var3), '--k');
% hold on
% plot (t+2, y(:,var4), '--r');
% hold on

% Main Simulation t=4-6 days
w(5)=0.0;
y0 = y(end,:);
tspan = [0 2]; 
rpar = [w;n;EC50];
params = {rpar,tau,ymax,speciesNames};
options = []; 
[t,y] = ode45(@ODEfun,tspan,y0,options,params); 
delta=diff(t);
dXX=10*diff(y(:,var1))./delta;
dYY=10*diff(y(:,var2))./delta;
dVq=-1*(dXX.^2+dYY.^2).*delta;
Vq3=cumsum ([Vq2(end);dVq]);
t3=t;
X3=y(:,var1);
Y3=y(:,var2);
C3=y(:,var3);
N3=y(:,var4);

% Plot variations of GATA4 and FGF8 in time
% plot (t+4, y(:,var3), '--k');
% hold on
% plot (t+4, y(:,var4), '--r');

% Main Simulation t=6-10 days
w(2)=1;
y0 = y(end,:);
tspan = [0 4]; 
rpar = [w;n;EC50];
params = {rpar,tau,ymax,speciesNames};
options = []; 
[t,y] = ode45(@ODEfun,tspan,y0,options,params); 
delta=diff(t);
dXX=10*diff(y(:,var1))./delta;
dYY=10*diff(y(:,var2))./delta;
dVq=-1*(dXX.^2+dYY.^2).*delta;
Vq4=cumsum ([Vq3(end);dVq]);
t4=t;
X4=y(:,var1);
Y4=y(:,var2);
C4=y(:,var3);
N4=y(:,var4);
TT=[t1;t2+2;t3+4;t4+6];
XXX=[X1;X2;X3;X4];
YYY=[Y1;Y2;Y3;Y4];
CCC=[C1;C2;C3;C4];
NNN=[N1;N2;N3;N4];
VQ=[Vq1;Vq2;Vq3;Vq4];

% Plot variations of GATA4 and FGF8 in time
% plot (t+6, y(:,var3), '--k');
% hold on
% plot (t+6, y(:,var4), '--r');
% hold on

GOAL=real([XXX YYY TT VQ CCC NNN]);

        

 
function dydt=ODEfun(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames]=params{:}; 
Act = 1; 
CSMES = 2; 
CSD6 = 3; 
BRM = 4; 
BRMA = 5; 
BRG1 = 6; 
POU = 7; 
BMP4 = 8; 
BMPS = 9; 
REST = 10; 
FGF8 = 11; 
GATA4 = 12; 
BMPSI = 13; 
CAR = 14; 
NER = 15; 
dydt = zeros(15,1); 
dydt(Act) = (rpar(1,1)*ymax(Act) - y(Act))/tau(Act); 
dydt(CSMES) = (rpar(1,3)*ymax(CSMES) - y(CSMES))/tau(CSMES); 
dydt(CSD6) = (rpar(1,2)*ymax(CSD6) - y(CSD6))/tau(CSD6); 
dydt(BRM) = (OR(rpar(1,4),inhib(y(BRG1),rpar(:,9)))*ymax(BRM) - y(BRM))/tau(BRM); 
dydt(BRMA) = (AND(rpar(:,6),act(y(CSMES),rpar(:,6)),act(y(BRM),rpar(:,6)))*ymax(BRMA) - y(BRMA))/tau(BRMA); 
dydt(BRG1) = (OR(act(y(CSMES),rpar(:,7)),inhib(y(BRM),rpar(:,8)))*ymax(BRG1) - y(BRG1))/tau(BRG1); 
dydt(POU) = (OR(AND(rpar(:,12),act(y(Act),rpar(:,12)),inhib(y(CSMES),rpar(:,12))),OR(inhib(y(BRMA),rpar(:,13)),inhib(y(BMPS),rpar(:,14))))*ymax(POU) - y(POU))/tau(POU); 
dydt(BMP4) = (rpar(1,5)*ymax(BMP4) - y(BMP4))/tau(BMP4); 
dydt(BMPS) = (OR(act(y(BMP4),rpar(:,10)),act(y(BMPS),rpar(:,11)))*ymax(BMPS) - y(BMPS))/tau(BMPS); %$$$$
dydt(REST) = (OR(AND(rpar(:,16),act(y(CSD6),rpar(:,16)),act(y(BMPS),rpar(:,16)),inhib(y(BMPSI),rpar(:,16))),AND(rpar(:,17),act(y(CSD6),rpar(:,17)),act(y(BRMA),rpar(:,17)),inhib(y(BMPSI),rpar(:,17))))*ymax(REST) - y(REST))/tau(REST); 
dydt(FGF8) = (OR(AND(rpar(:,18),act(y(POU),rpar(:,18)),inhib(y(BMPS),rpar(:,18)),inhib(y(REST),rpar(:,18))),OR(act(y(POU),rpar(:,19)),AND(rpar(:,23),act(y(FGF8),rpar(:,23)),inhib(y(GATA4),rpar(:,23)))))*ymax(FGF8) - y(FGF8))/tau(FGF8); 
dydt(GATA4) = (OR(AND(rpar(:,20),inhib(y(POU),rpar(:,20)),act(y(BMPS),rpar(:,20)),inhib(y(BMPSI),rpar(:,20))),OR(AND(rpar(:,21),act(y(BRMA),rpar(:,21)),act(y(BRG1),rpar(:,21)),inhib(y(BMPSI),rpar(:,21))),AND(rpar(:,22),inhib(y(FGF8),rpar(:,22)),act(y(GATA4),rpar(:,22)))))*ymax(GATA4) - y(GATA4))/tau(GATA4); 
dydt(BMPSI) = (AND(rpar(:,15),act(y(BRMA),rpar(:,15)),act(y(BMPS),rpar(:,15)))*ymax(BMPSI) - y(BMPSI))/tau(BMPSI); 
dydt(CAR) = (OR(AND(rpar(:,26),inhib(y(POU),rpar(:,26)),act(y(BMPS),rpar(:,26)),inhib(y(BMPSI),rpar(:,26))),AND(rpar(:,27),act(y(BRMA),rpar(:,27)),act(y(BRG1),rpar(:,27)),inhib(y(BMPSI),rpar(:,27))))*ymax(CAR) - y(CAR))/tau(CAR); 
dydt(NER) = (OR(AND(rpar(:,24),act(y(POU),rpar(:,24)),inhib(y(BMPS),rpar(:,24)),inhib(y(REST),rpar(:,24))),act(y(POU),rpar(:,25)))*ymax(NER) - y(NER))/tau(NER); 

% utility functions 
function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
    K = (beta - 1).^(1./n); 
    fact = w.*(beta.*x.^n)./(K.^n + x.^n); 
    if fact>w,                 % cap fact(x)<= 1 
        fact = w; 
    end
 
function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    finhib = rpar(1) - act(x,rpar);
 
function z = OR(x,y) 
% OR logic gate 
    z = x + y - x*y;
 
function z = AND(rpar,varargin) 
% AND logic gate, multiplying all of the reactants together 
    w = rpar(1); 
    if w == 0, 
        z = 0; 
    else 
        v = cell2mat(varargin); 
        z = prod(v)/w^(nargin-2);  
    end 
