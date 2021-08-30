%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deterministic rate equations
% dxdt = c1 + a1*(x^n)/(Kdxx^n + (x^n)) + (b1*(Kdyx^n))/(Kdyx^n + (y^n)) - (x*k1)            
% dydt = c1 + a2*(y^n)/(Kdyy^n + (y^n)) + (b2*(Kdxy^n))/(Kdxy^n + (x^n)) - (y*k2)
% created by Ali Khalilimeybodi 8/1/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine cell and BMP4 level
Brahma=1;% Set Brahma (1 for WT and 0 for KO cells)
BMP4L=0.25; % Set BMP4 level (0.25 for WT and 1 for High BMP4)
% Model
GOAL= NetfluxODE(Brahma,BMP4L);

% %    Making Video clip of Phase Planes
%    vvv = VideoWriter('KO.mp4','MPEG-4');
%    vvv.FrameRate=5;
%    vvv.Quality=100;
%    open(vvv);
%    iteration = 40;
%    FF(iteration) = struct('cdata',[],'colormap',[]);

for ii=1:10:1100 % phase plane of system from Day 0-10
% Tristable gene model parameters
a1 = 10*GOAL(ii,2);
a2 = 10*GOAL(ii,1);
b1=1;
b2=0.5;
n=2;
c1=0;
c2=0;
k1=1;
k2=1;
Kdxx=0.5;
Kdyx=0.5;
Kdyy=0.5;
Kdxy=0.5;
f = @(t,Y) [(a1*Y(1)^n)/(Kdxx^n + Y(1)^n) + b1*(Kdyx^n)/(Kdyx^n + Y(2)^n) - Y(1)*k1+ c1; (a2* Y(2)^n)/(Kdyy^n + Y(2)^n) + b2*(Kdxy^n)/(Kdxy^n + Y(1)^n) - Y(2)+ c2];
% for plotting the vector field using quiver
y1 = linspace(0,10,51);
y2 = linspace(0,10,51);
[x,y] = meshgrid(y1,y2);
u = zeros(size(x));
v = zeros(size(x));
t=0;
axes1 = axes('Parent',figure);

for i = 1:numel(x)
    % Calculating gradient value for a and b
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

quiver(x,y,u,v,1.5,'r');
hold on

grid on
% solve the equation for which gradient is zero and plotting it
syms x y

fimplicit(@(x,y)(a1* x^n)/(0.5^n + x^n) + b1*(0.5^n)/(0.5^n + y^n) - x,'b',[0 10 0 10]);
hold on
fimplicit(@(x,y)(a2* y^n)/(0.5^n + y^n) + b2*(0.5^n)/(0.5^n + x^n) - y,'k',[0 10 0 10]);
hold off
xlim(axes1,[0 10]);
ylim(axes1,[0 10]);
set(axes1,'XTickLabel',...
    {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'},...
    'YTickLabel',...
    {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
xlabel('GATA4');
ylabel('FGF8');
% %    Making Video clip of Phase Planes
%      FF(ii) = getframe(gcf);
%      frame= getframe(gcf);
%      writeVideo(vvv,frame);
end
%    close(vvv);
%    close all