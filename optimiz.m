clc
clear

P0=29*10^5;

ro_b=285.7;
vo=0.22;
phi=0.63;
Dp=0.001;
mu=2.31*10^(-5);
ro_o=21.9;
L=6;
bed = 3;


D= 0.5;
Ac=3.14*D*D/4;
u=vo/Ac;
% u = 1.1210;
G=ro_o*u;

beta=G*(1-phi)/(ro_o*Dp*phi^3)*(150*(1-phi)*mu/Dp + 1.75*G);

Fao=70.13;
k=1.034/36960320250;
% w=Ac*ro_b*z;
x_0=0;


Fa0 = 70.13; %meoh
Fb0= 70.87; %hcl
Fc0= 0;    %mecl
Fd0= 0;     %steam
theta = [Fa0/Fa0; Fb0/Fa0; Fc0/Fa0; Fd0/Fa0;];

cp_meoh=64;
cp_hcl=29.5;
cp_mecl=81.2;
cp_steam=37.6;
cpData = [cp_meoh; cp_hcl; cp_mecl; cp_steam];

delH_298=-34.7*1000;
delCp=cp_steam+cp_mecl-cp_hcl-cp_meoh;

thermSystem = [delH_298; delCp];

convConst = Ac*ro_b*k/Fao;
T0 = 573.15;


zspan = [0 2];
y0 = [P0; 0; T0];
% [z,y]=ode45(@(z, y)ReactorODE,zspan, y0, P0);
[z, y]=ode45(@(z, y)ReactorODE(z,y,P0, T0, beta, convConst, theta, cpData, thermSystem),zspan,y0);

% P = y(:,1)
figure(1)
plot (z,y(:,1)/101325);
xlabel('z');
ylabel('Pressure')
legend('P');

figure(2)
plot (z,y(:,2));
xlabel('z');
ylabel('Conversion');
legend('x');

figure(3)
plot (z,y(:,3));
xlabel('z');
ylabel('Temperature');
legend('T');


fprintf('P = %16.6f \n',y(length(y),1)/ 101325);
fprintf('x = %16.6f \n',y(length(y),2));
fprintf('T = %16.6f \n',y(length(y),3));





function Reactor_diff = ReactorODE(z,Yfuncvec, P0, T0, beta,convConst, theta, cpData, thermSystem)
   P = Yfuncvec(1);
   x = Yfuncvec(2);
   T = Yfuncvec(3);

   
   dPdz = 0 - beta*(P0/P);
   dxdz = convConst*(P/2.01)^(2)*((theta(2)-x))*((theta(1)-x));
   dTdz = (x* -(thermSystem(1)+thermSystem(2)*(T-298)) ...
        + theta(1)*cpData(1)*T0 + theta(2)*cpData(2)*T0 + theta(3)*cpData(3)*T0 + theta(4)*cpData(4)*T0 ...
        + x*thermSystem(2)*298) ...
        /(theta(1)*cpData(1)+theta(2)*cpData(2)+theta(3)*cpData(3)+theta(4)*cpData(4)+x*thermSystem(2));
%    dTdz = (x* -(thermSystem(1)+thermSystem(2)*(T-298)) ...
%         + theta(1)*cpData(1)*T0 + theta(2)*cpData(2)*T0 + theta(3)*cpData(3)*T0 + theta(4)*cpData(4)*T0 + x*thermSystem(2)*T0) ...
%         /(theta(1)*cpData(1) + theta(2)*cpData(2) + theta(3)*cpData(3) + theta(4)*cpData(4) + x*thermSystem(2));



Reactor_diff = [dPdz; dxdz; dTdz];
end



