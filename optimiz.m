clc
clear

global Po;
Po=29*10^5;
ro_b=285.7;
vo=0.22;
phi=0.63;
Dp=0.001;
mu=2.31*10^(-5);
ro_o=21.9;
L=6;
bed = 3;


D= 0.3;
Ac=3.14*D*D/4;
u=vo/Ac;
G=ro_o*u;

global beta;
global gamma;

beta=G*(1-phi)/(ro_o*Dp*phi^3)*(150*(1-phi)*mu/Dp + 1.75*G)
Fao=70.13;
k=1.034/36960320250;
% w=Ac*ro_b*z;
x_0=0;
gamma = Ac*ro_b*k/Fao


zspan = [0 1.2];
y0 = [Po; 0];
[z,y]=ode45(@ODEfun,zspan, y0);

% P = y(:,1)
figure(1)
plot (z,y(:,1)/101325);
xlabel('z');
ylabel('Pressure')
legend('A');

figure(2)
plot (z,y(:,2));
xlabel('z');
ylabel('Conversion');
legend('B');

fprintf('A = %16.6f \n',y(length(y),1));
fprintf('B = %16.6f \n',y(length(y),2));

function dYfuncvecdt = ODEfun(t,Yfuncvec)
   A = Yfuncvec(1);
   B = Yfuncvec(2);
   global beta;
   global gamma;
   global Po;
   dAdt = 0 - beta*(Po/A);
   dBdt = gamma*A*A*(1.01-B)*(1-B)/(2.01^2);
   dYfuncvecdt = [dAdt; dBdt];
end