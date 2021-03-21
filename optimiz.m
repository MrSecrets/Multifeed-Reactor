clc
clear

P0=29*10^5;

ro_b=285.7;
vo=0.22;
phi=0.63;
Dp=0.001;
mu=2.31*10^(-5);
ro_o=21.9;
k_conv=36960320250;
k_cons = 1.034;

D= 0.4;
Ac=pi*D*D/4;
u=vo/Ac;
% u = 1.1210;
G=ro_o*u;

beta=G*(1-phi)/(ro_o*Dp*phi^3)*(150*(1-phi)*mu/Dp + 1.75*G);

latentMeoh = 400*32.04;
cp_meoh=64;
cp_hcl=29.5;
cp_mecl=81.2;
cp_steam=37.6;
cpData = [cp_meoh; cp_hcl; cp_mecl; cp_steam];

delH_298=-34.7*1000;
delCp=cp_steam+cp_mecl-cp_hcl-cp_meoh;

thermSystem = [delH_298; delCp];

T0 = 573.15;
Tmeoh = 273.15+70;
Tout = T0;

global Tvalues
Tvalues= [];

Fa=70.13;
Fa0 = 0;
Fb0= 70.87; %hcl
Fc0= 0;    %mecl
Fd0= 0;     %steam
feedIn = [0.3 0.3 0.4];

feedIn = feedIn/sum(feedIn)

for i = 1:length(feedIn)
    
    
    Fa0 = Fa0 + Fa*feedIn(i); %meoh
    V = (Fa0+Fb0+Fc0+Fb0)*8.314*T0/P0;

    theta = [Fa0/Fa0; Fb0/Fa0; Fc0/Fa0; Fd0/Fa0;];
    
    
    Tin = (theta(1)*cpData(1)*Tmeoh + theta(2)*cpData(2)*Tout + theta(3)*cpData(3)*Tout + theta(4)*cpData(4)*Tout - latentMeoh*theta(1)) ...
        /(theta(1)*cpData(1)+theta(2)*cpData(2)+theta(3)*cpData(3)+theta(4)*cpData(4));
   
    
    zspan = [0 2];
    y0 = [P0; 0];
    
    if i == 1
        Tin = T0;
    end
   
    k = KCons(Tin)/k_conv;
    convConst = Ac*ro_b*k/Fa0;
    
    [z, y]=ode45(@(z, y)ReactorODE(z,y,P0, Tin, beta, convConst, theta, cpData, thermSystem),zspan,y0);
    
    P = y(:,1);
    x = y(:,2);

    figure(10*i+1)
    plot (z,y(:,1)/101325);
    xlabel('z');
    ylabel('Pressure')
    legend('P');

    figure(10*i+2)
    plot (z,y(:,2));
    xlabel('z');
    ylabel('Conversion');
    legend('x');

    termin = find(x>0.99);

    if isempty(termin)

      x_t = x(end);
      P0 = P(end);
      l = z(end);
      Tout = Tvalues(end);

    else          
      indi = termin(1);
      x_t = x(indi);
      P0 = P(indi);
      l = z(indi);
      Tout = Tvalues(indi);
    end

    Fb0 = Fb0 - x_t*Fa0;
    Fc0 = Fc0 + x_t*Fa0;
    Fd0 = Fc0 + x_t*Fa0;
    Fa0 = Fa0 - x_t*Fa0;
    
    Ac = V/u;
    D = sqrt(4*Ac/pi);
    
    fprintf('P = %16.6f \n',P0/ 101325);
    fprintf('x = %16.6f \n',x_t);
    fprintf('z = %16.6f \n',l);
    fprintf('D = %16.6f \n',D);
    fprintf('Tout = %16.6f \n',Tout-273.15);
    fprintf('Tin = %16.6f \n',Tin-273.15);
    
    Tvalues = [];
    
end

function k = KCons(T)

    T= T-273.15;
    
    if T<330
        k = 1.034;
    elseif T<360
        k = 2.202;
    elseif T<390
        k = 3.404;
    else 
        k = 7.714;
    end

end


function Reactor_diff = ReactorODE(z,Yfuncvec, P0, T0, beta,convConst, theta, cpData, thermSystem)
   P = Yfuncvec(1);
   x = Yfuncvec(2);
%    T = Yfuncvec(3);
   global Tvalues 
   
   dPdz = 0 - beta*(P0/P);
   dxdz = convConst*(P/2.01)^(2)*((theta(2)-x))*((theta(1)-x));
   dTdx = (-x*thermSystem(1) ...
        + theta(1)*cpData(1)*T0 + theta(2)*cpData(2)*T0 + theta(3)*cpData(3)*T0 + theta(4)*cpData(4)*T0 ...
        + x*thermSystem(2)*298) ...
        /(theta(1)*cpData(1)+theta(2)*cpData(2)+theta(3)*cpData(3)+theta(4)*cpData(4)+x*thermSystem(2));
    
   Tvalues = [Tvalues, dTdx];
    
Reactor_diff = [dPdz; dxdz];
end


