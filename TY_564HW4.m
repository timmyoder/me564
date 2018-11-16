clear all; close all; clc;
% Question 4-1

%lorenz parameters
sigma = 10;
beta = 8/3;
rho = 28;

% the three fixed points of lorenz
fp1 = [0, 0, 0];
fp2 = [sqrt(beta*(rho-1)), sqrt(beta*(rho-1)), (rho-1)];
fp3 = [-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), (rho-1);];

 %calculate eig/vectors at each step
fp = fp1; 

A1 = [
    -sigma, sigma, 0;
    rho-fp(3), -1, -fp(1);
    fp(2), fp(1), -beta;
    ];

fp = fp2; 

A2 = [
    -sigma, sigma, 0;
    rho-fp(3), -1, -fp(1);
    fp(2), fp(1), -beta;
    ];

fp = fp3; 

A3 = [
    -sigma, sigma, 0;
    rho-fp(3), -1, -fp(1);
    fp(2), fp(1), -beta;
    ];

[v1, d1] = eig(A1);
[v2, d2] = eig(A2);
[v3, d3] = eig(A3);



%% Question 4-2 - Lorenz's parameters (chaotic)
format short
%lorenz parameters
sigma = 10;
beta = 8/3;

% the three fixed points of lorenz
% fp1 = [0, 0, 0];
% fp2 = [sqrt(beta*(rho-1)), sqrt(beta*(rho-1)), (rho-1)];
% fp3 = [-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), (rho-1);];

dp1 = [0; 0; 0;];
dp2 = dp1;
dp3 = dp1;

for i = 1:10;

    rho = 5*i;
    fp1 = [0, 0, 0];
    fp2 = [sqrt(beta*(rho-1)), sqrt(beta*(rho-1)), (rho-1)];
    fp3 = [-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), (rho-1);];
    
    %calculate eigvalues at each step
    fp = fp1; 

    A1 = [
        -sigma, sigma, 0;
        rho-fp(3), -1, -fp(1);
        fp(2), fp(1), -beta;
        ];

    fp = fp2; 

    A2 = [
        -sigma, sigma, 0;
        rho-fp(3), -1, -fp(1);
        fp(2), fp(1), -beta;
        ];

    fp = fp3; 

    A3 = [
        -sigma, sigma, 0;
        rho-fp(3), -1, -fp(1);
        fp(2), fp(1), -beta;
        ];

    dp1 = [dp1, eig(A1)];
    dp2 = [dp2, eig(A2)];
    dp3 = [dp3, eig(A3)];


end

dp1 = dp1(:,2:end);
dp2 = dp2(:,2:end);
dp3 = dp3(:,2:end);

r = 5*[1:10];

dp1 = [r;
    dp1;
    ];

dp2 = [r;
    dp2;
    ];


dp3 = [r;
    dp3;
    ];
%%
text = ["Rho";
    "lambda 1";
    "lambda 2";
    "lambda 3";]
    

%%
% Initial condition
y0=[10; 10; 10];

sigma = 10;
beta = 8/3;
  
for i = 1:10;
    rho = 5*i;

    % Compute trajectory 
    dt =0.01;
    tspan=[0:dt:20]; 

%     Y(:,1)=y0;
%     yin = y0;
%     for i=1:tspan(2)/dt
%         time = i*dt;
%         yout = rk4singlestep(@(t,y)lorenz(t,y,sigma,beta,rho),dt,time,yin);
%         Y = [Y yout];
%         yin = yout;
%     end
%     
    figure
%         plot3(Y(1,:),Y(2,:),Y(3,:),'b')
%         hold on
        [t,y] = ode45(@(t,y)lorenz(t,y,sigma,beta,rho),tspan,y0);
        plot3(y(:,1),y(:,2),y(:,3),'r')
        title(['Rho = ', num2str(rho)])
        view([33,16])
    
end

%%
a4 = [0 1;
    -1 0;];

eig(a4)




%%
% % Exercise 4-4
clear all; close all;

dt = 0.1;
tspan = [0:dt:30];
y0 = [0.1; -1];
epv = [0.1, 1, 20];
 
for i = 1:3
    eps = epv(i);
    [t,y] = ode45(@(t,y)vanderpol(t,y,eps),tspan,y0); 

    figure
    plot(y(:,1),y(:,2))
    xlabel('y(t)')
    ylabel("y'(t)")
    title(['Eps = ', num2str(eps)])
    
end





%convert row vectors to column vectors
% tcolumn = tspan';
% Ycolumn = Y';


%use ode45 to generate solution to check against
% [t,y] = ode45(@vdp5,[0 50],[1; 1]);
% yODE45 = y;
% compare ode45 solution to rk4 method
% figure
% plot(t,y(:,1),'-o',t,y(:,2),'-o')

% figure
% plot(tcolumn,Ycolumn(:,1),tcolumn,Ycolumn(:,2))
% 
% 

%%RK4 function
function yout = rk4(fun,dt,tk,yk)

f1 = dt*fun(tk,yk);
f2 = dt*fun(tk+dt/2,yk+(1/2)*f1);
f3 = dt*fun(tk+dt/2,yk+(1/2)*f2);
f4 = dt*fun(tk+dt,yk+f3);

yout = yk + 1/6 *(f1+2*f2+2*f3+f4);

end


function dydt = vdp5(t,y)
%VDP1  Evaluate the van der Pol ODEs for mu = 5
%
%from matlab help file

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

dydt = [y(2); 5*(1-y(1)^2)*y(2)-y(1)];
end