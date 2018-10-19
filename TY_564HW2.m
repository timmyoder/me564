clear all; close all; clc

%% 2a
b = [0 1 0;
    0 0 1;
    2 1 -2;];
[v d]= eig(b);

y0A= [1; 1; 1];
y0B = [1; -1; 1];
y0C= [-1; -1; -1];

tspan = 0:.01:10;

[t,yA] = ode45(@(t,y)b*y,tspan,y0A);
[t,yB] = ode45(@(t,y)b*y,tspan,y0B);
[t,yC] = ode45(@(t,y)b*y,tspan,y0C);

figure(1)                    
plot(t,yA,'m','linewidth',2)
title('IC = [1 1 1]')
ylabel('x(t)')
xlabel('t')

figure(2)
plot(t,yB,'m','linewidth',2)
title('IC = [1 -1 1]')
ylabel('x(t)')
xlabel('t')

figure(3)
plot(t,yB,'m','linewidth',2)
title('IC = [-1 -1 -1] - starting on other eigvec')
ylabel('x(t)')
xlabel('t')

%% 2b

c = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    6 5 -5 -5;];
[vc dc]= eig(c)

y0D= [1; 1; 1; 1;];
y0E = [1; -1; 1; -1;];

tspan = 0:.01:10;

[t,yD] = ode45(@(t,y)c*y,tspan,y0D);
[t,yE] = ode45(@(t,y)c*y,tspan,y0E);

figure(4)                     
plot(t,yD,'m','linewidth',2)
title('IC = [1 1 1 1]')
ylabel('x(t)')
xlabel('t')

figure(5)
plot(t,yE,'m','linewidth',2)
title('IC = [1 -1 1 -1]')
ylabel('x(t)')
xlabel('t')

%% 3

e = 0:.005:.5;

w1 = 1;
w2 = 1.5;

ep = 0;

A = [0 0 1 0;
    0 0 0 1;
    -w1^2-ep ep 0 0;
    ep -w2^2-ep 0 0;];

eig3 = eig(A);

for i = 1:size(e,2)-1;
    eig_out = eig(A);
    ep = ep + .005;
    A = [0 0 1 0;
        0 0 0 1;
        -w1^2-ep ep 0 0;
        ep -w2^2-ep 0 0;];
    eig3 = [eig3 eig_out];
end


freq1 = imag(eig3(1,:));
freq2 = imag(eig3(3,:));


figure
plot(e,freq1, 'm', 'linewidth',2)
hold on
plot(e,freq2, 'b', 'linewidth',2)
title('Epsilon vs Frequency')
xlabel('Epsilon')
ylabel('Frequency')
legend('Frequency 1', 'Frequency 2', 'location', 'best')

figure
plot(e,abs(freq1-freq2), 'm','linewidth',2)
title('Epsilon vs Frequency Difference')
xlabel('Epsilon')
ylabel('Frequency 1 - Frequency 2')



