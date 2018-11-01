ah = [-0.01 0;
    1 -0.011;];
   

[v,d]=eig(ah)

y0A= [0.1;
    0];


tspan = 0:.1:1000;

[t,yA] = ode45(@(t,y)ah*y,tspan,y0A);

figure                
plot(t,yA(:,1),'m','linewidth',2), hold on
plot(t,yA(:,2),'b','linewidth',2)
title('IC = [.1 0]')
% ylabel('x(t)')
xlabel('t')
legend('x(t)', '(v(t)', 'location', 'best')


figure                
plot(yA(:,1),yA(:,2),'m','linewidth',2)
title('IC = [.1 0]')
xlabel('x(t)')
ylabel('v(t)')

T = .5*yA(:,2).^2;

figure                
plot(t,T,'m','linewidth',2)
title('IC = [.1 0]')
ylabel('Kinetic Energy')
xlabel('t')

%%
A = [0 1;
    -8 -6;];

[V,D] = eig(A)

%%
A = [0 1;
    -1 0;];

[V,D] = eig(A)
