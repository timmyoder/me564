clear all; close all;

a = [-2 0;
    0 -4;];

[v, d] = eig(a)

%%

dt = 0.01;
x = 0:dt:1;

for i = 1:5
    fun = cos(i*pi()*x);
    plot(x,fun)
    legend('m = 1','m = 2', 'm = 3', 'm = 4', 'm = 5', 'location', 'best')
    xlabel('x')
    hold on
end

%%

test = [1, 4;
    2, 6;
    3, 15;];

out = zeros(3,1);

for i = 1:3
    
v = cos(pi()*x*test(i,1));
w = cos(pi()*x*test(i,2));
Y = v.*w;

out(i) = trapz(Y);

end

%%
m=7;
n=4;

v = cos(pi()*x*m);
w = cos(pi()*x*n);
Y = v.*w;
trapz(Y)

%% 
m=2345
n=25

sin(pi()*(m+n))