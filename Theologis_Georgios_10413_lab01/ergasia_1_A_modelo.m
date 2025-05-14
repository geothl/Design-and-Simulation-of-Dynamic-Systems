clear;
tspan=[0 10];
x0=[0; 0;];
m=8.5;
b=0.65;
k=2;
N=101; %length of tnorm %tnorm=linspace(0,10,N);
tnorm=transpose(0:0.1:10); 
%N=length(tnorm);
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);%better ode simulation
[t,x]=ode45(@sim1_A,tnorm,x0,opts);
y=x(:,1);
l1=2;
l2=5;
sys1=tf([1 0],[1 l1 l2]);
sys2=tf(1,[1 l1 l2]);
u=10*cos(0.5*pi*tnorm)+3;
% sysreal=tf([1/8.5],[1 0.65/8.5 2/8.5]);
% y=lsim(sysreal,u,tnorm,0);   %%Μπορούσαμε να βρούμε το Displacement y(t)
%και από την συνάρτηση μεταφοράς του συστήμτατος καθώς είναι ΜΕΜΕ με
%ελάχιστη υλοποιήση.
z1=-lsim(sys1,y,tnorm,0);
z2=-lsim(sys2,y,tnorm,0);
z3=lsim(sys2,u,tnorm,0);
zeta=[z1 z2 z3];
oros1=0;
oros2=0;
for j=1:1:N
oros1=oros1+(1/N*transpose(zeta(j,:))*zeta(j,:));
oros2=oros2+(1/N*y(j)*transpose(zeta(j,:)));
end
thl=oros1^(-1)*oros2;
thr=thl+[l1;l2;0;];
mek=1/thr(3);
kek=mek*thr(2);
bek=mek*thr(1);
yek=transpose(transpose(thl)*transpose(zeta)); %yek is yhat 
figure;
subplot(2,2,1)
plot(tnorm,y);
title("Real Displacement y(t)")
xlabel("Time in seconds (s)")
ylabel("Displacement in meters (m)");
grid on

subplot(2,2,2)
plot(tnorm,yek);
title("Displacement according to Model $\hat{y}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Displacement in meters (m)");
grid on

subplot(2,2,[3 4])
plot(tnorm,y,'r--',tnorm,yek,'g*');
title("Real Displacement y(t) and Displacement according to Model $\hat{y}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Displacement in meters (m)");
legend({'y(t)','$\hat{y}(t)$'},'Interpreter','latex')
grid on

figure;
e=y-yek;
plot(tnorm,e);
title("Prediction Error e(t)=y(t)-$\hat{y}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Displacement in meters (m)");
grid on

% sysreal=tf([1/mek],[1 bek/mek kek/mek]);
% yok=lsim(sysreal,u,tnorm,0);   %%Μπορούσαμε να βρούμε το Displacement y(t)

