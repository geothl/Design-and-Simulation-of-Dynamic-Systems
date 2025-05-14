clear;
%N=10001; 
t=transpose(0:0.01:10);
% t=transpose(0:0.01:100);
N=length(t);
u1=3*sin(pi*t);
u2=2.5*ones(length(t),1);
vr=zeros(length(t));
vc=zeros(length(t));
[vr,vc]=v(t);
% vr=vr+normrnd(0,0.1,length(t),1);
% vc=vc+normrnd(0,0.1,length(t),1);
l1=2;
l2=5;
sys1=tf([1 0],[1 l1 l2]);
sys2=tf(1,[1 l1 l2]);
z1=lsim(sys1,vr,t,vr(1));
z2=lsim(sys2,vr,t,vr(1));
z3=lsim(sys1,vc,t,vc(1));
z4=lsim(sys2,vc,t,vc(1));
z5=lsim(-sys2,u1,t,u1(1));

%if more compact implementation is used
% z1=lsim(sys1,vr,t,vr(1));
% z2=lsim(sys2,vc,t,vc(1));
% z3=lsim(sys1,vc,t,vc(1));
% z4=lsim(sys2,vr,t,vr(1));
% zeta=[z1 z2 z3 z4 z5];

zeta=[z1 z2 z3 z4 z5];
oros1=0;
oros2=0;
for j=1:1:N
oros1=oros1+(1/N*transpose(zeta(j,:))*zeta(j,:));
oros2=oros2+(1/N*vc(j)*zeta(j,:));
end
thl=oros1^(-1)*transpose(oros2);
CR_inv=thl(1); %CR_inv=1/CR
CL_inv=0.5*(thl(2)+thl(5)); %CL_inv=1/CL
L_R=CR_inv/CL_inv;
CR=1/CR_inv;
CL=1/CL_inv;

vcek=transpose(transpose(thl)*transpose(zeta)); %vrek is vrhat 
vrek=u1+u2-vcek;

figure;
subplot(2,2,1)
plot(t,vr);
title("$V_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

subplot(2,2,2)
plot(t,vrek);
title("$\hat{V}_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

subplot(2,2,[3 4])
plot(t,vr,'r--',t,vrek,'g*');
title("Real measured $V_{R}(t)$ and Model $\hat{V}_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
legend({'$V_{R}(t)$','$\hat{V}_{R}(t)$'},'Interpreter','latex')
grid on

figure;
evr=vr-vrek;
plot(t,evr);
title("Prediction Error $e_{R}(t)$=$V_{R}(t)$-$\hat{V}_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

figure;
subplot(2,2,1)
plot(t,vc);
title("$V_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

subplot(2,2,2)
plot(t,vcek);
title("$\hat{V}_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)");
ylabel("Volt (V)");
grid on

subplot(2,2,[3 4])
plot(t,vc,'r--',t,vcek,'g*');
title("Real measured $V_{C}(t)$ and Model $\hat{V}_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
legend({'$V_{C}(t)$','$\hat{V}_{C}(t)$'},'Interpreter','latex')
grid on

figure;
evc=vc-vcek;
plot(t,evc);
title("Prediction Error $e_{C}(t)$=$V_{C}(t)$-$\hat{V}_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on



