clear;
t=transpose(0:0.01:10);
N=length(t);
M=5;%certain number of time moments where the error of measuremnt is huge
arrofmoments=0;
counter=0;
hr=zeros(N,1);
hc=zeros(N,1);
while(counter<M)
r = randi([1 N]);
if(~ismember(r, arrofmoments))
counter=counter+1;
hr(r)=normrnd(500,50);
hc(r)=normrnd(500,50);
arrofmoments(counter)=r;
end
end

% % %Alternative inmplementation of time moment choice for errors. If to be
% % %tested the upper while loop must be in comments and this while loop must
% % %be removed from comments.
% pe=0.005;
% while(counter<N)
% r=rand();
% counter=counter+1;
% if(r<pe)
% hr(counter)=normrnd(500,50);
% hc(counter)=normrnd(500,50);
% end
% end



u1=3*sin(pi*t);
u2=2.5*ones(length(t),1);
vr=zeros(length(t));
vc=zeros(length(t));
[vr,vc]=v(t);
vr_real=vr;
vc_real=vc;
vr=vr+hr;
vc=vc+hc;
l1=3;
l2=2;
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
title("$V_{R}(t) with error$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

subplot(2,2,2)
plot(t,vr_real,'r--',t,vrek,'g*');
title("Measured $V_{R}(t) without error$ and Model $\hat{V}_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
legend({'$V_{R}(t) no error$','$\hat{V}_{R}(t)$'},'Interpreter','latex')
grid on

subplot(2,2,[3 4])
plot(t,vr,'r--',t,vrek,'g*');
title("Measured $V_{R}(t) with error$ and Model $\hat{V}_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
legend({'$V_{R}(t)$','$\hat{V}_{R}(t)$'},'Interpreter','latex')
grid on

figure;
subplot(1,2,1)
evr=vr-vrek;
plot(t,evr);
title("Prediction Error $e_{R}(t)$=$V_{R}(t) $with error-$\hat{V}_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

subplot(1,2,2)
evr=vr_real-vrek;
plot(t,evr);
title("Prediction Error $e_{R}(t)$=$V_{R}(t)$without error-$\hat{V}_{R}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

figure;
subplot(2,2,1)
plot(t,vc);
title("$V_{C}(t) with error $",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

subplot(2,2,2)
plot(t,vc_real,'r--',t,vcek,'g*');
title("Real measured $V_{C}(t)$ without error and Model $\hat{V}_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
legend({'$V_{C}(t) no error$','$\hat{V}_{C}(t)$'},'Interpreter','latex')
grid on

subplot(2,2,[3 4])
plot(t,vc,'r--',t,vcek,'g*');
title("Real measured $V_{C}(t)$ with error and Model $\hat{V}_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
legend({'$V_{C}(t)$','$\hat{V}_{C}(t)$'},'Interpreter','latex')
grid on

figure;
subplot(1,2,1)
evc=vc-vcek;
plot(t,evc);
title("Prediction Error $e_{C}(t)$=$V_{C}(t)$ with error-$\hat{V}_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on

subplot(1,2,2)
evc=vc_real-vcek;
plot(t,evc);
title("Prediction Error $e_{C}(t)$=$V_{C}(t)$ without error-$\hat{V}_{C}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Volt (V)");
grid on



