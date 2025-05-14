clear;
%%%%Thema 1
tspan=transpose(0:0.01:20);
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);%better ode simulation
a=4;
b=1.5;
Gamma=[5; 4;];
thm=1;
[t,indentifiers] = ode45(@(t,dyn)  gradient_method_indentifier(dyn,5,a,b,Gamma,thm),tspan,[0 0 0 0 0],opts);
x=indentifiers(:,1);
ph1=indentifiers(:,2);
ph2=indentifiers(:,3);
th1=indentifiers(:,4);
th2=indentifiers(:,5);
ahat=thm*ones(length(t),1)-th1;
bhat=th2;
xhat=th1.*ph1+th2.*ph2;
ex=x-xhat;
figure;
plot(t,x,t,xhat);
title("$x(t)$ and $\hat{x}(t)$ of system with $u=5$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x(t)$','$\hat{x}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,ex);
title("$e_{x}=x(t)-\hat{x}(t)$ of system with $u=5$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,ahat,t,bhat,tect,a*ones(length(tect),1),'r.',tect,b*ones(length(tect),1),'.');
title("$\hat{a}(t)$ and $\hat{b}(t)$ of system with $u=5$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of a,b parameters");
legend({'$\hat{a}(t)$','$\hat{b}(t)$'},'Interpreter','latex')
grid on;


thm=1;
Gamma=[12; 14;];
[t,indentifiers] = ode45(@(t,dyn)  gradient_method_indentifier(dyn,5*sin(2*t),a,b,Gamma,thm),tspan,[0 0 0 0 0],opts);
x=indentifiers(:,1);
ph1=indentifiers(:,2);
ph2=indentifiers(:,3);
th1=indentifiers(:,4);
th2=indentifiers(:,5);
ahat=thm*ones(length(t),1)-th1;
bhat=th2;
xhat=th1.*ph1+th2.*ph2;
ex=x-xhat;
figure;
plot(t,x,t,xhat);
title("$x(t)$ and $\hat{x}(t)$ of system with $u=5\sin(2t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x(t)$','$\hat{x}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,ex);
title("$e_{x}=x(t)-\hat{x}(t)$ of system with $u=5\sin(2t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,ahat,t,bhat,tect,a*ones(length(tect),1),'r.',tect,b*ones(length(tect),1),'.');
title("$\hat{a}(t)$ and $\hat{b}(t)$ of system with $u=5\sin(2\pi t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of a,b parameters");
legend({'$\hat{a}(t)$','$\hat{b}(t)$'},'Interpreter','latex')
grid on;




%%%%%%%%THEMA 2
tspan=transpose(0:0.01:20);
a=2;
b=5;
Gamma=[2; 2;];
h0=0.5;
f=40;
[t,indentifiers] = ode45(@(t,dyn) Lyapunov_parallel_method_modelopoihsh(dyn,5*sin(2*t),a,b,Gamma,h0*sin(2*pi*f*t)),tspan,[0 0 0 0],opts);
x=indentifiers(:,1);
xhat=indentifiers(:,2);
ahat=indentifiers(:,3);
bhat=indentifiers(:,4);
ex=x-xhat;
exmeasured=x+h0*sin(2*pi*f*t)-xhat;
figure;
plot(t,x,t,xhat);
title(['$x(t)$,$\hat{x}(t)$ on Parallel Lyapunov Method with measurement noise h0=',num2str(h0),' f=',num2str(f),' Hz'],'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x(t)$','$\hat{x}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,exmeasured);
title(['$e_{x_{measured}}=x(t)+\eta(t)-\hat{x}(t)$ on Parallel Lyapunov Method with measurement noise h0=',num2str(h0),' f=',num2str(f),' Hz'],'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
plot(t,ex);
title(['$e_{x}=x(t)-\hat{x}(t)$ on Parallel Lyapunov Method with measurement noise h0=',num2str(h0),' f=',num2str(f),' Hz'],'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,ahat,t,bhat,tect,a*ones(length(tect),1),'r.',tect,b*ones(length(tect),1),'.');
title(['$\hat{a}(t)$ and $\hat{b}(t)$ on Parallel Lyapunov Method with measurement noise h0=',num2str(h0),' f=',num2str(f),' Hz'],'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of a,b parameters");
legend({'$\hat{a}(t)$','$\hat{b}(t)$'},'Interpreter','latex')
grid on;




a=2;
b=5;
Gamma=[2; 2;];
h0=0.5;
f=40;
thm=3;
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);%better ode simulation
[t,indentifiers] = ode45(@(t,dyn) Lyapunov_seriesparallel_method_modelopoihsh(dyn,5*sin(2*t),a,b,Gamma,h0*sin(2*pi*f*t),thm),tspan,[0 0 0 0],opts);
x=indentifiers(:,1);
xhat=indentifiers(:,2);
ahat=indentifiers(:,3);
bhat=indentifiers(:,4);
ex=x+h0-xhat;
exmeasured=x-xhat;
figure;
plot(t,x,t,xhat);
title("$x(t)$ and $\hat{x}(t)$ on Series-Parallel Lyapunov Method with measurement noise h0=0.5 f=40 Hz",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x(t)$','$\hat{x}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,ex);
title("$e_{x}=x(t)-\hat{x}(t)$ on Series-Parallel Lyapunov Method with measurement noise h0=0.5 f=40 Hz",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
plot(t,exmeasured);
title("$e_{x}=x_{measured}(t)-\hat{x}(t)$ on Series-Parallel Lyapunov Method with measurement noise h0=0.5 f=40 Hz",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,ahat,t,bhat,tect,a*ones(length(tect),1),'r.',tect,b*ones(length(tect),1),'.');
title("$\hat{a}(t)$ and $\hat{b}(t)$ on Series-Parallel Lyapunov Method with measurement noise h0=0.5 f=40 Hz",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of a,b parameters");
legend({'$\hat{a}(t)$','$\hat{b}(t)$'},'Interpreter','latex')
grid on;



%%%%%%%%%%%%THEMA 3




A=[-1 1;-4 0;];
B=[2;1;];
Gamma=[2; 2;];
[t,indentifiers] = ode45(@(t,dyn) Lyapunov_parallel_method_modelopoihsh_2nd_dimension(dyn,4*sin(pi*t)+2*sin(8*pi*t),A,B,Gamma),tspan,zeros(1,10),opts);
x1=indentifiers(:,1);
x2=indentifiers(:,2);
x1hat=indentifiers(:,3);
x2hat=indentifiers(:,4);
a11hat=indentifiers(:,5);
a12hat=indentifiers(:,6);
a21hat=indentifiers(:,7);
a22hat=indentifiers(:,8);
b1hat=indentifiers(:,9);
b2hat=indentifiers(:,10);
ex1=x1-x1hat;
ex2=x2-x2hat;

figure;
plot(t,x1,t,x1hat);
title("$x_{1}(t)$ and $\hat{x}_{1}(t)$ on Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x_{1}(t)$','$\hat{x}_{1}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,x2,t,x2hat);
title("$x_{2}(t)$ and $\hat{x}_{2}(t)$ on Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x_{2}(t)$','$\hat{x}_{2}(t)$'},'Interpreter','latex')
grid on;

figure;
plot(t,ex1,t,ex2);
title("$e_{x}=x(t)-\hat{x}(t)$ on Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
legend({'$e_{x1}(t)$','$e_{x2}(t)$'},'Interpreter','latex')
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,a11hat,t,a12hat,t,a21hat,t,a22hat);
title("$\hat{A}(t)$ on Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of A parameters");
legend({'$\hat{a}_{11}(t)$','$\hat{a}_{12}(t)$','$\hat{a}_{21}(t)$','$\hat{a}_{22}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,b1hat,t,b2hat);
title("$\hat{B}(t)$ on Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of B parameters");
legend({'$\hat{b}_{11}(t)$','$\hat{b}_{12}(t)$'},'Interpreter','latex')
grid on;



tspan=transpose(0:0.01:150);
A=[-1 1;-4 0;];
B=[2;1;];
Gamma=[2;2;];
Thm=12*[1 0;0 1];
[t,indentifiers] = ode45(@(t,dyn) Lyapunov_seriesparallel_method_modelopoihsh_2nd_dimension(dyn,4*sin(pi*t)+2*sin(8*pi*t),A,B,Gamma,Thm),tspan,zeros(1,10),opts);
x1=indentifiers(:,1);
x2=indentifiers(:,2);
x1hat=indentifiers(:,3);
x2hat=indentifiers(:,4);
a11hat=indentifiers(:,5);
a12hat=indentifiers(:,6);
a21hat=indentifiers(:,7);
a22hat=indentifiers(:,8);
b1hat=indentifiers(:,9);
b2hat=indentifiers(:,10);
ex1=x1-x1hat;
ex2=x2-x2hat;

figure;
plot(t,x1,t,x1hat);
title("$x_{1}(t)$ and $\hat{x}_{1}(t)$ on Series-Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x_{1}(t)$','$\hat{x}_{1}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,x2,t,x2hat);
title("$x_{2}(t)$ and $\hat{x}_{2}(t)$ on Series-Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x_{2}(t)$','$\hat{x}_{2}(t)$'},'Interpreter','latex')
grid on;

figure;
plot(t,ex1,t,ex2);
title("$e_{x}=x(t)-\hat{x}(t)$ on Series-Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
legend({'$e_{x1}(t)$','$e_{x2}(t)$'},'Interpreter','latex')
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,a11hat,t,a12hat,t,a21hat,t,a22hat);
title("$\hat{A}(t)$ on Series-Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of A parameters");
legend({'$\hat{a}_{11}(t)$','$\hat{a}_{12}(t)$','$\hat{a}_{21}(t)$','$\hat{a}_{22}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,b1hat,t,b2hat);
title("$\hat{B}(t)$ on Series-Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of B parameters");
legend({'$\hat{b}_{11}(t)$','$\hat{b}_{12}(t)$'},'Interpreter','latex')
grid on;


%%%%%%THEMA 4

tspan=transpose(0:0.01:100);

% 
f1=@(x) 0.5*sin(x)*x;
f2=@(x) -0.25*x^2;
th1=0.5;
th2=2;
Gamma=1*[10; 9;]; 
thm=2;

[t,indentifiers] = ode45(@(t,dyn) Lyapunov_seriesparallel_method_modelopoihsh_nonlinear(dyn,1.5*sin(2*pi*t)*exp(-3*t) ,th1,th2,Gamma,thm,f1),tspan,[0 0 0 0],opts);

x=indentifiers(:,1);
xhat=indentifiers(:,2);
th1hat=indentifiers(:,3);
th2hat=indentifiers(:,4);
ex=x-xhat;
figure;
plot(t,x,t,xhat);
title("$x(t)$ and $\hat{x}(t)$ on Series-Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x(t)$','$\hat{x}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,ex);
title("$e_{x}=x(t)-\hat{x}(t)$ on Series-Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,th1hat,t,th2hat,tect,th1*ones(length(tect),1),'r.',tect,th2*ones(length(tect),1),'.');
title("$\hat{\theta}_{1}(t)$ and $\hat{\theta}_{2}(t)$ on Series-Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of th1,th2 parameters");
legend({'$\hat{\theta}_{1}(t)$','$\hat{\theta}_{2}(t)$'},'Interpreter','latex')
grid on;


tspan=transpose(0:0.01:17);
th1=0.5;
th2=2;
Gamma=1*[10; 9;]; 
thm=2;
[t,indentifiers] = ode45(@(t,dyn2) Lyapunov_seriesparallel_method_modelopoihsh_nonlinear(dyn2,1.5*sin(2*pi*t)*exp(-3*t) ,th1,th2,Gamma,thm,f2),tspan,[0 0 0 0],opts);

x=indentifiers(:,1);
xhat=indentifiers(:,2);
th1hat=indentifiers(:,3);
th2hat=indentifiers(:,4);
ex=x-xhat;
figure;
plot(t,x,t,xhat);
title("$x(t)$ and $\hat{x}(t)$ on Series-Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x(t)$','$\hat{x}(t)$'},'Interpreter','latex')
grid on;
figure;
plot(t,ex);
title("$e_{x}=x(t)-\hat{x}(t)$ on Series-Parallel Lyapunov Method ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Error");
grid on;
figure;
tect=linspace(0,t(end),40);
plot(t,th1hat,t,th2hat,tect,th1*ones(length(tect),1),'r.',tect,th2*ones(length(tect),1),'.');
title("$\hat{\theta}_{1}(t)$ and $\hat{\theta}_{2}(t)$ on Series-Parallel Lyapunov Method",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Estimation of th1,th2 parameters");
legend({'$\hat{\theta}_{1}(t)$','$\hat{\theta}_{2}(t)$'},'Interpreter','latex')
grid on;

tspan=transpose(0:0.01:20.66);
[t,indentifiers] = ode45(@(t,dyn2) Lyapunov_seriesparallel_method_modelopoihsh_nonlinear(dyn2,1.5*sin(2*pi*t)*exp(-3*t) ,th1,th2,Gamma,thm,f2),tspan,[0 0 0 0],opts);
x=indentifiers(:,1);
figure;
plot(t,x)
ylabel('x(t)')
xlabel('t')
title('x(t) state')
