function plot_Thema_2(tspan,A,omega,k,lamda,thm,Gamma)
x0=0;
xd=@(m) A*cos(omega*m);
rho_0=2*abs(x0-xd(0))+1;
rho_inf=0.05;
rho=(rho_0-rho_inf)*exp(-lamda*tspan)+rho_inf;
[t,indentifiers] = Model_Simulation(A,k,thm,Gamma,lamda,omega,tspan,0,0,0);
x=indentifiers(:,1);
x_ek=indentifiers(:,2);
a_ek=indentifiers(:,3);
b_ek=indentifiers(:,4);
ex=x-x_ek;
xd_arr=xd(t);
j=(x-xd_arr)./transpose(rho);
epsilon=log((1+j)/(1-j));
u=-k*epsilon./b_ek;

figure;
plot(t,x,t,x_ek);
title("State $x(t)$ of real system and State $\hat{x}(t)$ of indentification system",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("State");
legend({'$x(t)$','$\hat{x}(t)$'},'Interpreter','latex')
grid on;


figure;
plot(t,a_ek,t,b_ek,t,0.75*ones(length(t)),'--',t,1.25*ones(length(t)),'--');
title("Estimated Parameters $\hat{a}(t)$ and $\hat{\beta}(t)$ ",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Parameter units");
legend({'$\hat{a}(t)$ ','$\hat{\beta}(t)$'},'Interpreter','latex')
grid on;


figure;
plot(t,ex);
title("Identification Error $e_{x}(t)=x(t)-\hat{x}(t)$",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Ιdentification Error");
grid on;


figure;
plot(t,x-xd_arr,t,rho,t,-rho);
title("Observation error $e_{x_{d}}(t)=x(t)-x_{d}(t)$ with its upper $\rho(t)$ and lower $-\rho(t)$ bound",'Interpreter','latex')
xlabel("Time in seconds (s)")
ylabel("Observation error");
legend({'$e_{x_{d}}(t)=x(t)-x_{d}(t)$',' $\rho(t)$',' $-\rho(t)$'},'Interpreter','latex')
grid on;

figure;
plot(t,u);
title("Control Signal u(t)")
xlabel("Time in seconds (s)")
ylabel("u(t)");
grid on;

end

