function [l2,am,bm]=Cross_Validation_Step(Training_Data,Test_Data,M,k,thm,Gamma,lamda,N,tspan)

elements_per_section=N/M;
index=1;
for h=1:elements_per_section:(N-elements_per_section)
    
for i=0:1:(elements_per_section-1)
    
A=Training_Data(1,h+i);   
omega=Training_Data(2,h+i);
B=Training_Data(3,h+i);   
C=Training_Data(4,h+i);   
D=Training_Data(5,h+i);   

[t,indentifiers] = Model_Simulation(A,k,thm,Gamma,lamda,omega,tspan,B,C,D);
a_sec_arr(i+1)=indentifiers(end,3);
b_sec_arr(i+1)=indentifiers(end,4);
end
a_arr(index)=mean(a_sec_arr);
b_arr(index)=mean(b_sec_arr);
index=index+1;
end

am=mean(a_arr);
bm=mean(b_arr);

for i=1:1:elements_per_section
A=Test_Data(1);   
omega=Test_Data(2);
B=Test_Data(3);   
C=Test_Data(4);
D=Test_Data(5); 

xd=@(m) A*cos(omega*m)*exp(-D*m)+B+C*exp(-(m-2)^2/2);
x0=0;
rho_0=2*abs(x0-xd(0))+1;
rho_inf=0.05;

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);%better ode simulation
[t,indentifiers] = ode45(@(t,dyn) Model_Test(dyn,am,bm,k,xd(t),(rho_0-rho_inf)*exp(-lamda*t)+rho_inf) ,tspan,[x0 x0],opts);
x=indentifiers(:,1);
xm=indentifiers(:,2);
em=x-xm;
l2=L2_norm_calculator(em,t);
l2_arr(i)=l2;
end
l2=mean(l2_arr);
end