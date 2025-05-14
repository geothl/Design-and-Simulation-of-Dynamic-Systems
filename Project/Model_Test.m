function dynamics=Model_Test(dyn,am,bm,k,xd,rho)
x=dyn(1);
xm=dyn(2);
j=(x-xd)/rho;
epsilon=log((1+j)/(1-j));
a=0.75;
b=1.25;
u=-k/bm*epsilon;
xdot=a*x+b*u;

jm=(xm-xd)/rho;
epsilonm=log((1+jm)/(1-jm));
um=-k/bm*epsilonm;
xmdot=am*xm+bm*um;
dynamics=[xdot;xmdot;];
end
