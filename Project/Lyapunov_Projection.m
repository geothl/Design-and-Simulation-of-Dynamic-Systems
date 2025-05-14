function dynamics = Lyapunov_Projection(dyn,a,b,thm,k,rho,Gamma,xd)
j=(dyn(1)-xd)/rho;
epsilon=log((1+j)/(1-j));
u=-k/dyn(4)*epsilon;
dynamics(1)=a*dyn(1)+b*u;
ex=dyn(1)-dyn(2);
dynamics(2)=dyn(3)*dyn(1)+dyn(4)*u+thm*ex;
adot_ek=ex*dyn(1);
bdot_ek=ex*u;
para_dot_ek=[adot_ek;bdot_ek];
dynamics(3:4)=Projection(dyn(3),dyn(4),para_dot_ek,Gamma);
dynamics=transpose(dynamics);
end
