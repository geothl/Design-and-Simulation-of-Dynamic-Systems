function dynamics = gradient_method_indentifier(dyn,u,a,b,Gamma,thm)
dynamics(1)=-a*dyn(1)+b*u;
dynamics(2)=-thm*dyn(2)+dyn(1);
dynamics(3)=-thm*dyn(3)+u;
ex=dyn(1)-dyn(4)*dyn(2)-dyn(5)*dyn(3);
dynamics(4)=Gamma(1)*dyn(2)*ex;
dynamics(5)=Gamma(2)*dyn(3)*ex;

dynamics=transpose(dynamics);
end