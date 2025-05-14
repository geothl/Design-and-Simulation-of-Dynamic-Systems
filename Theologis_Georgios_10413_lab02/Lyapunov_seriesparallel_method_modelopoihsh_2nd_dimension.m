function dynamics = Lyapunov_seriesparallel_method_modelopoihsh_2nd_dimension(dyn,u,A,B,Gamma,Thm)
xdot=A*[dyn(1);dyn(2);]+B*u;
dynamics(1)=xdot(1);
dynamics(2)=xdot(2);
A_ek=[dyn(5) dyn(6);dyn(7) dyn(8);];
B_ek=[dyn(9);dyn(10)];
xdot_ek=A_ek*[dyn(1); dyn(2);]+B_ek*u+Thm*([dyn(1);dyn(2);]-[dyn(3);dyn(4);]);
dynamics(3)=xdot_ek(1);
dynamics(4)=xdot_ek(2);
A_ek_dot=Gamma(1)*([dyn(1);dyn(2);]-[dyn(3);dyn(4);])*transpose([dyn(1);dyn(2);]);
B_ek_dot=Gamma(2)*([dyn(1);dyn(2);]-[dyn(3);dyn(4);])*u;
dynamics(5)=A_ek_dot(1,1);
dynamics(6)=A_ek_dot(1,2);
dynamics(7)=A_ek_dot(2,1);
dynamics(8)=A_ek_dot(2,2);

dynamics(9)=B_ek_dot(1);
dynamics(10)=B_ek_dot(2);

dynamics=transpose(dynamics);
end