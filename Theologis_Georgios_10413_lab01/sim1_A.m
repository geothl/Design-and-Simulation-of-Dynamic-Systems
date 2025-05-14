function xdot = sim1_A(t,y)
xdot=[y(2);-2/8.5*y(1)-0.65/8.5*y(2)+(10*cos(0.5*pi*t)+3)/8.5];
end