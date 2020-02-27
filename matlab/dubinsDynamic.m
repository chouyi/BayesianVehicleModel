function [next_x]=dubinsDynamic(theta,x,h)

    u_v=theta(1); 
    u_a=theta(2); 
    next_x(1,1)=x(1,1)+h*x(3,1)*cos(x(4,1))+(h^2/2)*u_v*cos(x(4,1));
    next_x(2,1)=x(2,1)+h*x(3,1)*sin(x(4,1))+(h^2/2)*u_v*sin(x(4,1));
    next_x(3,1)=x(3,1)+h*u_v;
    next_x(4,1)=x(4,1)+h*u_a;


end