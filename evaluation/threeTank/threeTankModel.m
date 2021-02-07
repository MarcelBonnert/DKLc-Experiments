function dx = threeTankModel(t,x,u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

A = 0.14^2*pi/4;
q1 = 2.222870181488818e-05;
q2 = 4.245214693023062e-05;
q3 = 2.217589250364239e-05;
g = 9.81;

dx = [  (u(1) - q1*sign(x(1)-x(2))*sqrt(2*g*abs(x(1)-x(2))))/A;
        (q1*sign(x(1)-x(2))*sqrt(2*g*abs(x(1)-x(2))) + q3*sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2))) - q2*sqrt(2*g*abs(x(2))))/A;
        (u(2) - q3*sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2))))/A];

end

