%brief: linearize dynamic equation of motion for 3-axis robot using
%Newton-Euler method.reference:https://zhuanlan.zhihu.com/p/109770180
close all
clear all
clc
%% ----------import robot model
%including initial homegenious transformation matrixs(HTM) and rotation directions
robot_3axis_model;

%% ----------step1:get HTM tree for joint postion of [q1 q2 q3]
syms q1 q2 q3 real;% joint positions
T01=T01_init*expMatrix(axis1,q1);
T12=T12_init*expMatrix(axis2,q2);
T23=T23_init*expMatrix(axis3,q3);
R01=T01(1:3,1:3);P01=T01(1:3,4);
R12=T12(1:3,1:3);P12=T12(1:3,4);
R23=T23(1:3,1:3);P23=T23(1:3,4);

%% ----------step2:using Neuton-Euler to calculate each link's velocity and acceleration in body frame
syms dq1 dq2 dq3 ddq1 ddq2 ddq3 real;%symbols for joint velocities and accelerations
syms g real
G=[0 0 -g]'; %gravity
omega_base=[0 0 0]';
omega0=omega_base;
%formulation(3.1),calculate angular velocity for each link
omega1=R01'*omega0+axis1*dq1; 
omega2=R12'*omega1+axis2*dq2;
omega3=R23'*omega2+axis3*dq3;

vel_base=[0 0 0]';
vel0=vel_base;
vel1=R01'*(vel0+Vec2so3(omega0)*P01);
vel2=R12'*(vel1+Vec2so3(omega1)*P12);
vel3=R23'*(vel2+Vec2so3(omega2)*P23);

%formulation(3.2),calculate angular acceleration for each link
ang_acc_base=[0 0 0]';
ang_acc0=ang_acc_base;
ang_acc1=R01'*ang_acc0+Vec2so3(omega1)*axis1*dq1+axis1*ddq1;
ang_acc2=R12'*ang_acc1+Vec2so3(omega2)*axis2*dq2+axis2*ddq2;
ang_acc3=R23'*ang_acc2+Vec2so3(omega3)*axis3*dq3+axis3*ddq3;
%formulation(3.2),calculate acceleration for each link
acc_base=-G;
acc0=acc_base;
acc1=R01'*(acc0+Vec2so3(ang_acc0)*P01+Vec2so3(omega0)*Vec2so3(omega0)*P01);
acc2=R12'*(acc1+Vec2so3(ang_acc1)*P12+Vec2so3(omega1)*Vec2so3(omega1)*P12);
acc3=R23'*(acc2+Vec2so3(ang_acc2)*P23+Vec2so3(omega2)*Vec2so3(omega2)*P23);

%% -----------step3: linearize dynamics equation of motion
%formulation (3.23)(3.24)
H1=getHi(ang_acc1,omega1,acc1);
H2=getHi(ang_acc2,omega2,acc2);
H3=getHi(ang_acc3,omega3,acc3);
%formulation (3.27)
A1=getAi(ang_acc1,omega1,acc1);
A2=getAi(ang_acc2,omega2,acc2);
A3=getAi(ang_acc3,omega3,acc3);
%formulation (3.28),the subscript of Yf in (3.28) is 6, but there is
%3,because the DOF of model here is 3
Yf3=[zeros(3,20) H3];
Yn3=[zeros(3,20) A3];
Yf2=[zeros(3,10) H2 zeros(3,10)]+R23*Yf3;                  %formulation (3.30)
Yn2=[zeros(3,10) A2 zeros(3,10)]+R23*Yn3+Vec2so3(P23)*R23*Yf3; %formulation (3.31)
Yf1=[H1 zeros(3,10) zeros(3,10)]+R12*Yf2;
Yn1=[A1 zeros(3,10) zeros(3,10)]+R12*Yn2+Vec2so3(P12)*R12*Yf2;

%% ----------step4:get coefficient matrix Y
%tau=YP,formulation(3.33)
Y(1,1:30)=axis1(1:3)'*Yn1;
Y(2,1:30)=axis2(1:3)'*Yn2;
Y(3,1:30)=axis3(1:3)'*Yn3;
Y=simplify(Y)

%verification with lagrange 
Ymatrix_lagrange
a=simplify(Y-Y_tmp)

function Ai=getAi(dotwi,wi,ai)
Ai(1:3,1:6)=getK(dotwi)+Vec2so3(wi)*getK(wi);
Ai(1:3,7:9)=-Vec2so3(ai);
Ai(1:3,10)=zeros(3,1);
end
function K=getK(vec3)
K=[vec3(1) vec3(2) vec3(3) 0       0       0;
    0      vec3(1) 0       vec3(2) vec3(3) 0;
    0      0       vec3(1) 0       vec3(2) vec3(3)];
end

function Hi=getHi(dotwi,wi,ai)
Hi(1:3,7:9)=Vec2so3(dotwi)+Vec2so3(wi)*Vec2so3(wi);
Hi(1:3,10)=ai;
Hi(1:3,1:6)=zeros(3,6);
end


function T=expMatrix(rotation_axis,angle)
 w=rotation_axis(1:3);
 W=Vec2so3(w);
 T(1:3,1:3)=eye(3)+sin(angle)*W+(1-cos(angle))*W^2;
 T(1:3,4)=[0 0 0]';
 T(4,1:4)=[0 0 0 1];
end

%brief: Get screw-symmetric matrix for vector 3*1  
function screw=Vec2so3(vec)
[row,col]=size(vec);
if row~=3 && col~=1
    error("worng,input of Vec2so3 must be 3*1 vector");
end
screw=[0     -vec(3)  vec(2);
       vec(3)   0    -vec(1);
      -vec(2)   vec(1)  0];
end


