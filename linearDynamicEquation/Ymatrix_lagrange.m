% close all
% clear all;
% clc;
m = sym('m',[3 1],'real');
hx = sym('hx',[3 1],'real');
hy = sym('hy',[3 1],'real');
hz = sym('hz',[3 1],'real');
Ixx = sym('Ixx',[3 1],'real');
Ixy = sym('Ixy',[3 1],'real');
Ixz = sym('Ixz',[3 1],'real');
Iyz = sym('Iyz',[3 1],'real');
Iyy = sym('Iyy',[3 1],'real');
Izz = sym('Izz',[3 1],'real');

syms q1 q2 q3 real;
syms dq1 dq2 dq3 real;
syms ddq1 ddq2 ddq3 real;
syms l1 l2 l3  g real;
G=[0 0 -g]'; %gravity
%dynamic parameters
Pi_1=[Ixx(1),Ixy(1),Ixz(1),Iyy(1),Iyz(1),Izz(1),m(1)*hx(1),m(1)*hy(1),m(1)*hz(1),m(1)]';
Pi_2=[Ixx(2),Ixy(2),Ixz(2),Iyy(2),Iyz(2),Izz(2),m(2)*hx(2),m(2)*hy(2),m(2)*hz(2),m(2)]';
Pi_3=[Ixx(3),Ixy(3),Ixz(3),Iyy(3),Iyz(3),Izz(3),m(3)*hx(3),m(3)*hy(3),m(3)*hz(3),m(3)]';
PI=[Pi_1', Pi_2', Pi_3']';

q = [q1;q2;q3];
dq = [dq1;dq2;dq3];
ddq = [ddq1;ddq2;ddq3];
%% ----------import robot model
robot_3axis_model;

%% ----------step1:get HTM tree for joint postion of [q1 q2 q3]
T01=T01_init*expMatrix(axis1,q1);
T12=T12_init*expMatrix(axis2,q2);
T23=T23_init*expMatrix(axis3,q3);
T02=T01*T12;
T03=T02*T23;
R01=T01(1:3,1:3);P01=T01(1:3,4);
R12=T12(1:3,1:3);P12=T12(1:3,4);
R23=T23(1:3,1:3);P23=T23(1:3,4);
R02 = T02(1:3,1:3);P02=T02(1:3,4);
R03 = T03(1:3,1:3);P03=T03(1:3,4);

%% ----------step2:calculate each link's velocity and acceleration in body frame
%get angular velocity
omega_base=[0 0 0]';
omega0=omega_base;
omega1=R01'*omega0+axis1(1:3)*dq1; 
omega2=R12'*omega1+axis2(1:3)*dq2;
omega3=R23'*omega2+axis3(1:3)*dq3;
%get linear velocity
vel1_in_base = jacobian(P01, q) * dq;
vel2_in_base = jacobian(P02, q) * dq;
vel3_in_base = jacobian(P03, q) * dq;
vel1=simplify(R01'*vel1_in_base);
vel2=simplify(R02'*vel2_in_base);
vel3=simplify(R03'*vel3_in_base);

%% --------- kinetic energy and potential energy
T_1 = getTi(omega1, vel1);
T_2 = getTi(omega2, vel2);
T_3 = getTi(omega3, vel3);
T = [T_1, T_2, T_3];%kinetic energy=T*PI

V_1 = getVi(G, R01, P01);
V_2 = getVi(G, R02, P02);
V_3 = getVi(G, R03, P03);
V = [V_1, V_2, V_3];%potential energy=V*PI

%lagrange
L=T-V;
T_vq = jacobian(T, dq);
Y1 = diff(T_vq, q(1))*dq(1) + diff(T_vq, q(2))*dq(2) + diff(T_vq, q(3))*dq(3)...
    + diff(T_vq, dq(1))*ddq(1) + diff(T_vq, dq(2))*ddq(2) + diff(T_vq, dq(3))*ddq(3);
Y2 = jacobian(T, q);
Y3 = jacobian(V, q);
Y = (Y1 - Y2 + Y3)';
Y_tmp = simplify(Y)
% 
% tau = Y*PI;
% tau = simplify(tau)
% tau3=tau(3,:)
%  H = simplify(jacobian(tau, ddq));
%  C = simplify(tau - H*ddq);
% save('RegressorMatrix','Y');
function T=expMatrix(rotation_axis,angle)
 w=rotation_axis(1:3);
 W=Vec2so3(w);
 T(1:3,1:3)=eye(3)+sin(angle)*W+(1-cos(angle))*W^2;
 T(1:3,4)=[0 0 0]';
 T(4,1:4)=[0 0 0 1];
end
function screw=Vec2so3(vec)
[row,col]=size(vec);
if row~=3 && col~=1
    error("worng,input of Vec2so3 must be 3*1 vector");
end
screw=[0     -vec(3)  vec(2);
       vec(3)   0    -vec(1);
      -vec(2)   vec(1)  0];
end

function T_hat = getTi(w,v)
K = [w(1) w(2) w(3)  0    0    0
      0   w(1)  0   w(2) w(3)  0
      0    0   w(1)  0   w(2) w(3)];
T_hat = [1/2.* w' * K,   v' * Vec2so3(w),  1/2.* v' * v];
end

%R代表从0到i的旋转矩阵,P代表关节i在世界坐标系下的位置
function V_hat = getVi(G,R,P)
V_hat = [0,0,0,0,0,0, -G' * R,  -G'*P];

end