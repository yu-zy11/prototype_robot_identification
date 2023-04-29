%build a robot mdoel for a simple robot of 3-axis,refer to robot_3axis.jpg
syms l1 l2 l3 real;%link length
%inital homegenious transformation matrix T_init
T01_init=eye(4);
T12_init=[1 0 0 0;0 1 0 0;0 0 1 -l1;0 0 0 1];
T23_init=[1 0 0 0;0 1 0 0;0 0 1 -l2;0 0 0 1];
T34_init=[1 0 0 0;0 1 0 0;0 0 1 -l3;0 0 0 1];
%rotation axis in body frame
axis1=[0 1 0]';
axis2=[1 0 0]';
axis3=[1 0 0]';