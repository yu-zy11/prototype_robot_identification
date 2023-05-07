%function：get minimal set of iertial parameters by  QR decomposition
%reference:https://zhuanlan.zhihu.com/p/549740247
close all
clear all
clc
n=20;
q1set=rand(n,1);q2set=rand(n,1);q3set=rand(n,1);
dq1set=rand(n,1);dq2set=rand(n,1);dq3set=rand(n,1);
ddq1set=rand(n,1);ddq2set=rand(n,1);ddq3set=rand(n,1);
Y=[]
for i=1:1:n
    q=[q1set(i),q2set(i),q3set(i)];
    dq=[dq1set(i),dq2set(i),dq3set(i)];
    ddq=[ddq1set(i),ddq2set(i),ddq3set(i)];
    Yi=Ymatrix_example(q,dq,ddq)
    Y=[Y;Yi];
    
end
%% ------step1 
[Q,R]=qr(Y);
% num_base_param=rank(R);
base_param_sequence=[];
none_base_param_sequence=[];
Zero=0.0000001;
%%  ------step2
%handle the influence of computational accuracy
for i=1:1:size(R,2)
    if(abs(R(i,i))<Zero)
        none_base_param_sequence=[none_base_param_sequence,i];
    else
        base_param_sequence=[base_param_sequence,i];
    end
end

Y1=[];
for i=1:1:size(base_param_sequence,2)
    num=base_param_sequence(i);
    Y1=[Y1,Y(:,num)];
end
Y2=[];
for i=1:1:size(none_base_param_sequence,2)
  num=none_base_param_sequence(i);
  Y2=[Y2, Y(:,num)];
end
base_param_num=size(Y1,2);
%% -------------step3
[Q,R]=qr([Y1 Y2]);
n_base_num=size(none_base_param_sequence,2)
% Q1=Q(1:base_param_num,1:base_param_num);
% Q2=Q(1:base_param_num,base_param_num+1:base_param_num+n_base_num);
R1=R(1:base_param_num,1:base_param_num);
R2=R(1:base_param_num,base_param_num+1:end);
beta=inv(R1)*R2;
%handle the influence of computational accuracy
for i=1:1:size(beta,1)
    for j=1:1:size(beta,2)
     if(abs(beta(i,j))<Zero)
         beta(i,j)=0;
     end
    end
end
test=Y2-Y1*beta

% tau=Y*P;
%% ------step4
%符号运算，计算基惯性參數，最小惯性參數集，不可辨识參數。
m = sym('m',[3 1],'real');
mx = sym('mx',[3 1],'real');
my = sym('my',[3 1],'real');
mz = sym('mz',[3 1],'real');
Ixx = sym('Ixx',[3 1],'real');
Ixy = sym('Ixy',[3 1],'real');
Ixz = sym('Ixz',[3 1],'real');
Iyz = sym('Iyz',[3 1],'real');
Iyy = sym('Iyy',[3 1],'real');
Izz = sym('Izz',[3 1],'real');
Param_all=[]
for i=1:1:size(m,1)
    Param_all=[Param_all,[Ixx(i) Ixy(i) Ixz(i) Iyy(i) Iyz(i) Izz(i) mx(i) my(i) mz(i) m(i)]];
end
base_param=[];
for i=1:1:size(base_param_sequence,2)
    n=base_param_sequence(i);
    base_param=[base_param Param_all(n)];
end
none_base_param=[];
for i=1:1:size(none_base_param_sequence,2)
    n=none_base_param_sequence(i);
    none_base_param=[none_base_param,Param_all(n)];
end
%
bata_param=beta*none_base_param';
sequence_not_identified=[];
for i=1:1:size(beta,2)
 if( sum(beta(:,i).^2)==0)
     sequence_not_identified=[sequence_not_identified,none_base_param_sequence(i)];
 end
end
param_not_identified=[];
for i=1:1:size(sequence_not_identified,2)
    n=sequence_not_identified(i);
    param_not_identified=[param_not_identified Param_all(n)];
end
%% 
base_parameter=base_param
min_set_parameter=(base_param'+beta*none_base_param')'
paramter_not_identified=param_not_identified
