%% 

clc
clear

syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 = [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1];
T12 = [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1];
T23 = [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
T34 = [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1];
T45 = [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1];
T56 = [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1];

T06=simplify(T01*T12*T23*T34*T45*T56);
r11=T06(1,1);
r12=T06(1,2);
r13=T06(1,3);
r14=T06(1,4);
r21=T06(2,1);
r22=T06(2,2);
r23=T06(2,3);
r24=T06(2,4);
r31=T06(3,1);
r32=T06(3,2);
r33=T06(3,3);
r34=T06(3,4);
r41=T06(4,1);
r42=T06(4,2);
r43=T06(4,3);
r44=T06(4,4);


Theta1=0.1 ;
Theta2=0.1 ;
Theta3=0.1 ;
Theta4=0.1 ;
Theta5=0.1 ;
Theta6=0.1 ;


%% Z-Y-Z euler angles
beta=atan2(sqrt(r31^2+r32^2),r33)
alpha=atan2(r23/sin(beta),r13/sin(beta))
gamma=atan2(r32/sin(beta),-r31/sin(beta))
eval(alpha)
eval(beta)
eval(gamma)
%% X-Y-Z fixed angles
beta_fixed=beta;
alpha_fixed=alpha;
gamma_fixed=gamma;
% equivalent Angle-Axis
theta_equivalent_angle=acos((r11+r22+r33-1)/2);
k_equivalent_axis=(1/(2*sin(theta_equivalent_angle)))*[r32-r23 ; r13-r31 ; r21-r12];
% Euler parameters
epsilon4=0.5*sqrt(1+r11+r22+r33);
epsilon1=(r32-r23)/(4*epsilon4);
epsilon2=(r13-r31)/(4*epsilon4);
epsilon3=(r21-r12)/(4*epsilon4);