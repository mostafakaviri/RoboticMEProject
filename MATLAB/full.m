% 0T6 

clc
clear


syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 = [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1];
T12 = [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1];
T23 = [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
T34 = [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1];
T45 = [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1];
T56 = [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1];

T02 = T01 * T12 ;

T03 = T01 * T12 * T23 ;

T04 = T01 * T12 * T23 * T34 ;

T05=  T01 * T12 * T23 * T34 * T45 ;

T06=  T01 * T12 * T23 * T34 * T45 * T56 

%% animation

clc
clear

syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 = [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1];
T12 = [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1];
T23 = [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
T34 = [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1];
T45 = [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1];
T56 = [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1];

T02 = T01 * T12 ;

T03 = T01 * T12 * T23 ;

T04 = T01 * T12 * T23 * T34 ;

T05=  T01 * T12 * T23 * T34 * T45 ;

T06=  T01 * T12 * T23 * T34 * T45 * T56 ;

% mkan A

Ax = T01(1,4)
Ay = T01(2,4)
Az = T01(3,4)

%makan B

Bx = T02(1,4)
By = T02(2,4)
Bz = T02(3,4)

%makan C

Cx = T03(1,4)
Cy = T03(2,4)
Cz = T03(3,4)

%makan D

Dx = T04(1,4)
Dy = T04(2,4)
Dz = T04(3,4)

%makan E

Ex = T05(1,4)
Ey = T05(2,4)
Ez = T05(3,4)

Theta1=0;
Theta4=0;
Theta5=0;
Theta6=0;

%  for Theta1 = -pi : .5 : pi
      for Theta2 = pi/2 : -.1 : -1.31
        for Theta3 = -pi/2 : .1 : 2.93
              
            
            AX=eval(Ax)
            AY=eval(Ay)
            AZ=eval(Az)
            
            A=[AX,AY,AZ]
            
            BX=eval(Bx)
            BY=eval(By)
            BZ=eval(Bz)
            
            B=[BX,BY,BZ]
            
            CX=eval(Cx)
            CY=eval(Cy)
            CZ=eval(Cz)
            
            C=[CX,CY,CZ]
            
            DX=eval(Dx)
            DY=eval(Dy)
            DZ=eval(Dz)
            
            D=[DX,DY,DZ]
            
            EX=eval(Ex)
            EY=eval(Ey)
            EZ=eval(Ez)
            
            E=[EX,EY,EZ]
            
%             plot3(EX,EY,EZ,'.')
            hold off
            
            animation(A,B,C,D,E)
            
       end
     end
% end

%% work space

clc
clear

syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 = [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1];
T12 = [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1];
T23 = [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
T34 = [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1];
T45 = [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1];
T56 = [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1];

T02 = T01 * T12 ;

T03 = T01 * T12 * T23 ;

T04 = T01 * T12 * T23 * T34 ;

T05=  T01 * T12 * T23 * T34 * T45 ;

T06=  T01 * T12 * T23 * T34 * T45 * T56 ;

wrist_position=T04(1:3,4)

i=1 ;

for Theta1=-pi:0.1:pi ;
    for Theta2=-pi:0.1:-0.261 ;
        for Theta3=-pi:0.1:0.35 ;
%             Theta1=0;
            evaluated_wrist_position=eval(wrist_position)
            plot3(evaluated_wrist_position(1),evaluated_wrist_position(2),evaluated_wrist_position(3),'.')
            hold on
            grid on 
        end
    end
end

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

% Z-Y-X euler angles
beta=atan2(-r31,sqrt(r11^2+r21^2));
if cos(beta)==0 && sin(beta)==1
    alpha=0;
    gamma=atan2(r12,r22);
elseif cos(beta)==0 && sin(beta)==-1
     alpha=0;
    gamma=-1*atan2(r12,r22);
else
alpha=atan2(r21/cos(beta),r11/cos(beta));
gamma=atan2(r32/cos(beta),r33/cos(beta));
end
% X-Y-Z fixed angles
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


%% Jacobian ALTERNATIV

syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 =simplify( [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1]) ;
T12 =simplify( [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1]) ;
T23 =simplify( [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]) ;
T34 =simplify( [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1]) ;
T45 =simplify( [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1]) ;
T56 =simplify( [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1]) ;

T02 =simplify( T01 * T12) ;

T03 =simplify( T01 * T12 * T23) ;

T04 =simplify( T01 * T12 * T23 * T34) ;

T05=simplify(  T01 * T12 * T23 * T34 * T45) ;

T06=simplify(  T01 * T12 * T23 * T34 * T45 * T56) ;

Z1=T01(1:3,3);

Z2=T02(1:3,3);

Z3=T03(1:3,3);

Z4=T04(1:3,3);

Z5=T05(1:3,3);

Z6=T06(1:3,3);

O1=T01(1:3,4);

O2=T02(1:3,4);

O3=T03(1:3,4);

O4=T04(1:3,4);

O5=T05(1:3,4);

O6=T06(1:3,4);

J1=simplify([cross(Z1,(O6-O1));Z1]);

J2=simplify([cross(Z2,(O6-O2));Z2]);

J3=simplify([cross(Z3,(O6-O3));Z3]);

J4=simplify([cross(Z4,(O6-O4));Z4]);

J5=simplify([cross(Z5,(O6-O5));Z5]);

J6=simplify([cross(Z6,(O6-O6));Z6]);

JALTER=[simplify(J1) simplify(J2) simplify(J3) simplify(J4) simplify(J5) simplify(J6)]

% Theta1=0.1;
% Theta2=0.1;
% Theta3=0.1;
% Theta4=0.1;
% Theta5=0.1;
% Theta6=0.1;
% 
% eval(JALTER)

%% velocity propagation


syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 =simplify( [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1]) ;
T12 =simplify( [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1]);
T23 =simplify( [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]);
T34 =simplify( [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1]);
T45 =simplify( [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1]);
T56 =simplify( [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1]);

P01=T01(1:3,4);

P12=T12(1:3,4);

P23=T23(1:3,4);

P34=T34(1:3,4);

P45=T45(1:3,4);

P56=T56(1:3,4);

R01=T01(1:3,1:3);

R12=T12(1:3,1:3);

R23=T23(1:3,1:3);

R34=T34(1:3,1:3);

R45=T45(1:3,1:3);

R56=T56(1:3,1:3);

R10=transpose(R01);

R21=transpose(R12);

R32=transpose(R23);

R43=transpose(R34);

R54=transpose(R45);

R65=transpose(R56);


%rotational velocity

syms Thetad1 Thetad2 Thetad3 Thetad4 Thetad5 Thetad6

w00=[0;0;0];

w11=R10*w00+[0;0;Thetad1];

w22=R21*w11+[0;0;Thetad2];

w33=R32*w22+[0;0;Thetad3];

w44=R43*w33+[0;0;Thetad4];

w55=R54*w44+[0;0;Thetad5];

w66=simplify(R65*w55+[0;0;Thetad6]);


%liniar velocity

v00=[0;0;0];

v11=R10*(v00+cross(w00,P01));

v22=R21*(v11+cross(w11,P12));

v33=R32*(v22+cross(w22,P23));

v44=R43*(v33+cross(w33,P34));

v55=R54*(v44+cross(w44,P45));

v66=R65*(v55+cross(w55,P56));

R06=R01 * R12 * R23 * R34 * R45 * R56 ;

v06=simplify( R06 * v66) ; 

w06=simplify(R06 * w66) ;

v06(1,1);
v06(2,1);
v06(3,1);

w06(1,1);
w06(2,1);
w06(3,1);



v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=v06(1,1);
v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=v06(2,1);
v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=v06(3,1);

Jv(1,1)=simplify(v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0));
Jv(2,1)=simplify(v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0));
Jv(3,1)=simplify(v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0));

Jv(1,2)=simplify(v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0));
Jv(2,2)=simplify(v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0));
Jv(3,2)=simplify(v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0));

Jv(1,3)=simplify(v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0));
Jv(2,3)=simplify(v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0));
Jv(3,3)=simplify(v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0));

Jv(1,4)=simplify(v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0));
Jv(2,4)=simplify(v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0));
Jv(3,4)=simplify(v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0));

Jv(1,5)=simplify(v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0));
Jv(2,5)=simplify(v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0));
Jv(3,5)=simplify(v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0));

Jv(1,6)=simplify(v1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1));
Jv(2,6)=simplify(v2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1));
Jv(3,6)=simplify(v3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1));



w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=w06(1,1);
w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=w06(2,1);
w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, Thetad1, Thetad2, Thetad3, Thetad4, Thetad5, Thetad6)=w06(3,1);

Jw(1,1)=simplify(w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0));
Jw(2,1)=simplify(w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0));
Jw(3,1)=simplify(w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 1, 0, 0, 0, 0, 0));

Jw(1,2)=simplify(w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0));
Jw(2,2)=simplify(w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0));
Jw(3,2)=simplify(w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 1, 0, 0, 0, 0));

Jw(1,3)=simplify(w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0));
Jw(2,3)=simplify(w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0));
Jw(3,3)=simplify(w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 1, 0, 0, 0));

Jw(1,4)=simplify(w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0));
Jw(2,4)=simplify(w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0));
Jw(3,4)=simplify(w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 1, 0, 0));

Jw(1,5)=simplify(w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0));
Jw(2,5)=simplify(w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0));
Jw(3,5)=simplify(w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 1, 0));

Jw(1,6)=simplify(w1(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1));
Jw(2,6)=simplify(w2(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1));
Jw(3,6)=simplify(w3(Theta1,Theta2 ,Theta3, Theta4, Theta5, Theta6, 0, 0, 0, 0, 0, 1));


JVELO=[Jv;Jw]

%% compair
Theta1 = 0.1;
Theta2 = 0.1;
Theta3 = 0.1;
Theta4 = 0.1;
Theta5 = 0.1;
Theta6 = 0.1;

 eval(JALTER) - eval(JVELO)

% simplify(det(JVELO));

%% %% singularity

% Jacobian ALTERNATIV

syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 =simplify( [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1]) ;
T12 =simplify( [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1]) ;
T23 =simplify( [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]) ;
T34 =simplify( [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1]) ;
T45 =simplify( [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1]) ;
T56 =simplify( [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1]) ;

T02 =simplify( T01 * T12) ;

T03 =simplify( T01 * T12 * T23) ;

T04 =simplify( T01 * T12 * T23 * T34) ;

T05=simplify(  T01 * T12 * T23 * T34 * T45) ;

T06=simplify(  T01 * T12 * T23 * T34 * T45 * T56) ;

Z1=T01(1:3,3);

Z2=T02(1:3,3);

Z3=T03(1:3,3);

Z4=T04(1:3,3);

Z5=T05(1:3,3);

Z6=T06(1:3,3);

O1=T01(1:3,4);

O2=T02(1:3,4);

O3=T03(1:3,4);

O4=T04(1:3,4);

O5=T05(1:3,4);

O6=T06(1:3,4);

J1=simplify([cross(Z1,(O6-O1));Z1]);

J2=simplify([cross(Z2,(O6-O2));Z2]);

J3=simplify([cross(Z3,(O6-O3));Z3]);

J4=simplify([cross(Z4,(O6-O4));Z4]);

J5=simplify([cross(Z5,(O6-O5));Z5]);

J6=simplify([cross(Z6,(O6-O6));Z6]);

JALTER=[simplify(J1) simplify(J2) simplify(J3) simplify(J4) simplify(J5) simplify(J6)]

det_JALTER = simplify(det(JALTER))

% Theta1=0.1;
% Theta2=0.1;
% Theta3=0.1;
% Theta4=0.1;
% Theta5=0;
% Theta6=0.1;
%  eval(det_JALTER)

% det_JJ_nT5 = 15*sin(Theta3) - 15*sin(Theta2) + 15*cos(Theta3)^2*sin(Theta2) + 79*cos(Theta2)*sin(Theta3) + 15*cos(Theta2)*cos(Theta3)*sin(Theta3)
%% deteval
% Theta5 = .2
% Theta3 = .3
% Theta2 = .4
% 
% eval(det_JJ1)

%% singular t3 t2
simplify(15*sin(Theta3) - 15*sin(Theta2) + 15*cos(Theta3)^2*sin(Theta2) + 79*cos(Theta2)*sin(Theta3) + 15*cos(Theta2)*cos(Theta3)*sin(The

%% oyler-newton

clc
clear


syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 =simplify( [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1]) ;
T12 =simplify( [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1]); 
T23 =simplify( [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]) ;
T34 =simplify( [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1]) ;
T45 =simplify( [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1]) ;
T56 =simplify( [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1]); 

T02 =simplify( T01 * T12) ;

T03 =simplify( T01 * T12 * T23) ;

T04 =simplify( T01 * T12 * T23 * T34) ;

T05=simplify(  T01 * T12 * T23 * T34 * T45) ;

T06=simplify(  T01 * T12 * T23 * T34 * T45 * T56) ;

%velocity propagation


syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 =simplify( [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1]) ;
T12 =simplify( [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1]);
T23 =simplify( [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]);
T34 =simplify( [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1]);
T45 =simplify( [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1]);
T56 =simplify( [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1]);

P01=T01(1:3,4);

P12=T12(1:3,4);

P23=T23(1:3,4);

P34=T34(1:3,4);

P45=T45(1:3,4);

P56=T56(1:3,4);

R01=T01(1:3,1:3);

R12=T12(1:3,1:3);

R23=T23(1:3,1:3);

R34=T34(1:3,1:3);

R45=T45(1:3,1:3);

R56=T56(1:3,1:3);

R10=transpose(R01);

R21=transpose(R12);

R32=transpose(R23);

R43=transpose(R34);

R54=transpose(R45);

R65=transpose(R56);


%rotational velocity

syms Thetad1 Thetad2 Thetad3 Thetad4 Thetad5 Thetad6 Thetadd1 Thetadd2 Thetadd3 Thetadd4 Thetadd5 Thetadd6

w00=[0;0;0];

w11=R10*w00+[0;0;Thetad1];

w22=R21*w11+[0;0;Thetad2];

w33=R32*w22+[0;0;Thetad3];

w44=R43*w33+[0;0;Thetad4];

w55=R54*w44+[0;0;Thetad5];

w66=simplify(R65*w55+[0;0;Thetad6]);

%linear velocity

v00=[0;0;0];

v11=R10*(v00+cross(w00,P01));

v22=R21*(v11+cross(w11,P12));

v33=R32*(v22+cross(w22,P23));

v44=R43*(v33+cross(w33,P34));

v55=R54*(v44+cross(w44,P45));

v66=R65*(v55+cross(w55,P56));

R06=R01 * R12 * R23 * R34 * R45 * R56 ;

v06=simplify( R06 * v66) ; 
w06=simplify(R06 * w66) ;
% w dot

w11_d = simplify(R10*[0;0;0]+R10*(cross([0;0;0],[0;0;Thetad1]))+[0;0;Thetadd1]);
w22_d = simplify(R21*w11+R21*(cross(w11,[0;0;Thetad2]))+[0;0;Thetadd2]);
w33_d = simplify(R32*w22+R32*(cross(w22,[0;0;Thetad3]))+[0;0;Thetadd3]);
w44_d = simplify(R43*w33+R43*(cross(w33,[0;0;Thetad4]))+[0;0;Thetadd4]);
w55_d = simplify(R54*w44+R54*(cross(w44,[0;0;Thetad5]))+[0;0;Thetadd5]);
w66_d = simplify(R65*w55+R65*(cross(w55,[0;0;Thetad6]))+[0;0;Thetadd6]);

% v dot

v11_d=R10*(cross([0;0;0],P01)+cross([0;0;0],cross([0;0;0],P01))+[0;0;9.81]);
v22_d=R21*(cross(w11_d,P12)+cross(w11,cross(w11,P12))+v11_d);
v33_d=R32*(cross(w22_d,P23)+cross(w22,cross(w22,P23))+v22_d);
v44_d=R43*(cross(w33_d,P34)+cross(w33,cross(w33,P34))+v33_d);
v55_d=R54*(cross(w44_d,P45)+cross(w44,cross(w44,P45))+v44_d);
v66_d=R65*(cross(w55_d,P56)+cross(w55,cross(w55,P56))+v55_d);

v1c1_d=cross(w11_d,[0.075;0;0])+cross(w11,cross(w11,[0.075;0;0]))+v11_d;
v2c2_d=cross(w11_d,[0.395;0;0])+cross(w11,cross(w11,[0.395;0;0]))+v22_d;
v3c3_d=cross(w22_d,[0.075;-0.125;0])+cross(w22,cross(w22,[0.075;-0.125;0]))+v33_d;
v4c4_d=cross(w33_d,[0;0;0.557])+cross(w33,cross(w33,[0;0;0.557]))+v44_d;
v5c5_d=cross(w44_d,[0;-0.05;0])+cross(w44,cross(w44,[0;-0.05;0]))+v55_d;
v6c6_d=cross(w55_d,[0;0;0])+cross(w55,cross(w55,[0;0;0]))+v66_d;

syms m6

m1 = 77.42;
m2 = 125.823;
m3 = 168.077;
m4 = 75.908;
m5 = 6.939;

F11 = v1c1_d * m1;
F22 = v2c2_d * m2;
F33 = v3c3_d * m3;
F44 = v4c4_d * m4;
F55 = v5c5_d * m5;
F66 = v6c6_d * m6;

syms I66

I11 = [(m1/12)*((20/1000)^2+(20/1000)^2),0,0;0,(m1/12)*((150/1000)^2+(20/1000)^2),0;0,0,(m1/12)*((150/1000)^2+(20/1000)^2)];
I22 = [(m2/12)*((20/1000)^2+(20/1000)^2),0,0;0,(m2/12)*((790/1000)^2+(20/1000)^2),0;0,0,(m2/12)*((790/1000)^2+(20/1000)^2)];
I33 = [(m3/12)*((20/1000)^2+(20/1000)^2),0,0;0,(m3/12)*((297.36/1000)^2+(20/1000)^2),0;0,0,(m3/12)*((297.36/1000)^2+(20/1000)^2)];
I44 = [(m4/12)*((603.25/1000)^2+(20/1000)^2),0,0;0,(m4/12)*((603.25/1000)^2+(20/1000)^2),0;0,0,(m4/12)*((20/1000)^2+(20/1000)^2)];
I55 = [(m5/12)*((100/1000)^2+(20/1000)^2),0,0;0,(m5/12)*((20/1000)^2+(20/1000)^2),0;0,0,(m5/12)*((100/1000)^2+(20/1000)^2)];

N11 = I11 * w11_d + cross(w11,I11*w11);
N22 = I22 * w22_d + cross(w22,I22*w22);
N33 = I33 * w33_d + cross(w33,I33*w33);
N44 = I44 * w44_d + cross(w44,I44*w44);
N55 = I55 * w55_d + cross(w55,I55*w55);
N66 = I66 * w66_d + cross(w66,I66*w66);

f66=0+F66;
f55=R56*f66+F55;
f44=R45*f55+F44;
f33=R34*f44+F33;
f22=R23*f33+F22;
f11=R12*f22+F11;

n66=N66;
n55=N55+R56*n66+cross([0;-0.05;0],F55)+cross([0;-0.1;0],R56*f66);
n44=N44+R45*n55+cross([0;0;.558],F44)+cross([0;0;0.859],R45*f55);
n33=N33+R34*n44+cross([0.075;-0.128;0],F33)+cross([0.15;-0.256;0],R34*f44);
n22=N22+R23*n33+cross([0.395;0;0],F22)+cross([0.79;0;0],R23*f33);
n11=N11+R12*n22+cross([0.075;0;0],F11)+cross([0.15;0;0],R12*f22);

taw6=transpose(n66)*[0;0;1];
taw5=transpose(n55)*[0;0;1];
taw4=transpose(n44)*[0;0;1];
taw3=transpose(n33)*[0;0;1];
taw2=transpose(n22)*[0;0;1];
taw1=transpose(n11)*[0;0;1];

%% taw matrix

taw =[vpa(simplify(taw1),3) ;vpa(simplify(taw2),3) ;vpa(simplify(taw3),3) ;vpa(simplify(taw4),3) ;vpa(simplify(taw5),3) ;vpa(simplify(taw6),3)]

%% G

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

g(1)=eval(taw1);
g(2)=eval(taw2);
g(3)=eval(taw3);
g(4)=eval(taw4);
g(5)=eval(taw5);
g(6)=eval(taw6);

G=transpose(g)


%% M

%satr1

Thetadd1=1;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

M(1,1)=eval(taw1)-G(1)

Thetadd1=0;
Thetadd2=1;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(1,2)=eval(taw1)-G(1)

Thetadd1=0;
Thetadd2=0;
Thetadd3=1;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(1,3)=eval(taw1)-G(1)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=1;
Thetadd5=0;
Thetadd6=0;

M(1,4)=eval(taw1)-G(1)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=1;
Thetadd6=0;

M(1,5)=eval(taw1)-G(1)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=1;

M(1,6)=eval(taw1)-G(1)

%satr2

Thetadd1=1;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

M(2,1)=eval(taw2)-G(2)

Thetadd1=0;
Thetadd2=1;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(2,2)=eval(taw2)-G(2)

Thetadd1=0;
Thetadd2=0;
Thetadd3=1;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(2,3)=eval(taw2)-G(2)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=1;
Thetadd5=0;
Thetadd6=0;

M(2,4)=eval(taw2)-G(2)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=1;
Thetadd6=0;

M(2,5)=eval(taw2)-G(2)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=1;

M(2,6)=eval(taw2)-G(2)

%satr3


Thetadd1=1;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

M(3,1)=eval(taw3)-G(3)

Thetadd1=0;
Thetadd2=1;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(3,2)=eval(taw3)-G(3)

Thetadd1=0;
Thetadd2=0;
Thetadd3=1;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(3,3)=eval(taw3)-G(3)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=1;
Thetadd5=0;
Thetadd6=0;

M(3,4)=eval(taw3)-G(3)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=1;
Thetadd6=0;

M(3,5)=eval(taw3)-G(3)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=1;

M(3,6)=eval(taw3)-G(3)

%satr4

Thetadd1=1;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

M(4,1)=eval(taw4)-G(4)

Thetadd1=0;
Thetadd2=1;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(4,2)=eval(taw4)-G(4)

Thetadd1=0;
Thetadd2=0;
Thetadd3=1;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(4,3)=eval(taw4)-G(4)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=1;
Thetadd5=0;
Thetadd6=0;

M(4,4)=eval(taw4)-G(4)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=1;
Thetadd6=0;

M(4,5)=eval(taw4)-G(4)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=1;

M(4,6)=eval(taw4)-G(4)

%satr5


Thetadd1=1;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

M(5,1)=eval(taw5)-G(5)

Thetadd1=0;
Thetadd2=1;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(5,2)=eval(taw5)-G(5)

Thetadd1=0;
Thetadd2=0;
Thetadd3=1;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(5,3)=eval(taw5)-G(5)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=1;
Thetadd5=0;
Thetadd6=0;

M(5,4)=eval(taw5)-G(5)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=1;
Thetadd6=0;

M(5,5)=eval(taw5)-G(5)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=1;

M(5,6)=eval(taw5)-G(5)

%satr6


Thetadd1=1;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

M(6,1)=eval(taw6)-G(6)

Thetadd1=0;
Thetadd2=1;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(6,2)=eval(taw6)-G(6)

Thetadd1=0;
Thetadd2=0;
Thetadd3=1;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

M(6,3)=eval(taw6)-G(6)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=1;
Thetadd5=0;
Thetadd6=0;

M(6,4)=eval(taw6)-G(6)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=1;
Thetadd6=0;

M(6,5)=eval(taw6)-G(6)

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=1;

M(6,6)=eval(taw6)-G(6)


%% C 
%satr1

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=1;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(1,1)=eval(taw1)-G(1)

Thetad1=0;
Thetad2=1;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(1,2)=eval(taw1)-G(1)

Thetad1=0;
Thetad2=0;
Thetad3=1;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(1,3)=eval(taw1)-G(1)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=1;
Thetad5=0;
Thetad6=0;

C(1,4)=eval(taw1)-G(1)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=1;
Thetad6=0;

C(1,5)=eval(taw1)-G(1)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=1;

C(1,6)=eval(taw1)-G(1)

%satr2

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=1;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(2,1)=eval(taw2)-G(2)

Thetad1=0;
Thetad2=1;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(2,2)=eval(taw2)-G(2)

Thetad1=0;
Thetad2=0;
Thetad3=1;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(2,3)=eval(taw2)-G(2)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=1;
Thetad5=0;
Thetad6=0;

C(2,4)=eval(taw2)-G(2)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=1;
Thetad6=0;

C(2,5)=eval(taw2)-G(2)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=1;

C(2,6)=eval(taw2)-G(2)

%satr3

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=1;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(3,1)=eval(taw3)-G(3)

Thetad1=0;
Thetad2=1;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(3,2)=eval(taw3)-G(3)

Thetad1=0;
Thetad2=0;
Thetad3=1;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(3,3)=eval(taw3)-G(3)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=1;
Thetad5=0;
Thetad6=0;

C(3,4)=eval(taw3)-G(3)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=1;
Thetad6=0;

C(3,5)=eval(taw3)-G(3)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=1;

C(3,6)=eval(taw3)-G(3)

%satr4

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=1;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(4,1)=eval(taw4)-G(4)

Thetad1=0;
Thetad2=1;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(4,2)=eval(taw4)-G(4)

Thetad1=0;
Thetad2=0;
Thetad3=1;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(4,3)=eval(taw4)-G(4)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=1;
Thetad5=0;
Thetad6=0;

C(4,4)=eval(taw4)-G(4)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=1;
Thetad6=0;

C(4,5)=eval(taw4)-G(4)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=1;

C(4,6)=eval(taw4)-G(4)

%satr5

Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=1;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(5,1)=eval(taw5)-G(5)

Thetad1=0;
Thetad2=1;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(5,2)=eval(taw5)-G(5)

Thetad1=0;
Thetad2=0;
Thetad3=1;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(5,3)=eval(taw5)-G(5)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=1;
Thetad5=0;
Thetad6=0;

C(5,4)=eval(taw1)-G(5)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=1;
Thetad6=0;

C(5,5)=eval(taw5)-G(5)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=1;

C(5,6)=eval(taw5)-G(5)

%satr6
Thetadd1=0;
Thetadd2=0;
Thetadd3=0;
Thetadd4=0;
Thetadd5=0;
Thetadd6=0;

Thetad1=1;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(6,1)=eval(taw6)-G(6)

Thetad1=0;
Thetad2=1;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(6,2)=eval(taw6)-G(6)

Thetad1=0;
Thetad2=0;
Thetad3=1;
Thetad4=0;
Thetad5=0;
Thetad6=0;

C(6,3)=eval(taw6)-G(6)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=1;
Thetad5=0;
Thetad6=0;

C(6,4)=eval(taw6)-G(6)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=1;
Thetad6=0;

C(6,5)=eval(taw6)-G(6)

Thetad1=0;
Thetad2=0;
Thetad3=0;
Thetad4=0;
Thetad5=0;
Thetad6=1;

C(6,6)=eval(taw6)-G(6)

clear Thetad1 Theatad2 Thetad3 Thetad4 Thetad5 Thetad6

syms Thetad1 Theatad2 Thetad3 Thetad4 Thetad5 Thetad6

C=C*[Thetad1;Thetad2;Thetad3;Thetad4;Thetad5;Thetad6]

%% lagrangian

clc
clear


syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 =simplify( [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1]) 
T12 =simplify( [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1])
T23 =simplify( [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1])
T34 =simplify( [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1])
T45 =simplify( [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1])
T56 =simplify( [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1])
Tc55=[1 0 0 0;0 1 0 -0.05;0 0 1 0;0 0 0 1];

Tc44=[1 0 0 0;0 1 0 0;0 0 1 0.55835;0 0 0 1];

Tc33=[1 0 0 0.075;0 1 0 -0.128;0 0 1 0;0 0 0 1];

Tc22=[1 0 0 0.395;0 1 0 0;0 0 1 0;0 0 0 1];

Tc11=[1 0 0 0.075;0 1 0 0;0 0 1 0;0 0 0 1];

T0C1=simplify(T01*Tc11);

T02 =simplify( T01 * T12); 

T0C2=simplify(T02*Tc22);

T03 =simplify( T01 * T12 * T23); 

T0C3=simplify(T03*Tc33);

T04 =simplify( T01 * T12 * T23 * T34);

T0C4=simplify(T04*Tc44);

T05=simplify(  T01 * T12 * T23 * T34 * T45); 

T0C5=simplify(T05*Tc55);

T06=simplify(  T01 * T12 * T23 * T34 * T45 * T56); 

T0C6=T06;

Z1=T01(1:3,3);

Z2=T02(1:3,3);

Z3=T03(1:3,3);

Z4=T04(1:3,3);

Z5=T05(1:3,3);

Z6=T06(1:3,3);

O1=T01(1:3,4);

Oc1=T0C1(1:3,4);

O2=T02(1:3,4);

Oc2=T0C2(1:3,4);

O3=T03(1:3,4);

Oc3=T0C3(1:3,4);

O4=T04(1:3,4);

Oc4=T0C4(1:3,4);

O5=T05(1:3,4);

Oc5=T0C5(1:3,4);

O6=T06(1:3,4);

Oc6=O6
% coulmn of jacobians of ee.

J1=[cross(Z1,(O6-O1));Z1];

J2=[cross(Z2,(O6-O2));Z2];

J3=[cross(Z3,(O6-O3));Z3];

J4=[cross(Z4,(O6-O4));Z4];

J5=[cross(Z5,(O6-O5));Z5];

J6=[cross(Z6,(O6-O6));Z6];

J=[J1 J2 J3 J4 J5 J6];

% Jc

zero=[0;0;0;0;0;0];

Jc6=J;

Jv6=Jc6(1:3,:);

Jw6=Jc6(4:6,:)

Jc51=[cross(Z1,(Oc5-O1));Z1];

Jc52=[cross(Z2,(Oc5-O2));Z2];

Jc53=[cross(Z3,(Oc5-O3));Z3];

Jc54=[cross(Z4,(Oc5-O4));Z4];

Jc55=[cross(Z5,(Oc5-O5));Z5];

Jc5=[Jc51 Jc52 Jc53 Jc54 Jc55 zero];

Jv5=Jc5(1:3,:);

Jw5=Jc5(4:6,:);

Jc41=[cross(Z1,(Oc4-O1));Z1];

Jc42=[cross(Z2,(Oc4-O2));Z2];

Jc43=[cross(Z3,(Oc4-O3));Z3];

Jc44=[cross(Z4,(Oc4-O4));Z4];

Jc4=[Jc41 Jc42 Jc43 Jc44 zero zero];

Jv4=Jc4(1:3,:);

Jw4=Jc4(4:6,:);

Jc31=[cross(Z1,(Oc3-O1));Z1];

Jc32=[cross(Z2,(Oc3-O2));Z2];

Jc33=[cross(Z3,(Oc3-O3));Z3];

Jc3=[Jc31 Jc32 Jc33 zero zero zero];

Jv3=Jc3(1:3,:);

Jw3=Jc3(4:6,:);

Jc21=[cross(Z1,(Oc2-O1));Z1];

Jc22=[cross(Z2,(Oc2-O2));Z2];

Jc2=[Jc21 Jc22 zero zero zero zero];

Jv2=Jc2(1:3,:);

Jw2=Jc2(4:6,:);

Jc11=[cross(Z1,(Oc1-O1));Z1];

Jc1=[Jc11 zero zero zero zero zero];

Jv1=Jc1(1:3,:);

Jw1=Jc1(4:6,:);

% M
syms m1 m2 m3 m4 m5 m6 I1 I2 I3 I4 I5 I6

R01=T01(1:3,1:3);

R02=T02(1:3,1:3);

R03=T03(1:3,1:3);

R04=T04(1:3,1:3);

R05=T05(1:3,1:3);

R06=T06(1:3,1:3);

M_theta1= m1*transpose(Jv1)*Jv1+transpose(Jw1)*R01*I1*transpose(R01)*Jw1

M_theta2= m2*transpose(Jv2)*Jv2+transpose(Jw2)*R02*I2*transpose(R02)*Jw2

M_theta3= m3*transpose(Jv3)*Jv3+transpose(Jw3)*R03*I3*transpose(R03)*Jw3

M_theta4= m4*transpose(Jv4)*Jv4+transpose(Jw4)*R04*I4*transpose(R04)*Jw4

M_theta5= m5*transpose(Jv5)*Jv5+transpose(Jw5)*R05*I5*transpose(R05)*Jw5

M_theta6= m6*transpose(Jv6)*Jv6+transpose(Jw6)*R06*I6*transpose(R06)*Jw6

M_theta=M_theta1+M_theta2+M_theta3+M_theta4+M_theta5+M_theta6

% C
 
% k=1 n1

C_111=0.5*(((diff(M_theta(1,1),Theta1)))+(diff(M_theta(1,1),Theta1))-(diff(M_theta(1,1),Theta1)))

C_121=0.5*(((diff(M_theta(1,2),Theta1)))+(diff(M_theta(1,1),Theta2))-(diff(M_theta(1,2),Theta1)))

C_211=C_121

C_131=0.5*(((diff(M_theta(1,3),Theta1)))+(diff(M_theta(1,1),Theta3))-(diff(M_theta(1,3),Theta1)))

C_311=C_131

C_141=0.5*(((diff(M_theta(1,4),Theta1)))+(diff(M_theta(1,1),Theta4))-(diff(M_theta(1,4),Theta1)))

C_411=C_141

C_151=0.5*(((diff(M_theta(1,5),Theta1)))+(diff(M_theta(1,1),Theta5))-(diff(M_theta(1,5),Theta1)))

C_511=C_151

C_161=0.5*(((diff(M_theta(1,6),Theta1)))+(diff(M_theta(1,1),Theta6))-(diff(M_theta(1,6),Theta1)))

C_611=C_161

% k=1 n2

C_221=0.5*(((diff(M_theta(1,2),Theta2)))+(diff(M_theta(1,2),Theta2))-(diff(M_theta(2,2),Theta1)))

C_231=0.5*(((diff(M_theta(1,3),Theta2)))+(diff(M_theta(1,2),Theta3))-(diff(M_theta(2,3),Theta1)))

C_321=C_231

C_241=0.5*(((diff(M_theta(1,4),Theta2)))+(diff(M_theta(1,2),Theta4))-(diff(M_theta(2,4),Theta1)))

C_421=C_241

C_251=0.5*(((diff(M_theta(1,5),Theta2)))+(diff(M_theta(1,2),Theta5))-(diff(M_theta(2,5),Theta1)))

C_521=C_251

C_261=0.5*(((diff(M_theta(1,6),Theta2)))+(diff(M_theta(1,2),Theta6))-(diff(M_theta(2,6),Theta1)))

C_621=C_261

% k=1 n3

C_331=0.5*(((diff(M_theta(1,3),Theta3)))+(diff(M_theta(1,3),Theta3))-(diff(M_theta(3,3),Theta1)))

C_341=0.5*(((diff(M_theta(1,4),Theta3)))+(diff(M_theta(1,3),Theta4))-(diff(M_theta(3,4),Theta1)))

C_431=C_341

C_351=0.5*(((diff(M_theta(1,5),Theta3)))+(diff(M_theta(1,3),Theta5))-(diff(M_theta(3,5),Theta1)))

C_531=C_351

C_361=0.5*(((diff(M_theta(1,6),Theta3)))+(diff(M_theta(1,3),Theta6))-(diff(M_theta(3,6),Theta1)))

C_631=C_361

% k=1 n4

C_441=0.5*(((diff(M_theta(1,4),Theta4)))+(diff(M_theta(1,4),Theta4))-(diff(M_theta(4,4),Theta1)))

C_451=0.5*(((diff(M_theta(1,5),Theta4)))+(diff(M_theta(1,4),Theta5))-(diff(M_theta(4,5),Theta1)))

C_541=C_451

C_461=0.5*(((diff(M_theta(1,6),Theta4)))+(diff(M_theta(1,4),Theta6))-(diff(M_theta(4,6),Theta1)))

C_641=C_461
% k=1 n5

C_551=0.5*(((diff(M_theta(1,5),Theta5)))+(diff(M_theta(1,5),Theta5))-(diff(M_theta(5,5),Theta1)))

C_561=0.5*(((diff(M_theta(1,6),Theta5)))+(diff(M_theta(1,5),Theta6))-(diff(M_theta(5,6),Theta1)))

C_651=C_561

% k=1 n6

C_661=0.5*(((diff(M_theta(1,6),Theta6)))+(diff(M_theta(1,6),Theta6))-(diff(M_theta(6,6),Theta1)))

% k=2 n1

C_112=0.5*(((diff(M_theta(2,1),Theta1)))+(diff(M_theta(2,1),Theta1))-(diff(M_theta(1,1),Theta2)))

C_122=0.5*(((diff(M_theta(2,2),Theta1)))+(diff(M_theta(2,1),Theta2))-(diff(M_theta(1,2),Theta2)))

C_212=C_122

C_132=0.5*(((diff(M_theta(2,3),Theta1)))+(diff(M_theta(2,1),Theta3))-(diff(M_theta(1,3),Theta2)))

C_312=C_132

C_142=0.5*(((diff(M_theta(2,4),Theta1)))+(diff(M_theta(2,1),Theta4))-(diff(M_theta(1,4),Theta2)))

C_412=C_142

C_152=0.5*(((diff(M_theta(2,5),Theta1)))+(diff(M_theta(2,1),Theta5))-(diff(M_theta(1,5),Theta2)))

C_512=C_152

C_162=0.5*(((diff(M_theta(2,6),Theta1)))+(diff(M_theta(2,1),Theta6))-(diff(M_theta(1,6),Theta2)))

C_612=C_162
% k=2 n2

C_222=0.5*(((diff(M_theta(2,2),Theta2)))+(diff(M_theta(2,2),Theta2))-(diff(M_theta(2,2),Theta2)))

C_232=0.5*(((diff(M_theta(2,3),Theta2)))+(diff(M_theta(2,2),Theta3))-(diff(M_theta(2,3),Theta2)))

C_322=C_232

C_242=0.5*(((diff(M_theta(2,4),Theta2)))+(diff(M_theta(2,2),Theta4))-(diff(M_theta(2,4),Theta2)))

C_422=C_242

C_252=0.5*(((diff(M_theta(2,5),Theta2)))+(diff(M_theta(2,2),Theta5))-(diff(M_theta(2,5),Theta2)))

C_522=C_252

C_262=0.5*(((diff(M_theta(2,6),Theta2)))+(diff(M_theta(2,2),Theta6))-(diff(M_theta(2,6),Theta2)))

C_622=C_262
% k=2 n3

C_332=0.5*(((diff(M_theta(2,3),Theta3)))+(diff(M_theta(2,3),Theta3))-(diff(M_theta(3,3),Theta2)))

C_342=0.5*(((diff(M_theta(2,4),Theta3)))+(diff(M_theta(2,3),Theta4))-(diff(M_theta(3,4),Theta2)))

C_432=C_342

C_352=0.5*(((diff(M_theta(2,5),Theta3)))+(diff(M_theta(2,3),Theta5))-(diff(M_theta(3,5),Theta2)))

C_532=C_352

C_362=0.5*(((diff(M_theta(2,6),Theta3)))+(diff(M_theta(2,3),Theta6))-(diff(M_theta(3,6),Theta2)))

C_632=C_362

% k=2 n4

C_442=0.5*(((diff(M_theta(2,4),Theta4)))+(diff(M_theta(2,4),Theta4))-(diff(M_theta(4,4),Theta2)))

C_452=0.5*(((diff(M_theta(2,5),Theta4)))+(diff(M_theta(2,4),Theta5))-(diff(M_theta(4,5),Theta2)))

C_542=C_452

C_462=0.5*(((diff(M_theta(2,6),Theta4)))+(diff(M_theta(2,4),Theta6))-(diff(M_theta(4,6),Theta2)))

C_642=C_462

% k=2 n5

C_552=0.5*(((diff(M_theta(2,5),Theta5)))+(diff(M_theta(2,5),Theta5))-(diff(M_theta(5,5),Theta2)))

C_562=0.5*(((diff(M_theta(2,6),Theta5)))+(diff(M_theta(2,5),Theta6))-(diff(M_theta(5,6),Theta2)))

C_652=C_562
% k=2 n6

C_662=0.5*(((diff(M_theta(2,6),Theta6)))+(diff(M_theta(2,6),Theta6))-(diff(M_theta(6,6),Theta2)))

% k=3 n1

C_113=0.5*(((diff(M_theta(3,1),Theta1)))+(diff(M_theta(3,1),Theta1))-(diff(M_theta(1,1),Theta3)))

C_123=0.5*(((diff(M_theta(3,2),Theta1)))+(diff(M_theta(3,1),Theta2))-(diff(M_theta(1,2),Theta3)))

C_213=C_123

C_133=0.5*(((diff(M_theta(3,3),Theta1)))+(diff(M_theta(3,1),Theta3))-(diff(M_theta(1,3),Theta3)))

C_313=C_133

C_143=0.5*(((diff(M_theta(3,4),Theta1)))+(diff(M_theta(3,1),Theta4))-(diff(M_theta(1,4),Theta3)))

C_413=C_143

C_153=0.5*(((diff(M_theta(3,5),Theta1)))+(diff(M_theta(3,1),Theta5))-(diff(M_theta(1,5),Theta3)))

C_513=C_153

C_163=0.5*(((diff(M_theta(3,6),Theta1)))+(diff(M_theta(3,1),Theta6))-(diff(M_theta(1,6),Theta3)))

C_613=C_163
% k=3 n2

C_223=0.5*(((diff(M_theta(3,2),Theta2)))+(diff(M_theta(3,2),Theta2))-(diff(M_theta(2,2),Theta3)))

C_233=0.5*(((diff(M_theta(3,3),Theta2)))+(diff(M_theta(3,2),Theta3))-(diff(M_theta(2,3),Theta3)))

C_323=C_233

C_243=0.5*(((diff(M_theta(3,4),Theta2)))+(diff(M_theta(3,2),Theta4))-(diff(M_theta(2,4),Theta3)))

C_423=C_243

C_253=0.5*(((diff(M_theta(3,5),Theta2)))+(diff(M_theta(3,2),Theta5))-(diff(M_theta(2,5),Theta3)))

C_523=C_253

C_263=0.5*(((diff(M_theta(3,6),Theta2)))+(diff(M_theta(3,2),Theta6))-(diff(M_theta(2,6),Theta3)))

C_623=C_263
% k=3 n3

C_333=0.5*(((diff(M_theta(3,3),Theta3)))+(diff(M_theta(3,3),Theta3))-(diff(M_theta(3,3),Theta3)))

C_343=0.5*(((diff(M_theta(3,4),Theta3)))+(diff(M_theta(3,3),Theta4))-(diff(M_theta(3,4),Theta3)))

C_433=C_343

C_353=0.5*(((diff(M_theta(3,5),Theta3)))+(diff(M_theta(3,3),Theta5))-(diff(M_theta(3,5),Theta3)))

C_533=C_353

C_363=0.5*(((diff(M_theta(3,6),Theta3)))+(diff(M_theta(3,3),Theta6))-(diff(M_theta(3,6),Theta3)))

C_633=C_363

% k=3 n4

C_443=0.5*(((diff(M_theta(3,4),Theta4)))+(diff(M_theta(3,4),Theta4))-(diff(M_theta(4,4),Theta3)))

C_453=0.5*(((diff(M_theta(3,5),Theta4)))+(diff(M_theta(3,4),Theta5))-(diff(M_theta(4,5),Theta3)))

C_543=C_453

C_463=0.5*(((diff(M_theta(3,6),Theta4)))+(diff(M_theta(3,4),Theta6))-(diff(M_theta(4,6),Theta3)))

C_643=C_463

% k=3 n5

C_553=0.5*(((diff(M_theta(3,5),Theta5)))+(diff(M_theta(3,5),Theta5))-(diff(M_theta(5,5),Theta3)))

C_563=0.5*(((diff(M_theta(3,6),Theta5)))+(diff(M_theta(3,5),Theta6))-(diff(M_theta(5,6),Theta3)))

C_653=C_563
% k=3 n6

C_663=0.5*(((diff(M_theta(3,6),Theta6)))+(diff(M_theta(3,6),Theta6))-(diff(M_theta(6,6),Theta3)))

% k=4 n1

C_114=0.5*(((diff(M_theta(4,1),Theta1)))+(diff(M_theta(4,1),Theta1))-(diff(M_theta(1,1),Theta4)))

C_124=0.5*(((diff(M_theta(4,2),Theta1)))+(diff(M_theta(4,1),Theta2))-(diff(M_theta(1,2),Theta4)))

C_214=C_124

C_134=0.5*(((diff(M_theta(4,3),Theta1)))+(diff(M_theta(4,1),Theta3))-(diff(M_theta(1,3),Theta4)))

C_314=C_134

C_144=0.5*(((diff(M_theta(4,4),Theta1)))+(diff(M_theta(4,1),Theta4))-(diff(M_theta(1,4),Theta4)))

C_414=C_144

C_154=0.5*(((diff(M_theta(4,5),Theta1)))+(diff(M_theta(4,1),Theta5))-(diff(M_theta(1,5),Theta4)))

C_514=C_154

C_164=0.5*(((diff(M_theta(4,6),Theta1)))+(diff(M_theta(4,1),Theta6))-(diff(M_theta(1,6),Theta4)))

C_614=C_164
% k=4 n2

C_224=0.5*(((diff(M_theta(4,2),Theta2)))+(diff(M_theta(4,2),Theta2))-(diff(M_theta(2,2),Theta4)))

C_234=0.5*(((diff(M_theta(4,3),Theta2)))+(diff(M_theta(4,2),Theta3))-(diff(M_theta(2,3),Theta4)))

C_324=C_234

C_244=0.5*(((diff(M_theta(4,4),Theta2)))+(diff(M_theta(4,2),Theta4))-(diff(M_theta(2,4),Theta4)))

C_424=C_244

C_254=0.5*(((diff(M_theta(4,5),Theta2)))+(diff(M_theta(4,2),Theta5))-(diff(M_theta(2,5),Theta4)))

C_524=C_254

C_264=0.5*(((diff(M_theta(4,6),Theta2)))+(diff(M_theta(4,2),Theta6))-(diff(M_theta(2,6),Theta4)))

C_624=C_264
% k=4 n3

C_334=0.5*(((diff(M_theta(4,3),Theta3)))+(diff(M_theta(4,3),Theta3))-(diff(M_theta(3,3),Theta4)))

C_344=0.5*(((diff(M_theta(4,4),Theta3)))+(diff(M_theta(4,3),Theta4))-(diff(M_theta(3,4),Theta4)))

C_434=C_344

C_354=0.5*(((diff(M_theta(4,5),Theta3)))+(diff(M_theta(4,3),Theta5))-(diff(M_theta(3,5),Theta4)))

C_534=C_354

C_364=0.5*(((diff(M_theta(4,6),Theta3)))+(diff(M_theta(4,3),Theta6))-(diff(M_theta(3,6),Theta4)))

C_634=C_364

% k=4 n4

C_444=0.5*(((diff(M_theta(4,4),Theta4)))+(diff(M_theta(4,4),Theta4))-(diff(M_theta(4,4),Theta4)))

C_454=0.5*(((diff(M_theta(4,5),Theta4)))+(diff(M_theta(4,4),Theta5))-(diff(M_theta(4,5),Theta4)))

C_544=C_454

C_464=0.5*(((diff(M_theta(4,6),Theta4)))+(diff(M_theta(4,4),Theta6))-(diff(M_theta(4,6),Theta4)))

C_644=C_464

% k=4 n5

C_554=0.5*(((diff(M_theta(4,5),Theta5)))+(diff(M_theta(4,5),Theta5))-(diff(M_theta(5,5),Theta4)))

C_564=0.5*(((diff(M_theta(4,6),Theta5)))+(diff(M_theta(4,5),Theta6))-(diff(M_theta(5,6),Theta4)))

C_654=C_564
% k=4 n6

C_664=0.5*(((diff(M_theta(4,6),Theta6)))+(diff(M_theta(4,6),Theta6))-(diff(M_theta(6,6),Theta4)))

% k=5 n1

C_115=0.5*(((diff(M_theta(5,1),Theta1)))+(diff(M_theta(5,1),Theta1))-(diff(M_theta(1,1),Theta5)))

C_125=0.5*(((diff(M_theta(5,2),Theta1)))+(diff(M_theta(5,1),Theta2))-(diff(M_theta(1,2),Theta5)))

C_215=C_125

C_135=0.5*(((diff(M_theta(5,3),Theta1)))+(diff(M_theta(5,1),Theta3))-(diff(M_theta(1,3),Theta5)))

C_315=C_135

C_145=0.5*(((diff(M_theta(5,4),Theta1)))+(diff(M_theta(5,1),Theta4))-(diff(M_theta(1,4),Theta5)))

C_415=C_145

C_155=0.5*(((diff(M_theta(5,5),Theta1)))+(diff(M_theta(5,1),Theta5))-(diff(M_theta(1,5),Theta5)))

C_515=C_155

C_165=0.5*(((diff(M_theta(5,6),Theta1)))+(diff(M_theta(5,1),Theta6))-(diff(M_theta(1,6),Theta5)))

C_615=C_165
% k=5 n2

C_225=0.5*(((diff(M_theta(5,2),Theta2)))+(diff(M_theta(5,2),Theta2))-(diff(M_theta(2,2),Theta5)))

C_235=0.5*(((diff(M_theta(5,3),Theta2)))+(diff(M_theta(5,2),Theta3))-(diff(M_theta(2,3),Theta5)))

C_325=C_235

C_245=0.5*(((diff(M_theta(5,4),Theta2)))+(diff(M_theta(5,2),Theta4))-(diff(M_theta(2,4),Theta5)))

C_425=C_245

C_255=0.5*(((diff(M_theta(5,5),Theta2)))+(diff(M_theta(5,2),Theta5))-(diff(M_theta(2,5),Theta5)))

C_525=C_255

C_265=0.5*(((diff(M_theta(5,6),Theta2)))+(diff(M_theta(5,2),Theta6))-(diff(M_theta(2,6),Theta5)))

C_625=C_265
% k=5 n3

C_335=0.5*(((diff(M_theta(5,3),Theta3)))+(diff(M_theta(5,3),Theta3))-(diff(M_theta(3,3),Theta5)))

C_345=0.5*(((diff(M_theta(5,4),Theta3)))+(diff(M_theta(5,3),Theta4))-(diff(M_theta(3,4),Theta5)))

C_435=C_345

C_355=0.5*(((diff(M_theta(5,5),Theta3)))+(diff(M_theta(5,3),Theta5))-(diff(M_theta(3,5),Theta5)))

C_535=C_355

C_365=0.5*(((diff(M_theta(5,6),Theta3)))+(diff(M_theta(5,3),Theta6))-(diff(M_theta(3,6),Theta5)))

C_635=C_365

% k=5 n4

C_445=0.5*(((diff(M_theta(5,4),Theta4)))+(diff(M_theta(5,4),Theta4))-(diff(M_theta(4,4),Theta5)))

C_455=0.5*(((diff(M_theta(5,5),Theta4)))+(diff(M_theta(5,4),Theta5))-(diff(M_theta(4,5),Theta5)))

C_545=C_455

C_465=0.5*(((diff(M_theta(5,6),Theta4)))+(diff(M_theta(5,4),Theta6))-(diff(M_theta(4,6),Theta5)))

C_645=C_465

% k=5 n5

C_555=0.5*(((diff(M_theta(5,5),Theta5)))+(diff(M_theta(5,5),Theta5))-(diff(M_theta(5,5),Theta5)))

C_565=0.5*(((diff(M_theta(5,6),Theta5)))+(diff(M_theta(5,5),Theta6))-(diff(M_theta(5,6),Theta5)))

C_655=C_565
% k=5 n6

C_665=0.5*(((diff(M_theta(5,6),Theta6)))+(diff(M_theta(5,6),Theta6))-(diff(M_theta(6,6),Theta5)))

% k=6 n1

C_116=0.5*(((diff(M_theta(6,1),Theta1)))+(diff(M_theta(6,1),Theta1))-(diff(M_theta(1,1),Theta6)))

C_126=0.5*(((diff(M_theta(6,2),Theta1)))+(diff(M_theta(6,1),Theta2))-(diff(M_theta(1,2),Theta6)))

C_216=C_126

C_136=0.5*(((diff(M_theta(6,3),Theta1)))+(diff(M_theta(6,1),Theta3))-(diff(M_theta(1,3),Theta6)))

C_316=C_136

C_146=0.5*(((diff(M_theta(6,4),Theta1)))+(diff(M_theta(6,1),Theta4))-(diff(M_theta(1,4),Theta6)))

C_416=C_146

C_156=0.5*(((diff(M_theta(6,5),Theta1)))+(diff(M_theta(6,1),Theta5))-(diff(M_theta(1,5),Theta6)))

C_516=C_156

C_166=0.5*(((diff(M_theta(6,6),Theta1)))+(diff(M_theta(6,1),Theta6))-(diff(M_theta(1,6),Theta6)))

C_616=C_166
% k=6 n2

C_226=0.5*(((diff(M_theta(6,2),Theta2)))+(diff(M_theta(6,2),Theta2))-(diff(M_theta(2,2),Theta6)))

C_236=0.5*(((diff(M_theta(6,3),Theta2)))+(diff(M_theta(6,2),Theta3))-(diff(M_theta(2,3),Theta6)))

C_326=C_236

C_246=0.5*(((diff(M_theta(6,4),Theta2)))+(diff(M_theta(6,2),Theta4))-(diff(M_theta(2,4),Theta6)))

C_426=C_246

C_256=0.5*(((diff(M_theta(6,5),Theta2)))+(diff(M_theta(6,2),Theta5))-(diff(M_theta(2,5),Theta6)))

C_526=C_256

C_266=0.5*(((diff(M_theta(6,6),Theta2)))+(diff(M_theta(6,2),Theta6))-(diff(M_theta(2,6),Theta6)))

C_626=C_266
% k=6 n3

C_336=0.5*(((diff(M_theta(6,3),Theta3)))+(diff(M_theta(6,3),Theta3))-(diff(M_theta(3,3),Theta6)))

C_346=0.5*(((diff(M_theta(6,4),Theta3)))+(diff(M_theta(6,3),Theta4))-(diff(M_theta(3,4),Theta6)))

C_436=C_346

C_356=0.5*(((diff(M_theta(6,5),Theta3)))+(diff(M_theta(6,3),Theta5))-(diff(M_theta(3,5),Theta6)))

C_536=C_356

C_366=0.5*(((diff(M_theta(6,6),Theta3)))+(diff(M_theta(6,3),Theta6))-(diff(M_theta(3,6),Theta6)))

C_636=C_366

% k=6 n4

C_446=0.5*(((diff(M_theta(6,4),Theta4)))+(diff(M_theta(6,4),Theta4))-(diff(M_theta(4,4),Theta6)))

C_456=0.5*(((diff(M_theta(6,5),Theta4)))+(diff(M_theta(6,4),Theta5))-(diff(M_theta(4,5),Theta6)))

C_546=C_456

C_466=0.5*(((diff(M_theta(6,6),Theta4)))+(diff(M_theta(6,4),Theta6))-(diff(M_theta(4,6),Theta6)))

C_646=C_466

% k=6 n5

C_556=0.5*(((diff(M_theta(6,5),Theta5)))+(diff(M_theta(6,5),Theta5))-(diff(M_theta(5,5),Theta6)))

C_566=0.5*(((diff(M_theta(6,6),Theta5)))+(diff(M_theta(6,5),Theta6))-(diff(M_theta(5,6),Theta6)))

C_656=C_566
% k=6 n6

C_666=0.5*(((diff(M_theta(6,6),Theta6)))+(diff(M_theta(6,6),Theta6))-(diff(M_theta(6,6),Theta6)))

% C final 

C(1)=C_111+C_121+C_131+C_141+C_151+C_161+C_211+C_221+C_231+C_241+C_251+C_261+C_311+C_321+C_331+C_341+C_351+C_361+C_411+C_421+C_431+C_441+C_451+C_461+C_511+C_521+C_531+C_541+C_551+C_561+C_611+C_621+C_631+C_641+C_651+C_661;
C(2)=C_112+C_122+C_132+C_142+C_152+C_162+C_212+C_222+C_232+C_242+C_252+C_262+C_312+C_322+C_332+C_342+C_352+C_362+C_412+C_422+C_432+C_442+C_452+C_462+C_512+C_522+C_532+C_542+C_552+C_562+C_612+C_622+C_632+C_642+C_652+C_662;
C(3)=C_113+C_123+C_133+C_143+C_153+C_163+C_213+C_223+C_233+C_243+C_253+C_263+C_313+C_323+C_333+C_343+C_353+C_363+C_413+C_423+C_433+C_443+C_453+C_463+C_513+C_523+C_533+C_543+C_553+C_563+C_613+C_623+C_633+C_643+C_653+C_663;
c(4)=C_114+C_124+C_134+C_144+C_154+C_164+C_214+C_224+C_234+C_244+C_254+C_264+C_314+C_324+C_334+C_344+C_354+C_364+C_414+C_424+C_434+C_444+C_454+C_464+C_514+C_524+C_534+C_544+C_554+C_564+C_614+C_624+C_634+C_644+C_654+C_664;
c(5)=C_115+C_125+C_135+C_145+C_155+C_165+C_215+C_225+C_235+C_245+C_255+C_265+C_315+C_325+C_335+C_345+C_355+C_365+C_415+C_425+C_435+C_445+C_455+C_465+C_515+C_525+C_535+C_545+C_555+C_565+C_615+C_625+C_635+C_645+C_655+C_665;
c(6)=C_116+C_126+C_136+C_146+C_156+C_166+C_216+C_226+C_236+C_246+C_256+C_266+C_316+C_326+C_336+C_346+C_356+C_366+C_416+C_426+C_436+C_446+C_456+C_466+C_516+C_526+C_536+C_546+C_556+C_566+C_616+C_626+C_636+C_646+C_656+C_666;


% G

% Pref = [0,0,0]

syms m6 gr

rc1=T0C1(1:3,4)

rc2=T0C2(1:3,4)

rc3=T0C3(1:3,4)

rc4=T0C4(1:3,4)

rc5=T0C5(1:3,4)

rc6=T0C6(1:3,4)

G=[0;0;-9.81]

m1 = 77.42;

m2 = 125.823;

m3 = 168.077;

m4 = 75.908;

m5 = 6.939;

transpose(G)

P1=-m1*transpose(G)*rc1

P2=-m2*transpose(G)*rc2

P3=-m3*transpose(G)*rc3

P4=-m4*transpose(G)*rc4

P5=-m5*transpose(G)*rc5

P6=-m6*transpose(G)*rc6

P = P1 + P2 + P3 + P4 + P5 + P6 
% k=1

g1=diff(P,Theta1)

% k=2

g2=diff(P,Theta2)

% k=3

g3=diff(P,Theta3)

% k=4 

g4=diff(P,Theta4)

% k=5 

g5=diff(P,Theta5)

% k=6

g6=diff(P,Theta6)

% taw 
syms Thetadd1 Thetadd2 Thetadd3 Thetadd4 Thetadd5 Thetadd6

syms Thetad1 Thetad2 Thetad3 Thetad4 Thetad5 Thetad6

taw1=(M_theta(1,1)*Thetadd1+M_theta(1,2)*Thetadd2+M_theta(1,3)*Thetadd3+M_theta(1,4)*Thetadd4+M_theta(1,5)*Thetadd5+M_theta(1,6)*Thetadd6)+(C_111*Thetad1*Thetad1+C_121*Thetad1*Thetad2+C_131*Thetad1*Thetad3+C_141*Thetad1*Thetad4+C_151*Thetad1*Thetad5+C_161*Thetad1*Thetad6+C_211*Thetad2*Thetad1+C_221*Thetad2*Thetad2+C_231*Thetad2*Thetad3+C_241*Thetad2*Thetad4+C_251*Thetad2*Thetad5+C_261*Thetad2*Thetad6+C_311*Thetad3*Thetad1+C_321*Thetad3*Thetad2+C_331*Thetad3*Thetad3+C_341*Thetad3*Thetad4+C_351*Thetad3*Thetad5+C_361*Thetad3*Thetad6+C_411*Thetad4*Thetad1+C_421*Thetad4*Thetad2+C_431*Thetad4*Thetad3+C_441*Thetad4*Thetad4+C_451*Thetad4*Thetad5+C_461*Thetad4*Thetad6+C_511*Thetad5*Thetad1+C_521*Thetad5*Thetad2+C_531*Thetad5*Thetad3+C_541*Thetad5*Thetad4+C_551*Thetad5*Thetad5+C_561*Thetad5*Thetad6+C_611*Thetad6*Thetad1+C_621*Thetad6*Thetad2+C_631*Thetad6*Thetad3+C_641*Thetad6*Thetad4+C_651*Thetad6*Thetad5+C_661*Thetad1*Thetad6)+g1

taw2=(M_theta(2,1)*Thetadd1+M_theta(2,2)*Thetadd2+M_theta(2,3)*Thetadd3+M_theta(2,4)*Thetadd4+M_theta(2,5)*Thetadd5+M_theta(2,6)*Thetadd6)+(C_112*Thetad1*Thetad1+C_122*Thetad1*Thetad2+C_132*Thetad1*Thetad3+C_142*Thetad1*Thetad4+C_152*Thetad1*Thetad5+C_162*Thetad1*Thetad6+C_212*Thetad2*Thetad1+C_222*Thetad2*Thetad2+C_232*Thetad2*Thetad3+C_242*Thetad2*Thetad4+C_252*Thetad2*Thetad5+C_262*Thetad2*Thetad6+C_312*Thetad3*Thetad1+C_322*Thetad3*Thetad2+C_332*Thetad3*Thetad3+C_342*Thetad3*Thetad4+C_352*Thetad3*Thetad5+C_362*Thetad3*Thetad6+C_412*Thetad4*Thetad1+C_422*Thetad4*Thetad2+C_432*Thetad4*Thetad3+C_442*Thetad4*Thetad4+C_452*Thetad4*Thetad5+C_462*Thetad4*Thetad6+C_512*Thetad5*Thetad1+C_522*Thetad5*Thetad2+C_532*Thetad5*Thetad3+C_542*Thetad5*Thetad4+C_552*Thetad5*Thetad5+C_562*Thetad5*Thetad6+C_612*Thetad6*Thetad1+C_622*Thetad6*Thetad2+C_632*Thetad6*Thetad3+C_642*Thetad6*Thetad4+C_652*Thetad6*Thetad5+C_662*Thetad1*Thetad6)+g2

taw3=(M_theta(3,1)*Thetadd1+M_theta(3,2)*Thetadd2+M_theta(3,3)*Thetadd3+M_theta(3,4)*Thetadd4+M_theta(3,5)*Thetadd5+M_theta(3,6)*Thetadd6)+(C_113*Thetad1*Thetad1+C_123*Thetad1*Thetad2+C_133*Thetad1*Thetad3+C_143*Thetad1*Thetad4+C_153*Thetad1*Thetad5+C_163*Thetad1*Thetad6+C_213*Thetad2*Thetad1+C_223*Thetad2*Thetad2+C_233*Thetad2*Thetad3+C_243*Thetad2*Thetad4+C_253*Thetad2*Thetad5+C_263*Thetad2*Thetad6+C_313*Thetad3*Thetad1+C_323*Thetad3*Thetad2+C_333*Thetad3*Thetad3+C_343*Thetad3*Thetad4+C_353*Thetad3*Thetad5+C_363*Thetad3*Thetad6+C_413*Thetad4*Thetad1+C_423*Thetad4*Thetad2+C_433*Thetad4*Thetad3+C_443*Thetad4*Thetad4+C_453*Thetad4*Thetad5+C_463*Thetad4*Thetad6+C_513*Thetad5*Thetad1+C_523*Thetad5*Thetad2+C_533*Thetad5*Thetad3+C_543*Thetad5*Thetad4+C_553*Thetad5*Thetad5+C_563*Thetad5*Thetad6+C_613*Thetad6*Thetad1+C_623*Thetad6*Thetad2+C_633*Thetad6*Thetad3+C_643*Thetad6*Thetad4+C_653*Thetad6*Thetad5+C_663*Thetad1*Thetad6)+g3

taw4=(M_theta(4,1)*Thetadd1+M_theta(4,2)*Thetadd2+M_theta(4,3)*Thetadd3+M_theta(4,4)*Thetadd4+M_theta(4,5)*Thetadd5+M_theta(4,6)*Thetadd6)+(C_114*Thetad1*Thetad1+C_124*Thetad1*Thetad2+C_134*Thetad1*Thetad3+C_144*Thetad1*Thetad4+C_154*Thetad1*Thetad5+C_164*Thetad1*Thetad6+C_214*Thetad2*Thetad1+C_224*Thetad2*Thetad2+C_234*Thetad2*Thetad3+C_244*Thetad2*Thetad4+C_254*Thetad2*Thetad5+C_264*Thetad2*Thetad6+C_314*Thetad3*Thetad1+C_324*Thetad3*Thetad2+C_334*Thetad3*Thetad3+C_344*Thetad3*Thetad4+C_354*Thetad3*Thetad5+C_364*Thetad3*Thetad6+C_414*Thetad4*Thetad1+C_424*Thetad4*Thetad2+C_434*Thetad4*Thetad3+C_444*Thetad4*Thetad4+C_454*Thetad4*Thetad5+C_464*Thetad4*Thetad6+C_514*Thetad5*Thetad1+C_524*Thetad5*Thetad2+C_534*Thetad5*Thetad3+C_544*Thetad5*Thetad4+C_554*Thetad5*Thetad5+C_564*Thetad5*Thetad6+C_614*Thetad6*Thetad1+C_624*Thetad6*Thetad2+C_634*Thetad6*Thetad3+C_644*Thetad6*Thetad4+C_654*Thetad6*Thetad5+C_664*Thetad1*Thetad6)+g4

taw5=(M_theta(5,1)*Thetadd1+M_theta(5,2)*Thetadd2+M_theta(5,3)*Thetadd3+M_theta(5,4)*Thetadd4+M_theta(5,5)*Thetadd5+M_theta(5,6)*Thetadd6)+(C_115*Thetad1*Thetad1+C_125*Thetad1*Thetad2+C_135*Thetad1*Thetad3+C_145*Thetad1*Thetad4+C_155*Thetad1*Thetad5+C_165*Thetad1*Thetad6+C_215*Thetad2*Thetad1+C_225*Thetad2*Thetad2+C_235*Thetad2*Thetad3+C_245*Thetad2*Thetad4+C_255*Thetad2*Thetad5+C_265*Thetad2*Thetad6+C_315*Thetad3*Thetad1+C_325*Thetad3*Thetad2+C_335*Thetad3*Thetad3+C_345*Thetad3*Thetad4+C_355*Thetad3*Thetad5+C_365*Thetad3*Thetad6+C_415*Thetad4*Thetad1+C_425*Thetad4*Thetad2+C_435*Thetad4*Thetad3+C_445*Thetad4*Thetad4+C_455*Thetad4*Thetad5+C_465*Thetad4*Thetad6+C_515*Thetad5*Thetad1+C_525*Thetad5*Thetad2+C_535*Thetad5*Thetad3+C_545*Thetad5*Thetad4+C_555*Thetad5*Thetad5+C_565*Thetad5*Thetad6+C_615*Thetad6*Thetad1+C_625*Thetad6*Thetad2+C_635*Thetad6*Thetad3+C_645*Thetad6*Thetad4+C_655*Thetad6*Thetad5+C_665*Thetad1*Thetad6)+g5

taw6=(M_theta(6,1)*Thetadd1+M_theta(6,2)*Thetadd2+M_theta(6,3)*Thetadd3+M_theta(6,4)*Thetadd4+M_theta(6,5)*Thetadd5+M_theta(6,6)*Thetadd6)+(C_116*Thetad1*Thetad1+C_126*Thetad1*Thetad2+C_136*Thetad1*Thetad3+C_146*Thetad1*Thetad4+C_156*Thetad1*Thetad5+C_166*Thetad1*Thetad6+C_216*Thetad2*Thetad1+C_226*Thetad2*Thetad2+C_236*Thetad2*Thetad3+C_246*Thetad2*Thetad4+C_256*Thetad2*Thetad5+C_266*Thetad2*Thetad6+C_316*Thetad3*Thetad1+C_326*Thetad3*Thetad2+C_336*Thetad3*Thetad3+C_346*Thetad3*Thetad4+C_356*Thetad3*Thetad5+C_366*Thetad3*Thetad6+C_416*Thetad4*Thetad1+C_426*Thetad4*Thetad2+C_436*Thetad4*Thetad3+C_446*Thetad4*Thetad4+C_456*Thetad4*Thetad5+C_466*Thetad4*Thetad6+C_516*Thetad5*Thetad1+C_526*Thetad5*Thetad2+C_536*Thetad5*Thetad3+C_546*Thetad5*Thetad4+C_556*Thetad5*Thetad5+C_566*Thetad5*Thetad6+C_616*Thetad6*Thetad1+C_626*Thetad6*Thetad2+C_636*Thetad6*Thetad3+C_646*Thetad6*Thetad4+C_656*Thetad6*Thetad5+C_666*Thetad1*Thetad6)+g6

taw=taw1 + taw2 + taw3 + taw4 +taw5 + taw6

%% staric force( Motor sugesstion )

% Jacobian ALTERNATIV

syms Theta1 Theta2 Theta3 Theta4 Theta5 Theta6

T01 =simplify( [cos(Theta1), -sin(Theta1), 0, 0; sin(Theta1),cos(Theta1), 0, 0; 0, 0, 1, 0;0, 0, 0, 1]) ;
T12 =simplify( [cos(Theta2), -sin(Theta2), 0, 0.15; 0,0, -1, 0; sin(Theta2), cos(Theta2), 0, 0; 0, 0, 0, 1]) ;
T23 =simplify( [cos(Theta3), -sin(Theta3), 0, 0.79; sin(Theta3),cos(Theta3), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]) ;
T34 =simplify( [cos(Theta4), -sin(Theta4), 0, 0.15; 0,0, -1, -0.86; sin(Theta4), cos(Theta4), 0, 0; 0, 0, 0, 1]) ;
T45 =simplify( [cos(Theta5), -sin(Theta5), 0, 0; 0,0, 1, 0; -sin(Theta5), -cos(Theta5), 0, 0; 0, 0, 0, 1]) ;
T56 =simplify( [cos(Theta6), -sin(Theta6), 0, 0; 0,0, -1, 0; sin(Theta6), cos(Theta6), 0, 0; 0, 0, 0, 1]) ;

T02 =simplify( T01 * T12) ;

T03 =simplify( T01 * T12 * T23) ;

T04 =simplify( T01 * T12 * T23 * T34) ;

T05=simplify(  T01 * T12 * T23 * T34 * T45) ;

T06=simplify(  T01 * T12 * T23 * T34 * T45 * T56) ;

Z1=T01(1:3,3);

Z2=T02(1:3,3);

Z3=T03(1:3,3);

Z4=T04(1:3,3);

Z5=T05(1:3,3);

Z6=T06(1:3,3);

O1=T01(1:3,4);

O2=T02(1:3,4);

O3=T03(1:3,4);

O4=T04(1:3,4);

O5=T05(1:3,4);

O6=T06(1:3,4);

J1=simplify([cross(Z1,(O6-O1));Z1]);

J2=simplify([cross(Z2,(O6-O2));Z2]);

J3=simplify([cross(Z3,(O6-O3));Z3]);

J4=simplify([cross(Z4,(O6-O4));Z4]);

J5=simplify([cross(Z5,(O6-O5));Z5]);

J6=simplify([cross(Z6,(O6-O6));Z6]);

JALTER=[simplify(J1) simplify(J2) simplify(J3) simplify(J4) simplify(J5) simplify(J6)];

Theta1=0;
Theta2=0;
Theta3=0;
Theta4=0;
Theta5=0;
Theta6=0;

F = [20*9.81;20*9.81;20*9.81;110;110;60]

taw_static=transpose(eval(JALTER))*F
