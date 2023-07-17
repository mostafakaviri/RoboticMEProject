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