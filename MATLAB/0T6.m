


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
