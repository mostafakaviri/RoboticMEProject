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

%% singular t3 t2
simplify(15*sin(Theta3) - 15*sin(Theta2) + 15*cos(Theta3)^2*sin(Theta2) + 79*cos(Theta2)*sin(Theta3) + 15*cos(Theta2)*cos(Theta3)*sin(The