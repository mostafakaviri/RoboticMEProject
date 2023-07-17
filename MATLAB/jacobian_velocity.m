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