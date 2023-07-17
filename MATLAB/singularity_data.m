clear;clc;close all;
syms Theta2 Theta3
j = 1;
for Theta2=-pi:0.1:-0.261
    for Theta3=-pi:0.1:1.35
        x=(15*sin(Theta3) - 15*sin(Theta2) + 15*cos(Theta3)^2*sin(Theta2) + 79*cos(Theta2)*sin(Theta3) + 15*cos(Theta2)*cos(Theta3)*sin(Theta3));
        
        if (x<=0.0001)
            s_data(j,1) = Theta2;
            s_data(j,2) = Theta3;
            j = j+1;

        end
    end
end
writematrix(s_data,'s_data_tab.txt','Delimiter','tab')