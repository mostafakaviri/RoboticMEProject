function a = animation( link1,link2,link3,link4,link5 )
a = 0;

        
        plot([link1(1) , link2(1)] , [link1(3) link2(3)],[link2(1) , link3(1)] , [link2(3) link3(3)],[link3(1) , link4(1)] , [link3(3) link4(3)],[link4(1) , link5(1)] , [link4(3) link5(3)])
       
        hold on
        xlim ([-2, 2]);
        ylim ([-2, 2]);
        
        circle(link1(1),link1(3), .01)
        
        circle(link2(1),link2(3), .01)
        
        circle(link3(1),link3(3), .01)
        
        circle(link4(1),link4(3), .01)
        
        circle(link5(1),link5(3), .01)
        
        
        
        pause(0.0000001)
        
end

