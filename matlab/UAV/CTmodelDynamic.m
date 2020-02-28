function [next_x]=CTmodelDynamic(W,x,T)


    WT=W*T;  
    if (W==0)
    CT_A=[1 1 0 -1; 0 cos(WT) 0 -sin(WT);...
                0 1 1 1; 0 sin(WT) 0 cos(WT)];
    else
    CT_A=[1 sin(WT)/W 0 -(1-cos(WT))/W; 0 cos(WT) 0 -sin(WT);...
                0 (1-cos(WT))/W 1 sin(WT)/W; 0 sin(WT) 0 cos(WT)];
    end
    next_x=CT_A * x; 


end