% Convert rotation matrix to quaternion
% Takes in rotation matrix, generates quaternion
% Caution: May not pick quaternion you want because have multiple solutions

function Q=R_to_Q(R)
    dx=1+R(1,1)-R(2,2)-R(3,3);
    dy=1-R(1,1)+R(2,2)-R(3,3);
    dz=1-R(1,1)-R(2,2)+R(3,3);
    ds=1+R(1,1)+R(2,2)+R(3,3);
    
    vector=[dx,dy,dz,ds];
    [~,index]=max(vector);
    
    switch(index)
        case 1
            vx=sqrt(dx)/2;
            vy=(R(1,2)+R(2,1))/(4*vx);
            vz=(R(1,3)+R(3,1))/(4*vx);
            s=(R(3,2)-R(2,3))/(4*vx);
        case 2
            vy=sqrt(dy)/2;
            vx=(R(1,2)+R(2,1))/(4*vy);
            vz=(R(2,3)+R(3,2))/(4*vy);
            s=(R(1,3)-R(3,1))/(4*vy);
        case 3
            vz=sqrt(dz)/2;
            vx=(R(1,3)+R(3,1))/(4*vz);
            vy=(R(2,3)+R(3,2))/(4*vz);
            s=(R(2,1)-R(1,2))/(4*vz);
        case 4
            s=sqrt(ds)/2;
            vx=(R(3,2)-R(2,3))/(4*s);
            vy=(R(1,3)-R(3,1))/(4*s);
            vz=(R(2,1)-R(1,2))/(4*s);
    end
    
    Q=[vx, vy, vz, s];
end