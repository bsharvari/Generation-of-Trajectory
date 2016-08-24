% Convert quaternion to rotation matrix
% Takes in a quaternion, returns corresponding rotation matrix

function R=Q_to_R(Q)
    % n-column
    R(1,1)=1-(2*(Q(2)*Q(2)+Q(3)*Q(3)));
    R(2,1)=2*(Q(1)*Q(2)+Q(4)*Q(3));
    R(3,1)=2*(Q(1)*Q(3)-Q(4)*Q(2));
    
    % o-column
    R(1,2)=2*(Q(1)*Q(2)-Q(4)*Q(3));
    R(2,2)=1-(2*(Q(1)*Q(1)+Q(3)*Q(3)));
    R(3,2)=2*(Q(2)*Q(3)+Q(4)*Q(1));
    
    % a-column
    R(1,3)=2*(Q(1)*Q(3)+Q(4)*Q(2));
    R(2,3)=2*(Q(2)*Q(3)-Q(4)*Q(1));
    R(3,3)=1-(2*(Q(1)*Q(1)+Q(2)*Q(2)));
end