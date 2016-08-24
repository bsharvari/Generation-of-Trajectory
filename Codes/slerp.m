% Spherical linear interpolation between two quaternions
% Returns a quaternion, takes in two quaternions and value of u (step-size)
function quat=slerp(q1, q2, u)
    cosine=q1*q2';
    sine=sqrt(1-cosine*cosine);
    theta=atan2(sine,cosine);
    
    if (sine<0.0001)
        a1=1-u;
        a2=u;
    else
        a1=sin((1-u)*theta)/sine;
        a2=sin(u*theta)/sine;
    end
    
    quat=(a1*q1)+(a2*q2);
end