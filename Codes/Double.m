function q_out=Double(q1,q2) 
    d2=2*(q1*q2');
    q_out=d2*(q2-q1);
end