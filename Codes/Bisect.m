function q_out=Bisect(q1,q2)
    q_out=q1+q2;
    d=sqrt(q_out*q_out');
    q_out=q_out./d;
end