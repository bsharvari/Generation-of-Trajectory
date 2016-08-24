% Project 1

% Robot

clc;
clear all;

file_name_robot = 'robot.key';
R = importdata(file_name_robot);
num_joints = 6;                                                            % Change the value of num_joints according to the problem statement
num_frames = R(1,2);
key_frames = R(1,1);

Pos_r = zeros(key_frames,num_joints);
D = zeros(key_frames,num_joints);

% Get the positions from the robot.key file
x = 0;
for n = 2 : num_joints : (size(R,1)-(num_joints-1))
    x = x + 1;
    p = 0;
    for i = n : n + ((num_joints/2)-1)
        for j = 1 : 2
            p = p + 1;
            Pos_r(x,p) = R(i,j);
        end
    end
end

% Get the velocities from the robot.key file
x = 0;
for n = (num_joints - 1) : num_joints : (size(R,1)-2)
    x = x + 1;
    p = 0;
    for i = n : n + ((num_joints/2)-1)
        for j = 1 : 2
            p = p + 1;
            D(x,p) = R(i,j);
        end
    end
end

% Compute the spline coefficients using hermite interpolation

% Consider the standard hermite spline equations 
% h1 = 2*u_3 -3*u_2 + 1;
% h2 = -2*u_3 +3*u_2;
% h3 = u_3 - 2*u_2 + u;
% h4 = u_3 - u_2;

spline_r = zeros(((key_frames-1)*num_joints),4);
s = 0;
for i = 1 : (key_frames-1)
    for j = 1 : num_joints
        s = s + 1;
        spline_r(s,4) = 2*Pos_r(i,j) - 2*Pos_r(i+1,j) + D(i,j) + D(i+1,j);
        spline_r(s,3) = -3*Pos_r(i,j) + 3*Pos_r(i+1,j) - 2*D(i,j) - D(i+1,j);
        spline_r(s,2) = D(i,j);
        spline_r(s,1) = Pos_r(i,j);
    end
end

% Compute the intermediate frames
inter_r = zeros(num_frames,num_joints);
delta_u = (key_frames - 1)/(num_frames - 1);
k = 0;
l = 0;

for u = 0 : delta_u : (key_frames-1)
    k = k + 1;    
    y = u;
    int_u = fix(u);                                                        % Scale u so that it is between 0 and 1 everytime 
    y = y - int_u;
    if int_u == key_frames - 1;
        j = size(spline_r,1)-(num_joints-1);
        l = l+1;                                                           % To get the last frame
        inter_r (k,l:l+(num_joints-1)) = spline_r(j:j+(num_joints-1),4)*y^3 + ...
            spline_r(j:j+(num_joints-1),3)*y^2 + spline_r(j:j+(num_joints-1),2)*y + ...
            spline_r(j:j+(num_joints-1),1);
        break;
    end
    for j = (1+int_u*num_joints) : (1+int_u)*num_joints                    % To get all the frames but the last
        l = l + 1;
        inter_r (k,l) = spline_r(j,4)*y^3 + spline_r(j,3)*y^2 + ...
            spline_r(j,2)*y + spline_r(j,1);
        if l == num_joints
            l = 0;
        end
    end
end

% Write to the robot.ang file

folder_robot = 'robot';
dlmwrite(fullfile(folder_robot,'robot.ang'),num_frames);
dlmwrite(fullfile(folder_robot,'robot.ang'),inter_r,'-append','delimiter','\t','precision',10);

% Object

file_name_object = 'object.key';
O = importdata(file_name_object);
T_o = zeros((key_frames*3),4);
x = 1;
y = 1;

% Get data from the file in a matrix
for n = 2 : num_joints : (size(O,1)-(num_joints-1))
    for i = n : n + (num_joints-1)
        for j = 1 : 2
            T_o(x,y) = O(i,j);
            y = y + 1;
            if y == (num_joints-1)
                x = x + 1;
                y = 1;
            end
        end
    end
end

% Get the positions from the transformation matrix
Pos_o = zeros(key_frames,3);
count = 1;
j = 1;
for  i = 1 : size(T_o,1)
    if count == 4
        count = 1;
        j = j + 1;
    end
    Pos_o(j,count) = T_o(i,4);
    count = count + 1;
end

% Compute the spline coefficients using hermite interpolation

% x = a1*u + a0;
% y = b1*u + b0;
% z = c1*u + c0;

spline_o = zeros(4,6);
for  i = 1 : 4
    spline_o(i,:) = [Pos_o(i,1) (Pos_o(i+1,1)-Pos_o(i,1)) Pos_o(i,2)...
            (Pos_o(i+1,2)-Pos_o(i,2)) Pos_o(i,3) (Pos_o(i+1,3)-Pos_o(i,3))];
end

% Compute the intermediste positions of the object
inter_p = zeros(3*num_frames,1);
k = 0;
for u = 0 : delta_u : (key_frames - 1)
    y = u;
    x = fix(u);
    y = y - x;
    if x == key_frames - 1                                                 % To get the last frame 
        for j = 1 : 2 : 5
            k = k + 1;
            inter_p(k,1) = spline_o(x,j) + spline_o(x,j+1)*y;
        end
        break;
    end
    for j = 1 : 2 : 5                                                      % To get all the frames but the last
        k = k + 1;
        inter_p(k,1) = spline_o(x+1,j) + spline_o(x+1,j+1)*y;
    end
end

% Convert the rotation matrices to equivalent quaternions
temp = zeros(3,4);
Q = zeros(key_frames,4);
k = 0;
for  i = 1 : 3 : (size(T_o,1)-2)
    k = k + 1;
    l = 0;
    for j = i : (i+2)
        l = l + 1;
        temp(l,:) = T_o(j,:);
    end
    Q(k,:) = R_to_Q(temp);
end

% Check whether the quaternion is correct
k = 0;
for i = 1 : 3 : (size(T_o,1)-2)
    x = 0;
    k = k + 1;
    l = 0;
    for j = i : (i+2)
        l = l + 1;
        x = x + T_o(j,4)*Q(k,l);
    end
    if x < 0
        Q(k,:) = [-1 -1 -1 1].*Q(k,:);
    end
end

% Calculation of intermediate frames
inter_q = Q_interpolation(Q,O(1,1),O(1,2));

% Convert the quaternions to respective rotation matrices
temp1 = zeros(3,3);
inter_o = zeros((num_frames*3),3);
l = 0;
for i = 2:(num_frames-1)
   temp1 = Q_to_R(inter_q(i,:));
   for j = 1:3
       l = l + 1;
       for k = 1:3
           inter_o(l,k) = temp1(j,k);
       end
   end
end

% Concatenate the intermediate rotation matrices and the position vectors
inter_object = horzcat(inter_o,inter_p);

% Write to the object.traj file

folder_object = 'object';
dlmwrite(fullfile(folder_object,'object.traj'),num_frames);
dlmwrite(fullfile(folder_object,'object.traj'),inter_object,'-append','delimiter','\t','precision',10);


