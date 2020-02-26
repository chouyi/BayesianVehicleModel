% Input: obstacle, omega, vx, vy, x, y
% Output: Yes/No

function collision = oracle(obs, w, vx, vy, x, y, pre_com)
    
    % Find a row number using omega and velocity
    v = sqrt(vx^2+vy^2);
    

    num_row = find(pre_com(:,1)==w & pre_com(:,2) <= v & pre_com(:,3) > v);
    max_x = 50;

    % Compute reachable set
    xlimit = [pre_com(num_row, 4) pre_com(num_row, 5)];
    ylimit = [pre_com(num_row, 6) pre_com(num_row, 7)];
    xbox = xlimit([1 1 2 2 1]);
    ybox = ylimit([1 2 2 1 1]);
    
    xlimit2 = [-max_x max_x];
    ylimit2 = [xlimit2(1)-pre_com(num_row, 9) xlimit2(1)-pre_com(num_row, 8) xlimit2(2)-pre_com(num_row, 8) xlimit2(2)-pre_com(num_row, 9)];
    xbox2 = xlimit2([1 1 2 2 1]);
    ybox2 = ylimit2([1 2 3 4 1]);
    
    xlimit3 = [-max_x max_x];
    ylimit3 = [-xlimit3(1)+pre_com(num_row, 10) -xlimit3(1)+pre_com(num_row, 11) -xlimit3(2)+pre_com(num_row, 11) -xlimit3(2)+pre_com(num_row, 10)];
    xbox3 = xlimit3([1 1 2 2 1]);
    ybox3 = ylimit3([1 2 3 4 1]);
   
    [xi,yi] = polyxpoly(xbox,ybox,xbox2,ybox2);
    [xi2,yi2] = polyxpoly(xbox,ybox,xbox3,ybox3);
    reach_x = [xi(1:2);xi2(1:2);xi(3:4);rot90(xi2(3:4),2)]';
    reach_y = [yi(1:2);yi2(1:2);yi(3:4);rot90(yi2(3:4),2);]';
    reach_set = polyshape(reach_x, reach_y);
    
    % Change coordinates frame
    t = atan2(vy,vx);
    
    Tz = [cos(t) -sin(t) 0 x; sin(t) cos(t) 0 y; 0 0 1 0; 0 0 0 1];

    for i = 1:size(reach_set.Vertices,1)
        new_coor = Tz * [reach_set.Vertices(i,:)'; 0; 1];
        reach_new(i,:) = [new_coor(1) new_coor(2)];
    end
    reach_new = polyshape(reach_new(:,1), reach_new(:,2));
    plot(reach_new); hold on;

    % Check intersection
    [in, out] = intersect(obs,reach_new);    
    if size(in.Vertices,1) ~= 0
        collision = 1;
    else
        collision = 0;
    end
end
% 
% omega: -0.171
% vx0 : 8 , 8.5
% 15.3673<= x <= 16.9797
% -3.85359<= y <= -1.69521
% 17.5254<= x - y <= 20.3703
% 11.5137<= x + y <= 15.2845