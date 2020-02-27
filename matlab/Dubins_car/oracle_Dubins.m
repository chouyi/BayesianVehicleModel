% Input: obstacle, u_v, u_psi, v, psi, x, y
% Output: Yes/No

function [collision, reach_new] = oracle_Dubins(obs, u_v, u_psi, v, psi, x, y, pre_com)
    
    % Find a row number using omega and velocity
    num_row = find(pre_com(:,1)==u_v & pre_com(:,2)==u_psi & pre_com(:,3) <= v & pre_com(:,4) > v);
    max_x = max(max(pre_com))+30;

    % Compute reachable set
    xlimit = [pre_com(num_row, 5) pre_com(num_row, 6)];
    ylimit = [pre_com(num_row, 7) pre_com(num_row, 8)];
    xbox = xlimit([1 1 2 2 1]);
    ybox = ylimit([1 2 2 1 1]);
    
    xlimit2 = [-max_x max_x];
    ylimit2 = [xlimit2(1)-pre_com(num_row, 10) xlimit2(1)-pre_com(num_row, 9) xlimit2(2)-pre_com(num_row, 9) xlimit2(2)-pre_com(num_row, 10)];
    xbox2 = xlimit2([1 1 2 2 1]);
    ybox2 = ylimit2([1 2 3 4 1]);
    
    xlimit3 = [-max_x max_x];
    ylimit3 = [-xlimit3(1)+pre_com(num_row, 11) -xlimit3(1)+pre_com(num_row, 12) -xlimit3(2)+pre_com(num_row, 12) -xlimit3(2)+pre_com(num_row, 11)];
    xbox3 = xlimit3([1 1 2 2 1]);
    ybox3 = ylimit3([1 2 3 4 1]);
   
    [xi,yi] = polyxpoly(xbox,ybox,xbox2,ybox2);
    [xi2,yi2] = polyxpoly(xbox,ybox,xbox3,ybox3);
    reach_x = [xi(1:2);xi2(1:2);xi(3:4);rot90(xi2(3:4),2)]';
    reach_y = [yi(1:2);yi2(1:2);yi(3:4);rot90(yi2(3:4),2);]';
    reach_set = polyshape(reach_x, reach_y);
    
    % Change coordinates frame
    t = psi;
    
    Tz = [cos(t) -sin(t) 0 x; sin(t) cos(t) 0 y; 0 0 1 0; 0 0 0 1];

    for i = 1:size(reach_set.Vertices,1)
        new_coor = Tz * [reach_set.Vertices(i,:)'; 0; 1];
        reach_new(i,:) = [new_coor(1) new_coor(2)];
    end
    reach_new = polyshape(reach_new(:,1), reach_new(:,2));

    % Check intersection
    [in, out] = intersect(obs,reach_new);    
    if size(in.Vertices,1) ~= 0
        collision = 1;
    else
        collision = 0;
    end
end
% u_v = -0.9, u_psi = -1.8, v0 in 10, 10.5
% -1.90898<= x <= 0.597625
% -9.08219<= y <= -7.37291
% 5.72717<= x - y <= 9.41657
% -10.8496<= x + y <= -6.91687