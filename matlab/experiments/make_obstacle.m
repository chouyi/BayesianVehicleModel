function obs_zone = make_obstacle(obs_center, obs_size, dist_min, dist_max)
    size = 2/obs_size;
    
    obs_coor = [-size -size size size; -size size size -size];
    
    t = -pi() + (pi()+pi())*rand(1,1);
    
    if dist_max == 0
        trans = -(size-0.1) + ((size-0.1)+(size-0.1))*rand(2,1);
    else
        trans = randi([dist_min dist_max],2,1);
    end
    
    for i = 1:4 
        Tz = [cos(t) -sin(t) 0 trans(1); sin(t) cos(t) 0 trans(2); 0 0 1 0; 0 0 0 1];
        new_coor = Tz * [obs_coor(1, i); obs_coor(2, i); 0; 1];
        obs_mod(i,:) = [new_coor(1) new_coor(2)];
    end
    obs_zone = polyshape(obs_mod(:,1)+obs_center(1), obs_mod(:,2)+obs_center(2));
end