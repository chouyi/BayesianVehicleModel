clear;
close; 
M = readmatrix('./data/UAV_data.txt');
obs = readmatrix('./data/obstacles.txt');
pre_com = csvread('./data/reach-sets-5-2.csv');

idx.xy = 3;         idx.xx = 4;         idx.xz = 5;
idx.vy = 6;         idx.vx = 7;         idx.vz = 8;

sim_start = 4280; 
sim_end = 4290;

x = M(sim_start:sim_end, idx.xx);
y = M(sim_start:sim_end, idx.xy);
vx = M(sim_start:sim_end, idx.vx);
vy = M(sim_start:sim_end, idx.vy);

obs_coor_x = [obs(1) obs(1) obs(2) obs(2)];
obs_coor_y = [obs(3) obs(4) obs(4) obs(3)];
obs_zone = polyshape(obs_coor_x,obs_coor_y);

plot(x(2:end), y(2:end), '-ob');
hold on;
plot(obs_zone);
hold on;

for i = 1:sim_end-sim_start

    w = atan2(vy(i+1),vx(i+1))- atan2(vy(i), vx(i));
    if w > pi
        w = w - ceil(w/pi)*pi;
    elseif w < -pi
        w = w + ceil(w/pi)*pi;
    end

    [m, idx] = min(abs(pre_com(:,1) - w));
    w = pre_com(idx,1);

    collision = oracle(obs_zone, w, vx(i+1), vy(i+1), x(i+1), y(i+1), pre_com);
    
    if collision==1
        fprintf('Collision: Yes\n'); 
    else
        fprintf('Collision: No\n'); 
    end
end
