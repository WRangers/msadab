% 尝试使用ga算法求解这个优化问题

%% 读取数据
% load 12_12_area_fireData
% load blo_fire_data_V2.mat

% DATAS = blo_fire_data_v2{4,1};
% A = DATAS{3}; % T0 = 100 才有�?
% 
% DATAS = blo_fire_data_v2{5,12};
% A = DATAS{3,1};

lat = A(:,1);
lon = A(:,2);
scan = A(:,4);
track = A(:,5);
global fireNum firePos vmax x y shortestpath t T0

fireNum = length(lat); % �?火点的个�?
firePos = [lon, lat]; % �?火点的经纬度
vmax = 1.2; % 无人机的飞行速度，km/min
T0 = 100; % �?大观察间�? 30 min
t = zeros(fireNum, 1); % 起飞点到达各个着火点的路程时�? @fit_ness
%%
y = zeros(fireNum, fireNum); % �?火点的路程时�?
for i = 1:fireNum
    for j = 1:fireNum
        if i == j
            y(i,j) = inf;
        else
            y(i,j) = spherediff(firePos(i,:), firePos(j,:));
        end
    end
end
y = y./ vmax;

x = (scan+track)*2/vmax; % �?火点的停留时�? 

% 各个�?火点为起点的�?短路�?
shortestpath = zeros(fireNum, fireNum);
for i = 1:fireNum
    shortestpath(:,i) = mydijkstra(i);
end
% cankao = y;
% for i = 1:fireNum
%     cankao(i,:) = cankao(i,:) + x';
% end
% disp('各个�?火点为起点的�?短路径：'); disp(shortestpath);
% disp("�?短距离的参�?�矩�?:"); disp(cankao);

posT0 = [mean(lon), mean(lat)]; % 实践证明 均�?�是�?接近�?
fun = @fitness;
LB = [min(lon), min(lat)];
UB = [max(lon), max(lat)];
T = ga(fun, 2, [], [], [], [], LB, UB);

%%
name = 'opentopomap';
url = 'a.tile.opentopomap.org';
displayName = 'Open Topo Map';
addCustomBasemap(name,url,'DisplayName',displayName,...
                'Attribution','Victoria, Austria');

%%
% % 显示火灾的分�?
% disp('�?火点的停留时�?:'); disp(x');
% disp('�?火点之间的路程时间：'); disp(y);
figure();
% geoplot(lat, lon, 'rs', 'markerfacecolor', 'r');
geoplot(lat, lon, 'r.');
hold on
for i = 1:fireNum
    text(lat(i), lon(i), string(i));
    sizetmp = (x(i) - min(x))/(max(x)-min(x))*6 + 6;
    if isnan(sizetmp)
        sizetmp = 6;
    end
    sizetmp = ceil(sizetmp);
    geoplot(lat(i), lon(i), 'marker', 's', 'markersize', sizetmp, ...
        'color', 'r', 'markerfacecolor', 'r');
end
%% ========================================
% geobasemap colorterrain 
geobasemap('opentopomap')
%% ========================================
% T = [mean(lon), mean(lat)]; % 起飞点位置设�?
[N, minid, vv] = fitness2(T);

% 绘制轨迹
if N == inf
    disp('无解！尝试修改T0或舍去距离较远的�?')
else
    hold on;
    geoplot(T(2), T(1), 'bs', 'markerfacecolor', 'b');

   
    points = find(vv==1);
    mypath = shortestpath(:, minid);
    for i = 1:length(points)
        if i < length(points)
            tmp = points(i):max(points(i+1)-1, points(i));
        else
            tmp = points(i):length(mypath);
        end
        dpnts = firePos(mypath(tmp), :);
        geoplot([T(2), dpnts(:,2)', T(2)], [T(1), dpnts(:,1)', T(1)]);
    end
%     disp('EOC到达各个顶点的�??:')
%     disp(t');
%     disp('选择的路�?:')
%     disp(shortestpath(:,minid)')
end
disp('SSA的个�?'); disp(N)
NNN = zhongjiNum;
disp('中继器的个数'); disp(NNN);
%% 根据经纬度计算两点之间的距离
function distance = spherediff(pa, pb)
% (经度, 纬度) (经度, 纬度)
R = 6371; % km
havLat = sin(deg2rad(pa(2)-pb(2)) / 2);
havLon = sin(deg2rad(pa(1)-pb(1)) / 2);
a = havLat*havLat + cos(deg2rad(pa(2)))*cos(deg2rad(pb(2)))*havLon*havLon;
distance = 2*R*atan2(sqrt(a),sqrt(1-a));
end

%% 根据Dijkstra算法求出某个起始点开始遍历所有点的最短路�?
function mypath = mydijkstra(sa)
% sa 起始点的序号，shortestpath 是最短路�?
global x y fireNum
vst = zeros(fireNum, 1); vst(sa) = 1;
mypath = zeros(fireNum, 1); mypath(1) = sa;
for i = 2:fireNum
    tmpCol = y(mypath(i-1), :) + x';
    tmp = min(tmpCol(vst == 0));
    tmpids = find(tmpCol == tmp);
    vst(tmpids(1)) = 1;
    mypath(i) = tmpids(1);
end
% disp("验证:"); disp(vst');
% disp("路径:"); disp(mypath');
end

%% ga的优化函�?
% function [N, minid, vv] = fitness(T)
function TS = fitness(T)
% T 是起飞点的经纬度
% 返回 无人机最少的无人机N
global fireNum x y firePos shortestpath t vmax T0

for i = 1:fireNum % 计算起飞点到达各个着火点的路程时�?
    t(i) = spherediff(T, firePos(i,:))/vmax;
end

Ns = ones(1, fireNum);
Ts = zeros(1, fireNum);
vst = zeros(fireNum, fireNum); vst(1,:) = ones(1, fireNum);

for i = 1:fireNum
    T = T0; % 剩余飞行时间
    curPath = shortestpath(:,i); % 待飞点的�?短路�?
    go = t(curPath(1)); back = t(curPath(1)); stop = y(curPath(1));
    if T0-go-back-stop < 0
        Ns(i) = inf;
        Ts(i) = inf;
    else
        j = 1;
        while j <= fireNum
            T = T-go-back-stop;
            Ts(i) = Ts(i) + go + stop;
            
            if T >= 0
                j = j + 1;
                
                T = T + back;
                % 下一个点的数�?
                if j <= fireNum
                    go = y(curPath(j-1),curPath(j));
                    back = t(curPath(j));
                    stop = x(curPath(j));
                end
            else
                vst(j,i) = 1;
                go = t(curPath(j));
                back = t(curPath(j));
                stop = x(curPath(j));
                T = T0;
                Ts(i) = Ts(i) + back;
                Ns(i) = Ns(i) + 1;
                if T0-go-back-stop < 0 % 无论如何都到达不了这个点（从EOC起飞后回不来�?
                    Ns(i) = fireNum + 10;
                    break;
                end
            end
        end
    end
end


N = min(Ns);
minid = find(Ts == min(Ts(Ns == min(Ns))));
TS = Ts(minid(1));
% vv = vst(:, minid);

end

%% 返回�?小的无人�?
function [N, minid, vv] = fitness2(T)
% T 是起飞点的经纬度
% 返回 无人机最少的无人机N
global fireNum x y firePos shortestpath t vmax T0

for i = 1:fireNum % 计算起飞点到达各个着火点的路程时�?
    t(i) = spherediff(T, firePos(i,:))/vmax;
end

Ns = ones(1, fireNum);
Ts = zeros(1, fireNum);
vst = zeros(fireNum, fireNum); vst(1,:) = ones(1, fireNum);

for i = 1:fireNum
    T = T0; % 剩余飞行时间
    curPath = shortestpath(:,i); % 待飞点的�?短路�?
    go = t(curPath(1)); back = t(curPath(1)); stop = y(curPath(1));
    if T0-go-back-stop < 0
        Ns(i) = inf;
        Ts(i) = inf;
    else
        j = 1;
        while j <= fireNum
            T = T-go-back-stop;
            Ts(i) = Ts(i) + go + stop;
            
            if T >= 0
                j = j + 1;
                
                T = T + back;
                % 下一个点的数�?
                if j <= fireNum
                    go = y(curPath(j-1),curPath(j));
                    back = t(curPath(j));
                    stop = x(curPath(j));
                end
            else
                vst(j,i) = 1;
                go = t(curPath(j));
                back = t(curPath(j));
                stop = x(curPath(j));
                T = T0;
                Ts(i) = Ts(i) + back;
                Ns(i) = Ns(i) + 1;
                if T0-go-back-stop < 0 % 无论如何都到达不了这个点（从EOC起飞后回不来�?
                    Ns(i) = fireNum + 10;
                    break;
                end
            end
        end
    end
end


N = min(Ns);
minid = find(Ts == min(Ts(Ns == min(Ns))));
minid = minid(1);
vv = vst(:, minid);

end

%% 计算中继器的个数
function n = zhongjiNum
global y t
Dis = y;
for i = 1:size(Dis,1)
    for j = 1:size(Dis,2)
        if (i == j)
         Dis(i,j) = 0;
        end
    end
end

Rmax = max([max(max(Dis)),t']);
R0 = 20;
if Rmax <= R0/sqrt(3)
    n = 1;
elseif Rmax <= R0*sqrt(3)
    n = 2;
else
    n = ceil((2/3 - sqrt(3)/(2*pi))*(Rmax/R0)^2)+1;
end

end
