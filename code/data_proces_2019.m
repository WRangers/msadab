%% Picking out the fire data within Victoria
data=csvread('../data/modis_2019_Australia.csv',1,0);

border=csvread('../data/victoria.csv');
xv=border(:,2);
yv=border(:,1);
xq=data(:,2);
yq=data(:,1);
[in, on]=inpolygon(xq,yq,xv,yv);
data=data(in,:);

%% confidence select
label=data(:,7)>=50;
data=data(label,:);

%% numbering
a=141;
b=-34;
data(:,8)=ceil((data(:,2)-a)*2); % block_x
data(:,9)=ceil(abs(data(:,1)-b)*2); % block_y
data(:,10)=11*(data(:,8)-1)+data(:,9); % block_num

%% divide the fires into blocks
bol=cell(11,18);
bol_tot=11*18;
for i=1:bol_tot
    label=data(:,10)==i;
    bol{i}=data(label,:);
end

%% sort out the fire of bolocks
blo_fire_data=cell(11,18);
blo_fire_num=zeros(11,18);
for k=1:bol_tot
    if bol{k}
        D=sortrows(bol{k},6);
        aqc_time=D(:,6);
        cut=zeros(length(aqc_time),1);
        for i=1:length(aqc_time)-1
           if abs(aqc_time(i+1)-aqc_time(i))>=100
              cut(i+1)=1; 
           end
        end
        id=cumsum(cut)+1;
        group_num=id(length(id));
        fire_data=cell(3,group_num);
        blo_fire_num(k)=group_num;
        for i=1:group_num
            label=id==i;
            a_fire_data=D(label,:);
            fire_data{3,i}=a_fire_data;
            fire_data{1,i}=mean(a_fire_data(:,1));
            fire_data{2,i}=mean(a_fire_data(:,2));
        end
        blo_fire_data{k}=fire_data;
    end
end

blo_fire_freq=blo_fire_num./sum(sum(blo_fire_num));

%% fire place
fire_info=[];
for i=1:bol_tot
    if ~isempty(blo_fire_data{i})
        fire_p=cell2mat(blo_fire_data{i}(1:2,:));
        fire_i=blo_fire_data{i}(3,:);
        n=length(fire_i);
        for f=1:n
            fire=cell2mat(fire_i(f));
            mean_area=mean(fire(:,4).*fire(:,5));
            mean_bri=mean(fire(:,3));
            fire_p(3,f)=mean_area;
            fire_p(4,f)=mean_bri;
        end
        fire_info=[fire_info,fire_p];
    end
end

%% plot the blocks and fire place
name = 'opentopomap';
url = 'a.tile.opentopomap.org';
displayName = 'Open Topo Map';
addCustomBasemap(name,url,'DisplayName',displayName,...
                'Attribution','Victoria, Austria');

close all
gx = geoaxes;
len=length(fire_info(4,:));
C=[fire_info(4,:)/max(fire_info(4,:));zeros(1,len);zeros(1,len)]';
geoscatter(gx,fire_info(1,:),fire_info(2,:),...
              fire_info(3,:)*30,C,'s','filled');
geobasemap('opentopomap')
% gx.LongtitudeLabel.String='Longtitude';
% gx.LatitudeLabel.String='Latitude';
% gx.Scalebar.Visible = 'on';
geolimits([-39.5 -34],[141 150]);
hold on;
lon=linspace(141,150,19);
lat=linspace(-34,-39.5,12);
for i=lon
    geoplot([-34 -39.5 ],[i i],'k');
end
for i=lat
    geoplot([i i],[141 150],'k');
end
hold off;
print('../figures/blocks','-dpdf','-r600');

%% fire severe map
close all
figure
geodensityplot(data(:,1),data(:,2),data(:,3),'FaceColor','interp');
geolimits([-39.5 -34],[141 150]);
print('../figures/severmap','-dpdf','-r600');

%% hotmap
close all
figure
heatmap(blo_fire_freq);
colormap(mycmap);
% set(gcf,'unit','centimeters','position',[0 0 28 16]);
print('../figures/hotmap','-dpdf','-r600');

%% a fire example
firexam=blo_fire_data{6,12}{3,1};
close all;
figure;
A=firexam(:,4).*firexam(:,5);
B=firexam(:,3);
len=length(B);
C=[B/max(B),zeros(len,1),zeros(len,1)];
geoscatter(firexam(:,1),firexam(:,2),...
              A*30,C,'s','filled');
geobasemap('opentopomap');
print('../figures/firexam','-dpdf','-r600');

%% visualization of ssa and repeaters
close all;
figure
cdata = SSA;
heatmap(cdata);
% axes euqal
set(gcf,'unit','centimeters','position',[0 0 20 10]);
print('../figures/ssan','-dpdf','-r600');

figure
cdata =roundn(ZJ,-3)
heatmap(cdata);
% axes euqal
set(gcf,'unit','centimeters','position',[0 0 20 10]);
print('../figures/repn','-dpdf','-r600');

%% calculate total number of drones to purchase
Ns=2*ceil(sum(sum(blo_fire_freq.*SSA)))
Nr=2*ceil(sum(sum(blo_fire_freq.*ZJ)))

%% 
close all;
name = 'opentopomap';
url = 'a.tile.opentopomap.org';
displayName = 'Open Topo Map';
addCustomBasemap(name,url,'DisplayName',displayName)
webmap opentopomap(opentopomap)
wmmarker(firexam(:,1),firexam(:,2))
% wmlimits([-33, -40], [140, 150])

%% 
figure
geoscatter(fire_origin(1,:),fire_origin(2,:));
geobasemap('opentopomap')