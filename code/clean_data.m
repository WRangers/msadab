%% Picking out the fire data within Victoria
border=csvread('../data/victoria.csv');
blo_fire_ynum=cell(1,19);
for y=1:19
    year=2000+y;
    fstr=['../data/modis_',num2str(year),'_Australia.csv'];
    
    data=csvread(fstr,1,0);
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
        end
    end
    
    blo_fire_ynum{y}=blo_fire_num;
end

%% reshape the data
blo_num=11*18;
blo_fire_num_ally=zeros(blo_num,19);
for i=1:blo_num
    for y=1:19
        blo_fire_num_ally(i,y)=blo_fire_ynum{y}(i);
    end
end

%% plot the data
close all
figure
plot(2001:2019,blo_fire_num_ally);
ylim([-1 6]);
xlabel('Year');
ylabel('Fire Frequency');
matlab2tikz('../figures/firenum_blo_cha.tex');

%% print the data
fileid=fopen('../data/firenum_year_blo.txt','w');
formatSpec='%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n';
fprintf(fileid,formatSpec,blo_fire_num_ally);
fclose(fileid);