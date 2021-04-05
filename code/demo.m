% load blo_fire_data.mat
ovo = [6,12,1;
    6,3,2;
    5,2,2;
    5,8,2;
    9,12,1;
    5,13,2
    ];
for iii = 1:size(ovo,1)
    DATAS = blo_fire_data_v2{ovo(iii,1), ovo(iii,2)};
    A = DATAS{3,ovo(iii,3)};
    cruise_line
    fstr=['../figures/example',num2str(iii)];
    print(fstr,'-dpdf','-r600');
end


