clear all
close all 
clc

for time_ind = 0:800
    filename = ['ExportCellProp_' num2str(time_ind) '.txt'];
    display(filename);
    temp = exist(filename);
    while(temp==0)
       pause(2);
       temp = exist(filename);
    end
   main_signaling(time_ind); 
end
   



