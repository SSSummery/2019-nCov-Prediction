%name:visualization(figure 4)
%generate figure 4 to visualize the to and alpha0
clear
clf
close all
clc
[k,tt]=meshgrid(0.29:0.001:0.35,20:0.01:22);
z=zeros(length(tt(:,1)),length(k(1,:)));
for k_index=1:length(k(1,:))
    for t_index=1:length(tt(:,1))
        z(t_index,k_index)=myFitCal(k(1,k_index),tt(t_index,1));
    end
end
clf
surfl(k,tt,z);
hold on;
mini=min(min(z));
[x,y]=find(z==mini);
plot3(k(1,y),tt(x,1),z(x,y),'r.','markersize',30)
shading interp;
material dull;
light('color',[1 1 1],'position',[0.35 22 250000000],'style','local') 
lighting gouraud;
text(k(1,y),tt(x,1),z(x,y),['(',num2str(tt(x,1)),',',num2str(k(1,y)),')'],'color','k')
title('{\itInversion Parameters Analysis}')
ylabel('days')
