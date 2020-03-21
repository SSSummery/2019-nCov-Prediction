# 2019-nCov-Prediction
#This Repository contains the code used in the essay ‘Prediction of 2019-nCov in Italy based on PSO and #inversion analysis’

#Summer here, I add the myFitCal.m/mypso.m/visualization.m/differentalpha.m code here.

%%%%%%%%%%%%%%%%%%%%%%%%%%begin%%%%%%%%%begin%%%%%%%%begin%%%%%%%%begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%name:myFitCal.m
%This code definite the fitness functions for PSO during the first and second period of our work. You may %need some time to figure out the note in the code, cause some plot command should be shield during PSO. %But if you want to see the detailed information and value, especially you want to get the Figure 5-6.

%Runge-Kutta 
%The Runge-Kutta iteration is in this function. By adopting this iteration, we can finally get the %prediction value and use them to get the LSE as fitness for PSO.
%%
function fitness = myFitCal(k,time) %period 1(2 parameters)
%function fitness = myFitCal(dk) %period 2(1 parameter)
actual=[1492;1247;778;769;587;667;146;561;240;238;250;78;93;97;53;62;14];
%3590;3497;2547;2651;2313;977;1797;
ac=flipud(actual);
acc=zeros(length(ac),1);
acc(1)=ac(1);
for i=1:length(ac)-1
    acc(i+1)=acc(i)+ac(i+1);
end
%%
c = 0.95;
d = 0.05;
e = 1 - d;
h = 0.1;
n =450;
t=0:h:n*h;
A = zeros(6,n+1);
A(:,1)= [60480000 1 0 0 0 0];
%%
for i = 1:n
%     if (i<370)
%         k=0.3166;
%     else
%         k=dk;
%     end
    s1sum=0;
    s2sum=0;
    s3sum=0;
    for j=1:i
            s1sum=s1sum+0.4812./(j*h).*exp(-(log(j*h) - log(4)).^2/1.3747)*A(2,i-j+1)*h;
            s2sum=s2sum+0.6847./(j*h).*exp(-(log(j*h) - log(10)).^2/0.6790)*A(3,i-j+1)*h;
            s3sum=s3sum+1.7469./(j*h).*exp(-(log(j*h) - log(12)).^2/0.1043)*A(4,i-j+1)*h;

    end
    %%
    y=A(:,i);
    k1 =zeros(6,1);%初始化列向量
    k1(1) = -y(1)*(y(2) + y(3))*k/(y(2)+y(3)+y(1));
    k1(2) = y(1)*(y(2)+y(3))*k/(y(2)+y(3)+y(1))-s1sum;
    k1(3) = s1sum-s2sum;
    k1(4) = c*s2sum-s3sum;
    k1(5) = d*s3sum + (1-c)*s2sum;
    k1(6) = e*s3sum;  
    %%
    y=A(:,i)+h*k1/2;
    k2=zeros(6,1);
    k2(1) = -y(1)*(y(2) + y(3))*k/(y(2)+y(3)+y(1));
    k2(2) = y(1)*(y(2)+y(3))*k/(y(2)+y(3)+y(1))-s1sum;
    k2(3) = s1sum-s2sum;
    k2(4) = c*s2sum-s3sum;
    k2(5) = d*s3sum + (1-c)*s2sum;
    k2(6) = e*s3sum;
    %%
    y=A(:,i)+h*k2/2;
    k3=zeros(6,1);
    k3(1) = -y(1)*(y(2) + y(3))*k/(y(2)+y(3)+y(1));
    k3(2) = y(1)*(y(2)+y(3))*k/(y(2)+y(3)+y(1))-s1sum;
    k3(3) = s1sum-s2sum;
    k3(4) = c*s2sum-s3sum;
    k3(5) = d*s3sum + (1-c)*s2sum;
    k3(6) = e*s3sum;
    %%
    y=A(:,i)+h*k3;
    k4=zeros(6,1);
    k4(1) = -y(1)*(y(2) + y(3))*k/(y(2)+y(3)+y(1));
    k4(2) = y(1)*(y(2)+y(3))*k/(y(2)+y(3)+y(1))-s1sum;
    k4(3) = s1sum-s2sum;
    k4(4) = c*s2sum-s3sum;
    k4(5) = d*s3sum + (1-c)*s2sum;
    k4(6) = e*s3sum;
    A(:,i+1)=A(:,i)+h*(k1+2*k2+2*k3+k4)/6; 
end
fitness=0;
%%
%choice 1 2.22-3.9
for dex=1:length(acc)
    ppp=10*(dex+time-1)+1;
    fitness=fitness+(A(4,ppp)-acc(dex))^2;
end
% plot(0:0.1:21+length(acc)-1,A(4,1:(21+length(acc)-1)*10+1),'-.','color','b');
% hold on
% plot(21:1:length(acc)+20,acc,'*','color','r');
% hold on
% plot(t,A(4,:));
%%
% plot(20:0.1:21+length(acc)-1,A(4,201:(21+length(acc)-1)*10+1),'-.','color','b');
% hold on
% plot(21:1:length(acc)+20,acc,'*','color','r');
% hold on
% plot(20:0.1:n*h,A(4,201:length(A(1,:))));
% text()
% ylim([0 inf]);
%%
%choice 1 3.10-3.16
% for dex=17:length(acc)
% 
%     ppp=10*(dex+20)+1;
%     fitness=fitness+(A(4,ppp)-acc(dex))^2;
% end
% plot(37:0.1:44,A(4,371:1:441),'-.','color','b');
% hold on
% plot(37:1:44,acc(17:length(acc)),'*','color','r')

%%
% xlabel('Days after the first patient')
% ylabel('Population')
% % legend('Prediction value',['Actual value',sprintf('\n'),'(22/02/2020-16/03/2020)'])
% legend(['Actual value',sprintf('\n'),'(22/02/2020-16/03/2020)'])
% title('{\itPrediction for confirmed cases (early period)}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end%%%end%%%end%%%end%%%end%%%end%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%begin%%%%%%%%%begin%%%%%%%%begin%%%%%%%%begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%name:mypso.m
%Then,here comes the PSO code, pay attention that you have to definite the myFitCal.m first in order to ensure you can use the fitness %function in the following code.

clf
clc
clear 
close all
E=0.0001;
maxnum=50;%
narvs=2;%pay attention that you are supposed to change this parameter during different periods
particlesize=50;%
c1=2;%
c2=2;%
w=0.6;%
vmax=1;%
v=rand(particlesize,narvs);%
x=rand(particlesize,narvs);%
%%
for i=1:particlesize	
    f(i)=myFitCal(x(i,1),x(i,2));	%change size with different peiods as well
end
personalbest_x=x;
personalbest_faval=f;
[globalbest_faval,i]=min(personalbest_faval);
globalbest_x=personalbest_x(i,:); 
k=1;
while (k<=maxnum)	
    for i=1:particlesize			
        f(i)=myFitCal(x(i,1),x(i,2));		%change size with different peiods as well
        if f(i)<personalbest_faval(i)		
            personalbest_faval(i)=f(i);			
            personalbest_x(i,:)=x(i,:);		
        end
    end
    [globalbest_faval,i]=min(personalbest_faval);	
    globalbest_x=personalbest_x(i,:);	
    for i=1:particlesize		
        v(i,:)=w*v(i,:)+c1*rand*(personalbest_x(i,:)-x(i,:))+c2*rand*(globalbest_x-x(i,:));		
        for j=1:narvs			
            if v(i,j)>vmax				
                v(i,j)=vmax;
            elseif v(i,j)<-vmax				
                v(i,j)=-vmax;            
            end
        end
        m=x(i,:);
        x(i,:)=m+v(i);
        if(x(i,1)>=1||x(i,1)<=0||x(i,2)>=1||x(i,2)<=0) %change size with different peiods as well
            x(i,:)=m;
        end
    end
    ff(k)=globalbest_faval;
    if globalbest_faval<E        
        break    
    end
%     figure(1)%
%     for i= 1:particlesize%       
%         plot(x(i,1),'*')%       
%     end
    k=k+1;
end
xbest=globalbest_x;


% hold on;
% figure(1)
% set(gcf,'color','white');
% plot(1:length(ff),ff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end%%%end%%%end%%%end%%%end%%%end%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%begin%%%%%%%%%begin%%%%%%%%begin%%%%%%%%begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end%%%end%%%end%%%end%%%end%%%end%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%begin%%%%%%%%%begin%%%%%%%%begin%%%%%%%%begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%name:differentalpha
%this code was built in order to show the bound of alpha (Figure 7)
%pay attention that you should cancle the shield of plot code in myFitCal.m

for i= 0.2:0.01:0.28
myFitCal(i);
hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end%%%end%%%end%%%end%%%end%%%end%%%%%%%%%%%%%%%%%%%%%%%%%%
