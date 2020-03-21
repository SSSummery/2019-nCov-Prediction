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
