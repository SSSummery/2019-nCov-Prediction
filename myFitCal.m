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
