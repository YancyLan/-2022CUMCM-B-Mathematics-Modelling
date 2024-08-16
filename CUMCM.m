// 附录
// 附录1
// 介绍： 支撑材料的文件列表
// 1.	Q1A1.m，第一题第（1）小问源程序
// 2.	Q1A2.m，第一题第（2）小问源程序
// 3.	Q1A31.m，第一题第（3）小问第一次迭代源程序
// 4.	Q1A32.m，第一题第（3）小问第二次迭代源程序
// 5.	Q1A33.m，第一题第（3）小问第三次迭代源程序
// 6.	Q1A3_calerror.m，第一题第（3）小问计算误差源程序
// 7.	plot_1.m，第一题第（3）小问作图代码
// 8.	plot2.m，第二题与稳健性分析作图代码
// 9.	第三问迭代过程与结果.xlsx
// 10.	iter.m，第一题第（3）小问有偏差位置移向理想位置计算函数

% 附录2
% 介绍： 第一题第（1）小问信息源回溯算法
t=0:2*pi/9:2*pi;
x=[0 100*cos(t)];y=[0 100*sin(t)];
k=5;%有偏差无人机的编号（FY00、FY03、FY01除外）,假设为5,已知编号为FY00、FY03、FY01
B=[x(k+1),y(k+1)];
l2=[];
l1=[[0,1];[0,2];[0,3];[0,4];[0,6];[0,7];[0,8];[0,9];[1,2];[1,3];[1,4];[1,6];[1,7];[1,8];[1,9];[2,3];[2,4];[2,6];[2,7];[2,8];[2,9];[3,4];[3,6];[3,7];[3,8];[3,9];[4,6];[4,7];[4,8];[4,9];[6,7];[6,8];[6,9];[7,8];[7,9];[8,9]];
for i=1:36
    M=l1(i,:);
    m=M(1)+1;n=M(2)+1;
    C=[x(m),y(m)];
    A=[x(n),y(n)];
    a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
    b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
    c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
    cosB=(a^2+c^2-b^2)/(2*a*c);
    l2=[l2 cosB];
end
alpha1=pi/18+0.1;%输入有偏差无人机接收到的FY00与未知编号无人机的夹角信息
alpha2=2*pi/9+0.001;%输入有偏差无人机接收到的FY00与未知编号无人机的夹角信息
alpha3=5*pi/18+0.07;%输入两架有偏差无人机接收到的FY00与未知编号无人机的夹角信息
C1=abs(l2-cos(alpha1))
C1_min=min(C1)
[row, col1]=find(C1==C1_min);
m1=l1(col1,:)%输出alpha1信号源两架无人机的编号（无顺序）
C2=abs(l2-cos(alpha2));
C2_min=min(C2);
C2_min=6.431705247449404e-04;%由min(C2)得出结果，忽略MATLAB的数字显示问题
[row, col2]=find(C2==C2_min);
m2=l1(col2,:)%输出alpha2信号源两架无人机的编号（无顺序）
C3=abs(l2-cos(alpha3));
C3_min=min(C3);
[row, col3]=find(C3==C3_min);
m3=l1(col3,:)%输出alpha3信号源两架无人机的编号（无顺序）


% 附录3
% 介绍： 第一题第（2）小问改进信息源回溯算法
t=0:2*pi/9:2*pi;
x=[0 100*cos(t)];y=[0 100*sin(t)];
k=3;%有偏差无人机的编号（FY00、FY01、FY0k除外）,在这里假设情况为3,已知编号FY00、FY01
B=[x(k+1),y(k+1)];
l3=[];l4=[];
l1=[[0,2];[0,4];[0,5];[0,6];[0,7];[0,8];[0,9]];
l2=[[1,2];[1,4];[1,5];[1,6];[1,7];[1,8];[1,9]]; %每架无人机有自己独一无二的l1与l2
for i=1:7
    M=l1(i,:);
    m=M(1)+1;n=M(2)+1;
    C=[x(m),y(m)];
    A=[x(n),y(n)];
    a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
    b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
    c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
    cosB=(a^2+c^2-b^2)/(2*a*c);
    angle=acos(cosB);
    reangle=round(angle*180/pi);
    l3=[l3 reangle];
end
for i=1:7
    M=l2(i,:);
    m=M(1)+1;n=M(2)+1;
    C=[x(m),y(m)];
    A=[x(n),y(n)];
    a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
    b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
    c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
    cosB=(a^2+c^2-b^2)/(2*a*c);
    angle=acos(cosB);
    reangle=round(angle*180/pi);
    l4=[l4 reangle];
end
alphax=7*pi/18+0.07;%输入有偏差无人机接收到的FY00与未知编号无人机的夹角信息
alphay=pi/9+0.005;%输入有偏差无人机接收到的FY01与未知编号无人机的夹角信息
alphax=round(alphax*180/pi/10)*10;alphay=round(alphay*180/pi/10)*10;
[row, X1]=find(l3==alphax);
[row, X2]=find(l3==alphay);
if isempty(X1)==0
    L4=l4(:,X1);
    if isempty(find(L4==alphay))==0
        c1=X1(find(L4==alphay));
        m1=l1(c1,:)%输出信号源FY00与未知编号无人机的编号（无顺序）
        m2=l2(c1,:)%输出信号源FY01与未知编号无人机的编号（无顺序）
    end
else
    L4=l4(:,X2);
    if isempty(find(L4==alphax))==0
        c2=X2(find(L4==alphax));
        m1=l1(c1,:)%输出信号源FY00与未知编号无人机的编号（无顺序）
        m2=l2(c1,:)%输出信号源FY01与未知编号无人机的编号（无顺序）
    end
end


% 附录4
% 介绍： 第一题第（3）小问作图代码
%理想位置
b=[0 0
    100	0
76.6044	64.2788
17.3648	98.4808
-50	86.6025
-93.9693	34.202
-93.9693	-34.202
-50	-86.6025
17.3648	-98.4808
76.6044	-64.2788];
scatter(b(:,1)',b(:,2)','ro')
hold on

aplha=0:pi/40:2*pi;
r=100;
x=r*cos(aplha);
y=r*sin(aplha);
plot(x,y,'k-');
axis equal
hold on

%初始位置
% c=[0	0
% 100	0
%     76.9623	63.1241
% 19.0442	110.369
% -52.1027	91.1609
% -92.0077	33.7429
% -105.2723	-38.2328
% -52.3889	-90.9967
% 17.3038	-96.4602
% 86.1478	-71.5721];
% scatter(c(:,1)',c(:,2)','b*')
% hold on

%第一次迭代后的不发信号无人机的位置
% c=[13.2557	96.5947
% -50.126	83.2211
% -92.8944	33.535
% -93.4266	-33.9651
% -51.3673	-88.7565
% 20.8285	-101.0852
% 80.9321	-61.3664];
% scatter(c(:,1)',c(:,2)','b*')
% hold on

%第一次迭代后的发信号无人机的位置
% c=[0	0
% 100	0
% 76.9623	63.1241];
% scatter(c(:,1)',c(:,2)','bv')
% hold on

%第二次迭代后的不发信号无人机的位置
% d1=[76.5466	64.313
% 17.1455	98.3329
% -49.9856	86.1775
% -93.4933	33.9298
% -50.0054	-85.9862
% 15.799	-92.9819
% 76.3405	-64.4687];
% scatter(d1(:,1)',d1(:,2)','b*')
% hold on
%第二次迭代后的发信号无人机的位置
% d2=[0 0
%     100 0
%     -93.426 -33.9651];
% scatter(d2(:,1)',d2(:,2)','bv')
% hold on

%第三次迭代后的发信号无人机的位置
d2=[0 0
    100 0
    17.1455 98.3329];
scatter(d2(:,1)',d2(:,2)','bv')
hold on

%第三次迭代后的不发信号无人机的位置
d1=[76.8107	64.1606
-49.9995	86.445
-93.9638	34.1988
-94.066	-34.2576
-50.0019	-86.853
17.6417	-98.6529
76.9563	-64.0831];
scatter(d1(:,1)',d1(:,2)','b*')
hold on

%无人机平均误差距离折线图
x=0:1:3;%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
 a=[5.128424326,2.307245791,0.855633382,0.175687179]; %a数据y值
 plot(x,a,'-*k'); %线性，颜色，标记
axis([0,3,0,6])  %确定x轴与y轴框图大小
set(gca,'XTick',[0:1:3]) 
set(gca,'YTick',[0:6:6]) 
legend('每架无人机平均误差距离'); 
xlabel('迭代次数')  %x轴坐标描述
ylabel('平均误差距离（m）') %y轴坐标描述

% 附录5
% 介绍： 第一题第（3）小问第一次迭代代码
clear;clc
t=0:2*pi/9:2*pi;
x=[0 100*cos(t)];y=[0 100*sin(t)];
C=[112*cos(80.21*pi/180),112*sin(80.21*pi/180)];
C=[105*cos(119.75*pi/180),105*sin(119.75*pi/180)];
C=[98*cos(159.86*pi/180),98*sin(159.86*pi/180)];
%C=[112*cos(199.96*pi/180),112*sin(199.96*pi/180)];
%C=[105*cos(240.7*pi/180),105*sin(240.7*pi/180)];
C=[98*cos(280.17*pi/180),98*sin(280.17*pi/180)];
%C=[112*cos(320.28*pi/180),112*sin(320.28*pi/180)];

A=[x(1),y(1)];B=[x(2),y(2)];D=[98*cos(40.1*pi/180),98*sin(40.1*pi/180)];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=70*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=68;c3=100;
syms x;
f=c1*sin(a1+x)/sin(a1)-c2*sin(a2+d-x)/sin(a2);
p0=0;
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x


%点6特殊情况
t=0:2*pi/9:2*pi;
x=[0 100*cos(t)];y=[0 100*sin(t)];
C=[112*cos(199.96*pi/180),112*sin(199.96*pi/180)];

A=[x(1),y(1)];B=[x(2),y(2)];D=[98*cos(40.1*pi/180),98*sin(40.1*pi/180)];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=40*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=68;c3=100;
syms x;
f=sin(pi-a2-d+x)/sin(pi-a1-x)-sin(a2)/sin(a1);
p0=0;
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x

% 附录6
介绍： 第一题第（3）小问第二次迭代代码
%点2、3、4、5特殊情况
C=[76.9623	63.1241];
%C=[13.2557	96.5947];
%C=[-50.126	83.2211];
%C=[-92.8944 33.535];

A=[0 0];B=[100 0];D=[-93.4266,-33.9651];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=200*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=100;c3=196.9616;
syms x;
f=sin(a2)/sin(a1)-sin(pi-d+x-a2)/sin(pi-x-a1);
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x


%7、8、9特殊情况
C=[-51.3673 -88.7565];
C=[20.8285	-101.0852];
C=[80.9321	-61.3664];

A=[0 0];B=[100 0];D=[-93.4266,-33.9651];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=160*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=100;c3=196.9616
syms x;
f=sin(a2)/sin(a1)-sin(pi-d+x-a2)/sin(pi-x-a1);
p0=0;
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x

function [p,k,Y]=NTM(f,p0,tol,maxK)
%p0表示迭代初始值
%f表示要求解的方程
%maxK表示规定的最大迭代次数
%tolr表示允许误差
%k表示最终迭代的次数
%p表示最终迭代的值
    syms x;
    P(1)=p0;
    k=2;
    df=diff(f);     %利用diff()函数计算f(x)的导数
   
    P(k)=P(k-1)-subs(f,x,P(k-1))/subs(df,x,P(k-1)); %第二次迭代的结果
    while k<=maxK
        err=abs(P(k)-P(k-1));    %err表示相邻的迭代值的差值
        if(err<tol)
            fprintf('迭代%d次即可满足允许误差值退出\n',k-1);
            break;
        end
        k=k+1;
        P(k)=P(k-1)-subs(f,x,P(k-1))/subs(df,x,P(k-1));  %迭代
    end      %共迭代了k-1次
    if(k-1==maxK) 
        disp("超过最大迭代次数！");
    end
    p=P(k); 
    k=k-1;
    Y=P;
end

% 附录7
% 介绍： 第一题第（3）小问第三次迭代代码
%点6、7特殊情况
C=[-93.4266 -33.9651];
C=[-50.0054 -85.9862];
A=[0 0];B=[100 0];D=[17.1455 98.3329];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=-80*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=100;c3=128.5575;
syms x;
f=sin(a2)/sin(a1)-sin(pi-x-a2)/sin(pi-d+x-a1);
%f=sin(a1)/sin(a2)-sin(pi-x-a1)/sin(pi-d+x-a2);具体情况调用不同函数
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x
%点2特殊情况
C=[76.5466 64.313];
A=[0 0];B=[100 0];D=[17.1455 98.3329];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=80*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=100;c3=128.5575;
syms x;
f=sin(a2)/sin(a1)-sin(pi-x-a2)/sin(pi-d+x-a1);
%f=sin(a1)/sin(a2)-sin(pi-x-a1)/sin(pi-d+x-a2);具体情况调用不同函数
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x
%点4、5特殊情况
C=[-49.9856 86.1775];
C=[-93.4933 33.9298];
A=[0 0];B=[100 0];D=[17.1455 98.3329];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=50*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=100;c3=128.5575;
syms x;
f=100*sin(a2)/(128.5575*sin(a1))-sin(pi-d+x-a2)/sin(pi-x-a1);
%f=sin(a1)/sin(a2)-sin(pi-x-a1)/sin(pi-d+x-a2);具体情况调用不同函数
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x
%点8、9特殊情况
C=[15.799 -92.9819];
C=[76.3405 -64.4687];
A=[0 0];B=[100 0];D=[17.1455 98.3329];
a=sqrt((D(1)-B(1))^2+(D(2)-B(2))^2);
b=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
cosC=(a^2+b^2-c^2)/(2*a*b);
d=50*pi/180
a=sqrt((C(1)-B(1))^2+(C(2)-B(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
c1=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa1=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((A(1)-C(1))^2+(A(2)-C(2))^2);
c=sqrt((A(1)-D(1))^2+(A(2)-D(2))^2);
c2=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa2=acos(cosC)
a=sqrt((C(1)-D(1))^2+(C(2)-D(2))^2);
b=sqrt((B(1)-C(1))^2+(B(2)-C(2))^2);
c=sqrt((B(1)-D(1))^2+(B(2)-D(2))^2)
c3=c;
cosC=(a^2+b^2-c^2)/(2*a*b);
aa3=acos(cosC);
m=[aa1,aa2,aa3];
cc=[c1,c2,c3];
a3=max(m);
c3=cc(find(m==a3));
xxx=find(m==a3);
[ad i]=setdiff(m,a3);
cc(xxx)=[];
m=m(sort(i))
a1=m(1);a2=m(2);
c1=100;c2=100;c3=128.5575;
syms x;
f=sin(a2)/sin(a1)-sin(pi-d+x-a2)/sin(pi-x-a1);
tol=0.0000001;
maxK=1000;
[p,k,Y]=NTM(f,p0,tol,maxK);
fprintf("迭代值如下：");
disp(Y);
x=Y(end)
y=d-x
function [p,k,Y]=NTM(f,p0,tol,maxK)
%p0表示迭代初始值
%f表示要求解的方程
%maxK表示规定的最大迭代次数
%tolr表示允许误差
%k表示最终迭代的次数
%p表示最终迭代的值
    syms x;
    P(1)=p0;
    k=2;
    df=diff(f);     %利用diff()函数计算f(x)的导数
   
    P(k)=P(k-1)-subs(f,x,P(k-1))/subs(df,x,P(k-1)); %第二次迭代的结果
    while k<=maxK
        err=abs(P(k)-P(k-1));    %err表示相邻的迭代值的差值
        if(err<tol)
            fprintf('迭代%d次即可满足允许误差值退出\n',k-1);
            break;
        end
        k=k+1;
        P(k)=P(k-1)-subs(f,x,P(k-1))/subs(df,x,P(k-1));  %迭代
    end      %共迭代了k-1次
    if(k-1==maxK) 
        disp("超过最大迭代次数！");
    end
    p=P(k); 
    k=k-1;
    Y=P;
end

% 附录8
% 介绍： 第一题第（3）小问计算误差代码（以第一次迭代为例）
a1=[76.6044 17.3648 -50 -93.9693    -93.9693    -50 17.3648 76.6044
64.2788 98.4808 86.6025 34.202  -34.202 -86.6025    -98.4808    -64.2788];
a2=[76.9623 13.2557 -50.126 -92.8944    -93.4266    -51.3673    20.8285 80.9321
63.1241 96.5947 83.2211 33.535  -33.9651    -88.7565    -101.0852   -61.3664];
m=a1-a2
l=[]
for i=1:8
    l=[l sqrt(m(1,i)^2+m(2,i)^2)]
end


% 附录9
% 介绍： 第一题第（3）问有偏差位置移向理想位置计算函数
%(x0,y0)是理想位置的坐标
%（x1，y1）是假的实际位置
%（x2，y2）是实际位置
%（x3,y3）是迭代出的下一次的位置
function [x3,y3]=iter(x0,y0,x1,y1,x2,y2)
x3=x2+x0-x1;
y3=y2+y0-y1;
end


% 附录10
% 介绍： 第二题与稳健性分析作图代码
%第二题
%%
a=[0,0
    25*sqrt(3),25
    50*sqrt(3),50
    75*sqrt(3),75
    100*sqrt(3),100
    0,50
    25*sqrt(3),75
    50*sqrt(3),100
    75*sqrt(3),125
    0,100
    25*sqrt(3),125
    50*sqrt(3),150
    0,150
    25*sqrt(3),175
    0,200
    ]
scatter(a(:,1),a(:,2),'b*')
hold on
for i=1:max(size(a(:,1)))
    c = num2str(i-1);
    text(a(i,1),a(i,2),c);
end
a=[0,0
    25*sqrt(3),25
    50*sqrt(3),50
    75*sqrt(3),75
    100*sqrt(3),100
    0,50
    75*sqrt(3),125
    0,100
    50*sqrt(3),150
    0,150
    25*sqrt(3),175
    0,200];
b=[25*sqrt(3),75
   50*sqrt(3),100
   25*sqrt(3),125];
scatter(a(:,1),a(:,2),'b*')
hold on
scatter(b(:,1),b(:,2),'r*')
hold on
x=25*sqrt(3);y=75;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
hold on
x=50*sqrt(3);y=100;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
hold on
x=25*sqrt(3);y=125;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
%%
a=[0,0
    25*sqrt(3),25
    50*sqrt(3),50
    75*sqrt(3),75
    100*sqrt(3),100
    0,50
    75*sqrt(3),125
    0,100
    50*sqrt(3),150
    0,150
    25*sqrt(3),175
    0,200];
b=[25*sqrt(3),75
   50*sqrt(3),100
   25*sqrt(3),125];
scatter(a(:,1),a(:,2),'b*')
hold on
scatter(b(:,1),b(:,2),'r*')
hold on
x=25*sqrt(3);y=75;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
hold on
x=50*sqrt(3);y=100;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
hold on
x=25*sqrt(3);y=125;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
%%
a=[0,0
   50*sqrt(3),50
   50*sqrt(3),150
   0,200
   0,100
   100*sqrt(3),100];
b=[0,150
   25*sqrt(3),175
   25*sqrt(3),125
   0,50
   25*sqrt(3),25
   25*sqrt(3),75
   75*sqrt(3),75
   75*sqrt(3),125
   50*sqrt(3),100
   
    ];
scatter(a(:,1),a(:,2),'b*')
hold on
scatter(b(:,1),b(:,2),'r*')
hold on
x=0;y=150;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
hold on
x=0;y=50;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
hold on
x=75*sqrt(3);y=75;r=50;
rectangle('Position', [x-r,y-r,2*r,2*r], 'Curvature', [1 1],'EdgeColor', 'r');
%%
a=[0,0
    25*sqrt(3),25
    50*sqrt(3),50
    75*sqrt(3),75
    100*sqrt(3),100
    0,50
    75*sqrt(3),125
    0,100
    50*sqrt(3),150
    0,150
    25*sqrt(3),175
    0,200
    -25*sqrt(3),-25
    -25*sqrt(3),25
    -25*sqrt(3),75
    -25*sqrt(3),125
    -25*sqrt(3),175
    -25*sqrt(3),225
    ]
b=[50*sqrt(3),100
   25*sqrt(3),75
   25*sqrt(3),125  
    ];
scatter(a(:,1),a(:,2),'b*')
hold on
scatter(b(:,1),b(:,2),'r*')
hold on

%稳健性分析可视化
box=[0.573 1.1459 1.7189 2.2918 2.8648 3.4377 4.0107 4.5837 5.1566 5.7296
1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	2	2];
scatter(box(1,:),box(2,:),'k*')
hold on
scatter(box(1,:),box(3,:),'kv')
hold on
scatter(box(1,:),box(4,:),'ko')
hold on

