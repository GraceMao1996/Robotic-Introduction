%% ��D-H����������������ģ��
L1=Link([0 0 0 -pi/2 0]);
L2=Link([-pi/2 100 0 pi/2 0]);
L3=Link([0 100 0 0 0]);
L4=Link([pi/2 0 100 -pi/2 0]);
L5=Link([0 0 0 pi/2 0]);
L6=Link([0 0 0  0 0]);
Robot=SerialLink([L1,L2,L3,L4,L5,L6]);
Robot.name='Grace��s Robot--DH';
Robot.display();
Robot.plot([0 0 1 0 0 0]);
% % %% ��matlab���ż��������˶�ѧ�⣬���ù�������֤
% syms theta1 theta2 theta3 theta4 theta5 theta6;
% alpha=[Robot.links.alpha];
% d=[Robot.links.d];
% a=[Robot.links.a];
% theta=[theta1 theta2 theta3 theta4 theta5 theta6];
% T=eye(4);
% for i=1:6
%     A=[cos(theta(i)) -sin(theta(i))*cos(alpha(i)) sin(theta(i))*sin(alpha(i)) 0;
%        sin(theta(i)) cos(theta(i))*cos(alpha(i)) -cos(theta(i))*sin(alpha(i)) 0;
%        0 sin(alpha(i)) cos(alpha(i)) d(i);
%        0   0   0  1];
%    T=T*A;
% end
% T=simplify(T)%�������˶�ѧ��
% zyd=Robot.fkine(theta)%������������˶�ѧ��
%% �ؽڿռ�㵽��켣�滮����ζ���ʽ�켣�滮
theta0=[0 0 0 0 0 0];
theta1=[pi/6 pi/3 pi/2 -pi/2 -pi/3 -pi/6];t=6;
c=zeros(6,1);q=zeros(61,6);dq=zeros(61,6);ddq=zeros(61,6);
A=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0;
   1 t t*t t^3 t^4 t^5;
   0 1 2*t 3*t*t 4*t^3 5*t^4;
   0 0  2  6*t 12*t*t 20*t^3];

for k=1:6
    B=[theta0(k);0;0;theta1(k);0;0];
    c=A^(-1)*B;
    for i=1:61
        t=(i-1)*0.1;
        q(i,k)=c(1)+c(2)*t+c(3)*t^2+c(4)*t^3+c(5)*t^4+c(6)*t^5;
        dq(i,k)=c(2)+2*c(3)*t+3*c(4)*t^2+4*c(5)*t^3+5*c(6)*t^4;
        ddq(i,k)=2*c(3)+6*c(4)*t+12*c(5)*t^2+20*c(6)*t^3;
    end
end
Robot.plot(q);   
figure;
subplot(3,1,1);plot(q);title('��ζ���ʽλ��');
subplot(3,1,2);plot(dq);title('�ٶ�');
subplot(3,1,3);plot(ddq);title('���ٶ�');
% %% �ؽڿռ�㵽��켣�滮�������߹��ɵ������˶��켣
clear q;
theta0=[0 0 0 0 0 0]*pi/180;v=[5 5 6 10 10 7]*pi/180;
theta1=[20 20 -30 40 -40 -35]*pi/180;t=6;
N=600;q=zeros(N,6);dq=zeros(N,6);ddq=zeros(N,6);
for i=1:6
    if theta1(i)-theta0(i)<0
        v(i)=-v(i);
    end
    tb=fix((theta0(i)-theta1(i)+v(i)*t)/v(i)*100)/100;
    thea=theta0(i)+v(i)*tb/2;
    for j=1:N
        te=(j-1)*t/N;
        if te<=tb
            q(j,i)=theta0(i)+1/2*v(i)/tb*te^2;
            dq(j,i)=v(i)/tb*te;
            ddq(j,i)=v(i)/tb;
        else
            if te>tb && te<(t-tb)
               q(j,i)=q(tb*100,i)+v(i)*(te-tb);
               dq(j,i)=v(i);
               ddq(j,i)=0;
            else
                q(j,i)=theta1(i)-(v(i)/tb*(t-te)^2)/2;
                dq(j,i)=v(i)/tb*(t-te);
                ddq(j,i)=-v(i)/tb;
            end
        end
    end
end
d=q(1:10:end,:);
figure;
Robot.plot(d);   
figure;
subplot(3,1,1);plot(q);title('������λ��');
subplot(3,1,2);plot(dq);title('�ٶ�');
subplot(3,1,3);plot(ddq);title('���ٶ�');
%% �ڹ����ռ�ʵ��ֱ���˶�
%���ó�ʼ��Ϊ��50��100,100�����յ�Ϊ��100��150,100��
sp=[50 100 100];ep=[100 150 100];
N=60;q=zeros(N,6);
figure;
for i=1:N
    temp=(ep-sp)/N*i+sp;
    
    plot3(temp(1),temp(2),temp(3),'*r');
    hold on;
    q(i,:)=ikine(Robot,SE3(temp));
end
Robot.plot(q);
plot(q);
%% �ڿռ�ʵ��Բ���˶�
n=[0 0 1]; %������n
r=25; %Բ�İ뾶Ϊ1
c=[75 125 100]; %Բ�ĵ�����
theta=(0:2*pi/100:2*pi)'; %theta�Ǵ�0��2*pi
a=cross(n,[1 0 0]); %n��i��ˣ���ȡa����
q=zeros(100,6);
if ~any(a) %���aΪ����������n��j���
    a=cross(n,[0 1 0]);
end
b=cross(n,a); %��ȡb����
a=a/norm(a); %��λ��a����
b=b/norm(b); %��λ��b����

c1=c(1)*ones(size(theta,1),1);
c2=c(2)*ones(size(theta,1),1);
c3=c(3)*ones(size(theta,1),1);

x=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%Բ�ϸ����x����
y=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%Բ�ϸ����y����
z=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);%Բ�ϸ����z����
figure;
plot3(x,y,z);
hold on;
for i=1:length(x)
    q(i,:)=ikine(Robot,SE3([x(i) y(i) z(i)]));
end

Robot.plot(q);
%% �����˶���Χ
p1=[50 100 100];p2=[150 300 270];p3=[300 200 400];
p=[p1;p2;p3];
fanwei=[-160 150;
    -150 75;
    -170 155;
    -180 180;
    -120 130;
    -180 180]*pi/180;
for i=1:3
    point=Robot.ikine(SE3(p(i,:)));
    p(i,:)
    if isempty(point)
        fprintf('can not reach!\n');
    else
        for i=1:6
            if point(i)<fanwei(i,1)||point(i)>fanwei(i,2)
                 fprintf('can not reach!\n');
             return
            end
        end
     fprintf('the point can reach!\n');
    end
end
%% ����ĩ��ִ������λ����
p1=[50 100 100];
clear q;
q=ikine(Robot,SE3(p1));
J=Robot.jacob0(q);
B=[0.5 0.5 0.5 0.5 0.5 0.5]'*pi/180;
D=J*B;
fprintf('��е���ڵ㣨%d,%d,%d����x,y,z��Ķ�λ���ֱ�Ϊ��%f,%f,%f\n',...
        p1(1),p1(2),p1(3),D(1),D(2),D(3));