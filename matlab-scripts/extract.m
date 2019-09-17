clc
clear
close all

cd /home/hsharsh/fnm

xfull = load('nodes.inp');
xfull = xfull(:,2:3);
connfull = load('elements.inp');
connfull = connfull(:,2:5);
un1 = load('un1.inp');

% plot(x(:,1),x(:,2),'b.');
% axis equal
% grid on

index = [
0		3
0		64
0		79
0		80
0		81
0		329
0		330
0		348
0		349
0		367
0		368
0		386
0		387
0		405
0		406
0		510
0		511
0		512
0		513
0		514
0		529
0		530
0		531
0		532
0		533
];

si = [25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 3, 4, 5, 1, 2, 14, 12, 10, 8, 6, 15, 13, 11, 9, 7];
nodes = index(si,2);

x = xfull(nodes+1,:);
un = [un1(nodes+1); un1(nodes*2+1)];

conn = [
1 2 7 6
2 3 8 7
3 4 9 8
4 5 10 9
6 7 12 11
7 8 13 12
8 9 14 13
9 10 15 14
11 12 17 16
12 13 18 17
13 14 19 18
14 15 20 19
16 17 22 21
17 18 23 22
18 19 25 23
19 20 25 24
];

qn = [0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 0 0 0]';

qn = [0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 1 1 1 0 0 0 0 0 0]';

xgp = sqrt(3/5)*[-1 0 1];
wgp = [5 8 5]/9;
ngp = length(xgp);
E = 1;
nu = 0;


j_int = 0;
conn = conn([6 7 10 11]',:);

for lmn = 1:length(conn)

    nodes = conn(lmn,:);
    xvec = x(nodes,1);
    yvec = x(nodes,2);
    dof =  [
        nodes(1)*2-1 nodes(1)*2 ...
        nodes(2)*2-1 nodes(2)*2 ...
        nodes(3)*2-1 nodes(3)*2 ...
        nodes(4)*2-1 nodes(4)*2 ...
    ];
    u = un(dof);
    q = qn(nodes);
    for i = 1:ngp
        for j=1:ngp
            r = xgp(i); s = xgp(j);
            B0 = [ -(1-s)  (1-s) (1+s) -(1+s);
                 -(1-r)  -(1+r)  (1+r)  (1-r)]/4;
            B1 = [1 0 0 0; 0 0 0 1; 0 0.5 0.5 0];

            jac(:,1) = (B0*xvec);
            jac(:,2) = (B0*yvec);

            B2(1:2,1:2) = inv(jac);
            B2(3:4,3:4) = inv(jac);

            B3 = zeros(4,8);
            B3(1:2, 1:2:end) = B0;
            B3(3:4, 2:2:end) = B0;

            B = B1*B2*B3;
            Bu = B2*B3;
            Bq = jac*B0;

            D = E/((1+nu)*(1-2*nu))*[1-nu 0 0; 0 1-nu 0; 0 0 1-2*nu];

            du = Bu*u;
            dq = Bq*q;
            strain = B*u;
            stress = D*strain;

            w = 0.5*(dot(stress,strain));
            j_int = j_int + (stress(1)*du(1)*dq(1) + stress(3)*du(3)*dq(1) + stress(3)*du(1)*dq(2) + stress(2)*du(3)*dq(2) - w*dq(1) )*wgp(i)*wgp(j);
        end
    end
end

disp(['J-integral: ',num2str(j_int)]);
% si = [
%     25 24 19 20
%     24 23 18 19
%     23 22 17 18
%     22 21 16 17
%     20 19 4 3
%     19 18 5 4
%     18 17 1 5
%     17 16 2 1
%     3 4 12 14
%     4 5 10 12
%     5 1 8 10
%     1 2 6 8
%     14 12 13 15
%     12 10 11 13
%     10 8 9 11
%     8 6 7 9
% ];
% for i = 1:length(si)
%     for j = 1:4
%         elem(i,j) = index(si(i,j),2);
%     end
% end


% for i = 1:length(elem)
%     for j = 1:4
%         plot([x(elem(i,mod((j-1),4)+1)+1,1),x(elem(i,mod((j),4)+1)+1,1)],[x(elem(i,mod((j-1),4)+1)+1,2),x(elem(i,mod((j),4)+1)+1,2)]);
%         pause(0.5)
%     end
% end


% plot(x(index(:,2)+1,1),str(index(:,2)+1),'r')
% hold on
% grid on
% fplot(@(r) (str(index(2,2)+1)*sqrt(0.45-0.35))/sqrt(r-0.35),[0.35 1],'k')
% plot(abaqus(:,1)+0.4,abaqus(:,2),'b')
% legend('Code','1/sqrt(r)','Abaqus')
% xlabel('x-distance (absolute)')
% ylabel('sig_y_y')
