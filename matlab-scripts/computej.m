clc
clear
close all

cd /home/hsharsh/fnm

x = load('nodes.inp');
x = x(:,2:3);
conn = load('elements.inp');
conn = conn(:,2:5);
un = load('un1.inp');


plot(x(:,1),x(:,2),'b.');
axis equal
grid on

elem = [201, 202, 221, 222, 241, 242, 261, 262, 411, 412, 413, 414, 431, 432, 433, 434];
hold on;
% for i = 1:length(elem)
%     xmid = 0; ymid = 0;
%     for j = 1:4
%         xmid = xmid+x(conn(elem(i),j),1);
%         ymid = ymid+x(conn(elem(i),j),2);
%         plot(x(conn(elem(i),j),1),x(conn(elem(i),j),2),'x');
%         pause
%     end
%     xmid = xmid/4;
%     ymid = ymid/4;
%     plot(xmid,ymid,'r.');
% end

xgp = sqrt(3/5)*[-1 0 1];
wgp = [5 8 5]/9;
ngp = length(xgp);
E = 1;
nu = 0;

qn = [
    0 0 1 1
    0 0 0 1
    1 1 1 1
    1 0 0 1
    1 1 1 1
    1 0 0 1
    1 1 0 1
    1 0 0 0
    0 0 1 1
    1 1 1 1
    1 1 1 1
    1 1 1 0
    0 0 1 0
    0 1 1 0
    0 1 1 0
    0 1 0 0
];

j_int = 0;

for lmn = 1:length(elem)

    nodes = conn(elem(lmn),:);
    xvec = x(nodes,1);
    yvec = x(nodes,2);
    xmid = sum(xvec)/4;
    ymid = sum(yvec)/4;

    dof =  [
        nodes(1)*2-1 nodes(1)*2 ...
        nodes(2)*2-1 nodes(2)*2 ...
        nodes(3)*2-1 nodes(3)*2 ...
        nodes(4)*2-1 nodes(4)*2 ...
    ];
    u = un(dof)';
    q = qn(lmn,:)';
    plot(xmid,ymid,'r.');

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

% nset = [];
% n = 0;
% for lmn = 1:length(elem)
%     nodes = conn(elem(lmn),:);
%     for i = 1:4
%         if(~ismember(nodes(i),nset))
%             n = n + 1
%             nset(n) = nodes(i);
%         end
%     end
% end
