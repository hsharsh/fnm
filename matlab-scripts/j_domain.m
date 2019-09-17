clc
clear
close all

cd /home/hsharsh/fnm

x = load('nodes.inp');
x = x(:,2:3);
conn = load('elements.inp');
conn = conn(:,2:5);

hold on;
plot(x(:,1),x(:,2),'b.');
axis equal
grid on

inner = [6 7 10 11 ];
outer = [1 2 3 4 5 8 9 12 13 14 15 16 ];
domain = [387 388 389 407 408 409 427 428 429];
domain = domain -1;
% plot(x(inner,1),x(inner,2),'rx');
% plot(x(outer,1),x(outer,2),'gx');

for i = 1:4
    plot(x(conn(domain,i),1),x(conn(domain,i),2), 'bo');
end
