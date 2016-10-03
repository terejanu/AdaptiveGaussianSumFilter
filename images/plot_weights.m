w2=[];
for i = 1 : 101
    w2 = [w2 F2.w{i}'];
end;

w3 = F4.w{1}';
for i = 2 : 101
    w3 = [w3 F4.w{i}];
end;

figure; hold on;
plot(time.tspan,w2(1,:),'b','LineWidth',2);
plot(time.tspan,w2(2,:),'b','LineWidth',2);
plot(time.tspan,w3(1,:),'r','LineWidth',2);
plot(time.tspan,w3(2,:),'r','LineWidth',2);
xlabel('time [sec]','FontSize',14);
ylabel('weights','FontSize',14);
hold off;