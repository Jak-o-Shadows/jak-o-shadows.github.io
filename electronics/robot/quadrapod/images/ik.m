
%1-09-2017


close all
clear
clc



%%



%%


circle = @(centre, radius) plot(centre(1) + radius*cos(linspace(0, 2*pi, 50)), centre(2) +radius*sin(linspace(0, 2*pi, 50)), '--');
plotLink = @(p1, p2, col) plot([p1(1), p2(1)], [p1(2), p2(2)], col);

O = [0;0];
A = [25;-20];

OD = 200;
OC = 60;
BC = 45;
AB = 50;
DH = 140;

c = @(x) cosd(x);
s = @(x) sind(x);


desired = [0; -(OD + DH)*0.8];

forward = @(th3, th4) [-c(th3)*OD + DH*c(180-th3-th4); -s(th3)*OD - DH*s(180-th3-th4)];

yPosList = -(75:5:(DH+OD)*0.95);
servoAngles = NaN(length(yPosList), 2);
endPosList = NaN(length(yPosList), 2);

for i=1:length(yPosList);
    desired = [0; yPosList(i)];

    func = @(th3, th4) sqrt(sum((desired - [-c(th3)*OD + DH*c(180-th3-th4); -s(th3)*OD - DH*s(180-th3-th4)]).^2));
    funcSolve = @(x) func(x(1), x(2));

    options = optimset('MaxFunEvals', 10000);

    Aeqn = [];
    beqn = [];
    Aeqeqn = [];
    Beqeqn = [];
    test = fmincon(funcSolve, [0, 90], Aeqn, beqn, Aeqeqn, Beqeqn, [0, 0], [180, 180])



    %hence find the point B
    C = @(th3) O + OC*[-c(th3); -s(th3)];
    C = C(test(1));
    %intersection of two circles to get point B
    d = sqrt(sum((A-C).^2));
    a = (BC^2 - AB^2 + d^2)/(2*d);
    h = sqrt(BC^2 - a^2);
    P2 = C + a*(A-C)/d;

    B_1 = [P2(1) + h*(A(2) - C(2))/d; 
           P2(2) - h*(A(1) - C(1))/d];
    B_2 = [P2(1) - h*(A(2) - C(2))/d; 
           P2(2) + h*(A(1) - C(1))/d];



    % figure(1);
    % hold on;
    % plot([O(1), C(1)], [O(2), C(2)], 'k-');
    % text(O(1), O(2), 'O');
    % text(A(1), A(2), 'A');
    % text(C(1), C(2), 'C');
    % text(B_1(1), B_1(2), 'B_1');
    % text(B_2(2), B_2(2), 'B_2');
    % 
    % plot([C(1), B_1(1)], [C(2), B_1(2)], 'r-');
    % plot([A(1), B_1(1)], [A(2), B_1(2)], 'r-');
    % 
    % plot([C(1), B_2(1)], [C(2), B_2(2)], 'g-');
    % plot([A(1), B_2(1)], [A(2), B_2(2)], 'g-');
    % 
    % circle(C, BC);
    % circle(A, AB);
    % 
    % text(P2(1), P2(2), 'P_2')

    %work out angle theta_1
    theta_1_limits = [-90, 90];
    B = B_1;
    theta_1 = atan2d(-(B(1)-A(1)), -(B(2)-A(2)));
    if ((theta_1 <= theta_1_limits(2)) && (theta_1 >= theta_1_limits(1)))
        %accept
        fprintf('Theta_1 is %f deg\n', theta_1);
        fprintf('\t Point B_1\n')
    else
        B = B_2;
        theta_1 = atan2d(-(B(1)-A(1)), -(B(2)-A(2)));
        fprintf('Theta_1 is %f deg\n', theta_1);
        fprintf('\t Point B_2\n')
    end


    %do it for leg 2!
    J2x = 40;
    J2y = 30;
    L1 = 140;
    L2 = 60;
    L3 = 45;
    L4 = 50;

    %use local coordinates
    D = [0;0];
    E = [-J2x; J2y];

    %hence find the point B
    G = @(th4) L4*[-c(th4); -s(th4)];
    G = G(test(2));
    %intersection of two circles to get point B
    d = sqrt(sum((E-G).^2));
    a = (L3^2 - L2^2 + d^2)/(2*d);
    h = sqrt(L3^2 - a^2);
    P2 = G + a*(E-G)/d;

    F_1 = [P2(1) + h*(E(2) - G(2))/d; 
           P2(2) - h*(E(1) - G(1))/d];
    F_2 = [P2(1) - h*(E(2) - G(2))/d; 
           P2(2) + h*(E(1) - G(1))/d];



    % figure(2);
    % hold on;
    % plot([D(1), G(1)], [D(2), G(2)], 'k-');
    % text(D(1), D(2), 'D');
    % text(E(1), E(2), 'E');
    % text(G(1), G(2), 'G');
    % text(F_1(1), F_1(2), 'F_1');
    % text(F_2(2), F_2(2), 'F_2');
    % 
    % plot([G(1), F_1(1)], [G(2), F_1(2)], 'r-');
    % plot([E(1), F_1(1)], [E(2), F_1(2)], 'r-');
    % 
    % plot([G(1), F_2(1)], [G(2), F_2(2)], 'g-');
    % plot([E(1), F_2(1)], [E(2), F_2(2)], 'g-');
    % 
    % 
    % circle(G, L3);
    % circle(E, L2);
    % 
    % text(P2(1), P2(2), 'P_2')

    %work out angle theta_1
    theta_2_limits = [-270, -90];
    F = F_1;
    theta_2 = atan2d(-(F(1)-E(1)), -(F(2)-E(2)));
    if ((theta_2 <= theta_2_limits(2)) && (theta_2 >= theta_2_limits(1)))
        %accept
        fprintf('Theta_2 is %f deg\n', theta_2);
        fprintf('\t Point F_1\n')
    else
        F = F_2;
        theta_2 = atan2d(-(F(1)-E(1)), -(F(2)-E(2)));
        fprintf('Theta_2 is %f deg\n', theta_2);
        fprintf('\t Point F_2\n')
    end

    servoAngles(i, :) = [theta_1, theta_2];
    

    %hence plot on same graph
    D = (C-O)*(OD/OC) + O; %easier rescaling than doing trig


    figure(3);
    cla;
    axis equal;
    hold on;
    text(O(1), O(2), 'O');
    text(A(1), A(2), 'A');
    text(C(1), C(2), 'C');
    text(B(1), B(2), 'B');
    plotLink(O, C, 'k-');
    plotLink(O, D, 'k-');
    plotLink(C, B, 'k--');
    plotLink(A, B, 'k-.');


    text(D(1), D(2), 'D');

    G = D + [c(test(1)), -s(test(1)); s(test(1)), c(test(1))]*G;
    text(G(1), G(2), 'G');
    plotLink(D, G, 'b-');

    E = D + [c(test(1)), -s(test(1)); s(test(1)), c(test(1))]*E;
    text(E(1), E(2), 'E');

    F = D + [c(test(1)), -s(test(1)); s(test(1)), c(test(1))]*F;
    text(F(1), F(2), 'F');
    plotLink(E, F, 'b-.');
    plotLink(F, G, 'b--');

    H = (G-D)*(L1/L4) + D;
    text(H(1), H(2), 'H');
    plotLink(D, H, 'b-');
    
    endPosList(i, :) = [H(1), H(2)];
    plot(endPosList(:, 1), endPosList(:, 2), 'go');
    xlabel('x (mm)');
    ylabel('y (mm)');
    
    print(gcf(), sprintf('%d.png', i), '-dpng', '-r300');


end

figure;
hold on;
plot(servoAngles(:, 1), 'x-');
plot(servoAngles(:, 2), 'o-');
legend('\theta_1', '\theta_2');
xlabel('Position Index');
ylabel('Angle (deg)');
print(gcf(), 'angles.png', '-dpng', '-r300');