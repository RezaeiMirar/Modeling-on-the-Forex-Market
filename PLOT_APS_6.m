function [y5, s_1_1, s_1_2, THETA_1, s_2_1, s_2_2, THETA_2, lambda, a4] ...
    = PLOT_APS_6(m_1, M, lambda, tt, b, S1, S2, S1_1, x, a1, a0, r) % PLOT_APS_4
%--------------------------------  (a_1,w_1)  -----------------------------
[solution_1,y2,Maximum] = SELECT_7(m_1, M, lambda, tt, b, x, r, S1_1, S2);
while Maximum > 2
    lambda = lambda - 0.1;
    [solution_1,y2,Maximum] = SELECT_7(m_1, M, lambda, tt, b, x, r,S1_1, S2);
end

%---------------------------------  Theta_1  ------------------------------
fix_1 = solution_1(2) * tt / 4;
Theta_1(1:m_1*M) = [ linspace(-fix_1(1),fix_1(1),M), linspace(-fix_1(2),fix_1(2),M), ...
    linspace(-fix_1(3),fix_1(3),M), linspace(-fix_1(4),fix_1(4),M), ...
    linspace(-fix_1(5),fix_1(5),M), linspace(-fix_1(6),fix_1(6),M) ];

k2 = 0;
for i = 1 : m_1*M
    k2 = k2 + 1;
    Index_3(k2,1) = Theta_1(i);
    SUM_1(k2) = 0;
    for j = 1 : m_1
        S32(j) = abs( S1(b(j), tt(j), a1, a0) - solution_1(1) * sin(solution_1(2)*tt(j)+Index_3(k2,1)) );
        SUM_1(k2) = SUM_1(k2) + S32(j);
    end
end
Minimum_1 = min(SUM_1);
[Index_6,Index_7] = find(SUM_1==Minimum_1);
THETA_1 = Index_3(Index_7(1),1); % Solution Set

%--------------------------------  (a_2,w_2)  -----------------------------
%disp('w_2 :')
a4 = 1;
S3=@(y, t) y-a1*t-a0-solution_1(1)*sin(solution_1(2)*t+THETA_1);

[solution_22,Maximum2] = SELECT_6(a4, solution_1, M, S3, S2, b, tt, m_1, r, x, THETA_1);
k5 = 1;
solution_2(:,k5) = solution_22';
Maximum(k5) = Maximum2;
while (Maximum2 < 5) && (a4 < 20)
    a4 = a4 + 1;
    [solution_22,Maximum2] = SELECT_6(a4, solution_1, M, S3, S2, b, tt, m_1, r, x, THETA_1);
    k5 = k5 + 1;
    solution_2(:,k5) = solution_22';
    Maximum(k5) = Maximum2;
end
if k5 == 1
    k6 = k5;
% elseif (k5 > 1) && (Maximum2 == 5)
%     k6 = k5;
else
    k6 = k5-1;
end
%---------------------------------  Theta_2  ------------------------------
fix_2 = solution_2(2,k6) * tt / 4;
Theta_2(1:m_1*M) = [ linspace(-fix_2(1),fix_2(1),M), linspace(-fix_2(2),fix_2(2),M), ...
    linspace(-fix_2(3),fix_2(3),M), linspace(-fix_2(4),fix_2(4),M), ...
    linspace(-fix_2(5),fix_2(5),M), linspace(-fix_2(6),fix_2(6),M) ];

k4 = 0;
for i = 1 : m_1*M
    k4 = k4 + 1;
    Index_11(k4,1) = Theta_2(i);
    SUM_3(k4) = 0;
    for j = 1 : m_1
        S34(j) = abs( S3(b(j), tt(j)) - solution_2(1,k6) * sin(solution_2(2,k6)*tt(j)+Index_11(k4,1)) );
        SUM_3(k4) = SUM_3(k4) + S34(j);
    end
end
Minimum_3 = min(SUM_3);
[Index_12,Index_13] = find(SUM_3==Minimum_3);
THETA_2 = Index_11(Index_13(1),1); % Solution Set

y5 = r + solution_1(1)*sin(solution_1(2)*x+THETA_1) + solution_2(1,k6)*sin(solution_2(2,k6)*x+THETA_2);

s_1_1 = solution_1(1);
s_1_2 = solution_1(2);
s_2_1 = solution_2(1,k6);
s_2_2 = solution_2(2,k6);

end

