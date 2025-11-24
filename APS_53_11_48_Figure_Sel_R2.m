% Anti Parallel Solving (APS) Method; Average Points; High accuracy
clear all
%clc
%close all
k = 2;
m_1 = 3 * k;
M = 5000; % Resolution
m22 = 99999;%1759, 3500, 4127, 7913, 8972, 9536
NSZO = 0; % To make the sequence; Sequence Index of Zero and One
NSZO_R1 = 0; % For function, RR=1
NSZO_R2 = 0; % For function, RR=1
NSZO_R3 = 0; % For function, RR=3
NSZO_Rev_R1 = 0; % Reverse, RR=1
NSZO_Rev_R2 = 0; % Reverse, RR=2
NSZO_Rev_R3 = 0; % Reverse, RR=3
ii3 = 0;
disp('--------------------------------')
kk = input('Please enter a number from 1-22: '); % High=A(:,2); Low=A(:,3); Close=A(:,4)

tic
[A, Text, Row, PipValue, Digit, Dist_HC, Stop_HD] = Select_Excel(kk);
Typ_Pri = ( A(1:m22,2) + A(1:m22,3) + A(1:m22,4) ) / 3; % Typical Price
EMA_11 = A(1,6); % 12 candels
EMA_22 = A(1,7); % 26 candels
MM = size(Typ_Pri,1);
[EXPRICE, Spark, Delay] = EXPORT_Typical_5(MM, Typ_Pri, A(1:m22,8), EMA_11, EMA_22);
b_T = EXPRICE(:,1); % Extreme prices
tt_T = EXPRICE(:,3); % Discrete time
for i = 1 : size(EXPRICE,1)
    TIME(i,1) = Row(EXPRICE(i,3),1); % Continuous time
end
SIZE = length(b_T);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
figure()
plot(A(1:m22,8),Typ_Pri,'k')
hold on
plot(EXPRICE(:,3),EXPRICE(:,1),'mo','linewidth',2)
legend('Price','Extreme')
hold on
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%---------------------------------  ii = 1  -------------------------------
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if SIZE < m_1 + 2
    disp('The number of extreme points is not enough. It is')
    disp(SIZE)
else
for i = 1 : m_1 + 1
    Sol_Total(1,i*2-1) = tt_T(i); % Discrete time
    Sol_Total(1,i*2) = b_T(i); % Price
    Data_DT2(1,i) = TIME(i);
    Time_Spark(1,i) = Spark(i);
    Time_Delay(1,i) = Delay(i);
end

tt = tt_T(1:6);
b = b_T(1:6);

x = A(tt(1) : tt(6),8); % Discrete time
y = Typ_Pri(tt(1):tt(6)); % Typical Price
[a1,a0] = regression(x,y);
%--------------------------------  (a_1,w_1)  -----------------------------
ttt(1:m_1-1) = ( tt(1:m_1-1) + tt(2:m_1) ) / 2;
bb(1:m_1-1) = ( b(1:m_1-1) + b(2:m_1) ) / 2;

S1=@(y, t, m, n) y-m*t-n;
S2=@(t, w) sin(w * t);
S1_1=@(y, t) S1(y, t, a1, a0);

r = a1 * x + a0;

lambda = 2;
[solution_1,Maximum] = SELECT_5(m_1, M, lambda, ttt, bb, x, r, S1_1, S2);
while Maximum > 2
    lambda = lambda - 0.1;
    [solution_1,Maximum] = SELECT_5(m_1, M, lambda, ttt, bb, x, r,S1_1, S2);
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if solution_1(2) < 10^(-4)
    HasMid = 0;
    lambda = 2;
    [y5, s_1_1, s_1_2, THETA_1, s_2_1, s_2_2, THETA_2, lambda, a4] = ...
        PLOT_APS_6(m_1, M, lambda, tt, b, S1, S2, S1_1, x, a1, a0, r);
    Sol_Total(1,15:23) = [s_1_1, s_1_2, THETA_1, s_2_1, s_2_2, THETA_2, lambda, a4-1, HasMid];
else
    HasMid = 1;
    %---------------------------------  Theta_1  --------------------------
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
            S32(j) = abs( S1(b(j), tt(j), a1, a0) - solution_1(1) * ...
                sin(solution_1(2)*tt(j)+Index_3(k2,1)) );
            SUM_1(k2) = SUM_1(k2) + S32(j);
        end
    end
    Minimum_1 = min(SUM_1);
    [Index_6,Index_7] = find(SUM_1==Minimum_1);
    THETA_1 = Index_3(Index_7(1),1); % Solution Set
    
    %--------------------------------  (a_2,w_2)  -------------------------
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
    else
        k6 = k5-1;
    end
    %---------------------------------  Theta_2  --------------------------
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
            S34(j) = abs( S3(b(j), tt(j)) - solution_2(1,k6) * ...
                sin(solution_2(2,k6)*tt(j)+Index_11(k4,1)) );
            SUM_3(k4) = SUM_3(k4) + S34(j);
        end
    end
    Minimum_3 = min(SUM_3);
    [Index_12,Index_13] = find(SUM_3==Minimum_3);
    THETA_2 = Index_11(Index_13(1),1); % Solution Set
    
    y5 = r + solution_1(1)*sin(solution_1(2)*x+THETA_1) + ...
        solution_2(1,k6)*sin(solution_2(2,k6)*x+THETA_2);
    
    Sol_Total(1,15:23) = [solution_1(1), solution_1(2), THETA_1, solution_2(1,k6), ...
        solution_2(2,k6), THETA_2, lambda, a4-1, HasMid];
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Sign Regression  >>>>>>>>>>>>>>>>>>>>>>>>>>
Sig_Reg(1) = sign(a1);

Slope(1,1) = ( b_T(3) - b_T(1) ) / ( tt_T(3) - tt_T(1) ); % 1_3
i3 = 0;
i4 = 0;
Slo_2 = {};
Slo_7 = zeros();
if sign(Slope(1,1)) == Sig_Reg(1)
    Slo_5(1) = 1; % Additional
    for i = 1 : 6
        if i == 1
            continue
        elseif i == 3
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(1,1) * ( tt_T(i) - tt_T(1) ) + b_T(1);
            if y_prim <= b_T(i) && b_T(1) < b_T(2) % min
                Slo_1(1,i3) = 1;
                Slo_6(1,1) = 1001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(1) > b_T(2) % max
                Slo_1(1,i3) = 1;
                Slo_6(1,1) = 1002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(1,i3) = 0;
                Slo_6(1,1) = 1005; % Additional
            end
        end
    end
    Slo_3(1,1) = CONT_R_2(Slo_1(1,:), size(Slo_1,2));
    if Slo_3(1,1) == 4 && Sig_Reg(1) == -1
        i4 = i4 + 1;
        Slo_2(1,1) = {'1_3'};
        Slo_2(1,2) = {'4'};
        Slo_2(1,3) = {'Descending'};
        Slo_2(1,4) = {'1'};
        Slo_2(1,5) = Slo_9;
        Slo_2(1,6) = Slo_10;
        Slo_2(1,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,1));
        Slo_7(i4,2) = 1;
    elseif Slo_3(1,1) == 4 && Sig_Reg(1) == 1
        i4 = i4 + 1;
        Slo_2(1,1) = {'1_3'};
        Slo_2(1,2) = {'4'};
        Slo_2(1,3) = {'Ascending'};
        Slo_2(1,4) = {'2'};
        Slo_2(1,5) = Slo_9;
        Slo_2(1,6) = Slo_10;
        Slo_2(1,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,1));
        Slo_7(i4,2) = 1;
    end
end
%--------------------------------------------------------------------------
Slope(1,2) = ( b_T(4) - b_T(2) ) / ( tt_T(4) - tt_T(2) ); % 2_4
i3 = 0;
if sign(Slope(1,2)) == Sig_Reg(1)
    Slo_5(2) = 2; % Additional
    for i = 2 : 6
        if i == 2
            continue
        elseif i == 4
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(1,2) * ( tt_T(i) - tt_T(2) ) + b_T(2);
            if y_prim >= b_T(i) && b_T(1) < b_T(2) % max
                Slo_1(2,i3) = 1;
                Slo_6(1,2) = 2001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(1) > b_T(2) % min
                Slo_1(2,i3) = 1;
                Slo_6(1,2) = 2002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(2,i3) = 0;
                Slo_6(1,2) = 2005; % Additional
            end
        end
    end
    Slo_3(1,2) = CONT_R_2(Slo_1(2,:), size(Slo_1,2));
    if Slo_3(1,2) == 3 && Sig_Reg(1) == -1
        i4 = i4 + 1;
        Slo_2(2,1) = {'2_4'};
        Slo_2(2,2) = {'3'};
        Slo_2(2,3) = {'Descending'};
        Slo_2(2,4) = {'1'};
        Slo_2(2,5) = Slo_9;
        Slo_2(2,6) = Slo_10;
        Slo_2(2,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,2));
        Slo_7(i4,2) = 2;
    elseif Slo_3(1,2) == 3 && Sig_Reg(1) == 1
        i4 = i4 + 1;
        Slo_2(2,1) = {'2_4'};
        Slo_2(2,2) = {'3'};
        Slo_2(2,3) = {'Ascending'};
        Slo_2(2,4) = {'2'};
        Slo_2(2,5) = Slo_9;
        Slo_2(2,6) = Slo_10;
        Slo_2(2,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,2));
        Slo_7(i4,2) = 2;
    end
end
%--------------------------------------------------------------------------
Slope(1,3) = ( b_T(5) - b_T(3) ) / ( tt_T(5) - tt_T(3) ); % 3_5
i3 = 0;
if sign(Slope(1,3)) == Sig_Reg(1)
    Slo_5(3) = 3; % Additional
    for i = 3 : 6
        if i == 3
            continue
        elseif i == 5
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(1,3) * ( tt_T(i) - tt_T(3) ) + b_T(3);
            if y_prim <= b_T(i) && b_T(1) < b_T(2) % min
                Slo_1(3,i3) = 1;
                Slo_6(1,3) = 3001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(1) > b_T(2) % max
                Slo_1(3,i3) = 1;
                Slo_6(1,3) = 3002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(3,i3) = 0;
                Slo_6(1,3) = 3005; % Additional
            end
        end
    end
    Slo_3(1,3) = CONT_R_2(Slo_1(3,:), size(Slo_1,2));
    if Slo_3(1,3) == 2 && Sig_Reg(1) == -1
        i4 = i4 + 1;
        Slo_2(3,1) = {'3_5'};
        Slo_2(3,2) = {'2'};
        Slo_2(3,3) = {'Descending'};
        Slo_2(3,4) = {'1'};
        Slo_2(3,5) = Slo_9;
        Slo_2(3,6) = Slo_10;
        Slo_2(3,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,3));
        Slo_7(i4,2) = 3;
    elseif Slo_3(1,3) == 2 && Sig_Reg(1) == 1
        i4 = i4 + 1;
        Slo_2(3,1) = {'3_5'};
        Slo_2(3,2) = {'2'};
        Slo_2(3,3) = {'Ascending'};
        Slo_2(3,4) = {'2'};
        Slo_2(3,5) = Slo_9;
        Slo_2(3,6) = Slo_10;
        Slo_2(3,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,3));
        Slo_7(i4,2) = 3;
    end
end
%--------------------------------------------------------------------------
Slope(1,4) = ( b_T(6) - b_T(4) ) / ( tt_T(6) - tt_T(4) ); % 4_6
i3 = 0;
if sign(Slope(1,4)) == Sig_Reg(1)
    Slo_5(4) = 4; % Additional
    for i = 4 : 6
        if i == 4
            continue
        elseif i == 6
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(1,4) * ( tt_T(i) - tt_T(4) ) + b_T(4);
            if y_prim >= b_T(i) && b_T(1) < b_T(2) % max
                Slo_1(4,i3) = 1;
                Slo_6(1,4) = 4001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(1) > b_T(2) % min
                Slo_1(4,i3) = 1;
                Slo_6(1,4) = 4002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(4,i3) = 0;
                Slo_6(1,4) = 4005; % Additional
            end
        end
    end
    Slo_3(1,4) = CONT_R_2(Slo_1(4,:), size(Slo_1,2));
    if Slo_3(1,4) == 1 && Sig_Reg(1) == -1
        i4 = i4 + 1;
        Slo_2(4,1) = {'4_6'};
        Slo_2(4,2) = {'1'};
        Slo_2(4,3) = {'Descending'};
        Slo_2(4,4) = {'1'};
        Slo_2(4,5) = Slo_9;
        Slo_2(4,6) = Slo_10;
        Slo_2(4,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,4));
        Slo_7(i4,2) = 4;
    elseif Slo_3(1,4) == 1 && Sig_Reg(1) == 1
        i4 = i4 + 1;
        Slo_2(4,1) = {'4_6'};
        Slo_2(4,2) = {'1'};
        Slo_2(4,3) = {'Ascending'};
        Slo_2(4,4) = {'2'};
        Slo_2(4,5) = Slo_9;
        Slo_2(4,6) = Slo_10;
        Slo_2(4,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,4));
        Slo_7(i4,2) = 4;
    end
end
%--------------------------------------------------------------------------
Slope(1,5) = ( b_T(5) - b_T(1) ) / ( tt_T(5) - tt_T(1) ); % 1_5
i3 = 0;
if sign(Slope(1,5)) == Sig_Reg(1)
    Slo_5(5) = 1; % Additional
    for i = 1 : 6
        if i == 1
            continue
        elseif i == 5
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(1,5) * ( tt_T(i) - tt_T(1) ) + b_T(1);
            if y_prim <= b_T(i) && b_T(1) < b_T(2) % min
                Slo_1(5,i3) = 1;
                Slo_6(1,5) = 5001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(1) > b_T(2) % max
                Slo_1(5,i3) = 1;
                Slo_6(1,5) = 5002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(5,i3) = 0;
                Slo_6(1,5) = 5005; % Additional
            end
        end
    end
    Slo_3(1,5) = CONT_R_2(Slo_1(5,:), size(Slo_1,2));
    if Slo_3(1,5) == 4 && Sig_Reg(1) == -1
        i4 = i4 + 1;
        Slo_2(5,1) = {'1_5'};
        Slo_2(5,2) = {'4'};
        Slo_2(5,3) = {'Descending'};
        Slo_2(5,4) = {'1'};
        Slo_2(5,5) = Slo_9;
        Slo_2(5,6) = Slo_10;
        Slo_2(5,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,5));
        Slo_7(i4,2) = 5;
    elseif Slo_3(1,5) == 4 && Sig_Reg(1) == 1
        i4 = i4 + 1;
        Slo_2(5,1) = {'1_5'};
        Slo_2(5,2) = {'4'};
        Slo_2(5,3) = {'Ascending'};
        Slo_2(5,4) = {'2'};
        Slo_2(5,5) = Slo_9;
        Slo_2(5,6) = Slo_10;
        Slo_2(5,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,5));
        Slo_7(i4,2) = 5;
    end
end
%--------------------------------------------------------------------------
Slope(1,6) = ( b_T(6) - b_T(2) ) / ( tt_T(6) - tt_T(2) ); % 2_6
i3 = 0;
if sign(Slope(1,6)) == Sig_Reg(1)
    Slo_5(6) = 2; % Additional
    for i = 2 : 6
        if i == 2
            continue
        elseif i == 6
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(1,6) * ( tt_T(i) - tt_T(2) ) + b_T(2);
            if y_prim >= b_T(i) && b_T(1) < b_T(2) % max
                Slo_1(6,i3) = 1;
                Slo_6(1,6) = 6001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(1) > b_T(2) % min
                Slo_1(6,i3) = 1;
                Slo_6(1,6) = 6002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(6,i3) = 0;
                Slo_6(1,6) = 6005; % Additional
            end
        end
    end
    Slo_3(1,6) = CONT_R_2(Slo_1(6,:), size(Slo_1,2));
    if Slo_3(1,6) == 3 && Sig_Reg(1) == -1
        i4 = i4 + 1;
        Slo_2(6,1) = {'2_6'};
        Slo_2(6,2) = {'3'};
        Slo_2(6,3) = {'Descending'};
        Slo_2(6,4) = {'1'};
        Slo_2(6,5) = Slo_9;
        Slo_2(6,6) = Slo_10;
        Slo_2(6,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,6));
        Slo_7(i4,2) = 6;
    elseif Slo_3(1,6) == 3 && Sig_Reg(1) == 1
        i4 = i4 + 1;
        Slo_2(6,1) = {'2_6'};
        Slo_2(6,2) = {'3'};
        Slo_2(6,3) = {'Ascending'};
        Slo_2(6,4) = {'2'};
        Slo_2(6,5) = Slo_9;
        Slo_2(6,6) = Slo_10;
        Slo_2(6,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(1,6));
        Slo_7(i4,2) = 6;
    end
end
%--------------------------------------------------------------------------

if Slo_7 == 0
    Trend_Num_1(1,1) = {'0'};
    Point_Num_1(1,1) = {'0'};
    Trend_1(1,1) = {'No Terend'};
    Trend_Mac_1(1,1) = {'0'};
    Up_Down_1(1,1) = {'No'};
    Up_Down_Mac_1(1,1) = {'0'};
    Min_Max_1(1,1) = {'No'};
else
    Slo_4 = min(Slo_7(:,1));
    ind_slo = find( Slo_7(:,1) == Slo_4 );
    Slo_8 = Slo_7(ind_slo,2);
    
    for i = 1 : length(ind_slo)
        plot(x,Slope(1,Slo_8) * (x-tt_T(Slo_5(Slo_8))) ...
            + b_T(Slo_5(Slo_8)),'b') % Additional
        hold on
    end
    
    Trend_Num_1(1,1) = Slo_2(Slo_8,1);
    Point_Num_1(1,1) = Slo_2(Slo_8,2);
    Trend_1(1,1) = Slo_2(Slo_8,3);
    Trend_Mac_1(1,1) = Slo_2(Slo_8,4);
    Up_Down_1(1,1) = Slo_2(Slo_8,5);
    Up_Down_Mac_1(1,1) = Slo_2(Slo_8,6);
    Min_Max_1(1,1) = Slo_2(Slo_8,7);
end

plot(x,r,'g')
hold on

% Slo_2
% Slo_7(:,1)
% Slo_7(:,2)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%---------------------------------  ii = 2  -------------------------------
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

tt = tt_T(2:7);
b = b_T(2:7);

x = A(tt(1) : tt(6),8); % Discrete time
y = Typ_Pri(tt_T(2):tt_T(7)); % Typical Price

a1_Previous = a1;
a0_Previous = a0;
[a1,a0] = regression(x,y);
%--------------------------------  (a_1,w_1)  -----------------------------
ttt(1:m_1-1) = ( tt(1:m_1-1) + tt(2:m_1) ) / 2;
bb(1:m_1-1) = ( b(1:m_1-1) + b(2:m_1) ) / 2;

S1=@(y, t, m, n) y-m*t-n;
S2=@(t, w) sin(w * t);
S1_1=@(y, t) S1(y, t, a1, a0);

r = a1 * x + a0;

lambda = 2;
[solution_1,Maximum] = SELECT_5(m_1, M, lambda, ttt, bb, x, r, S1_1, S2);
while Maximum > 2
    lambda = lambda - 0.1;
    [solution_1,Maximum] = SELECT_5(m_1, M, lambda, ttt, bb, x, r,S1_1, S2);
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if solution_1(2) < 10^(-4)
    HasMid = 0;
    lambda = 2;
    [y5, s_1_1, s_1_2, THETA_1, s_2_1, s_2_2, THETA_2, lambda, a4] = ...
        PLOT_APS_6(m_1, M, lambda, tt, b, S1, S2, S1_1, x, a1, a0, r);
    Sol_Total(1,24:32) = [s_1_1, s_1_2, THETA_1, s_2_1, s_2_2, THETA_2, lambda, a4-1, HasMid];
else
    HasMid = 1;
    %---------------------------------  Theta_1  --------------------------
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
            S32(j) = abs( S1(b(j), tt(j), a1, a0) - solution_1(1) * ...
                sin(solution_1(2)*tt(j)+Index_3(k2,1)) );
            SUM_1(k2) = SUM_1(k2) + S32(j);
        end
    end
    Minimum_1 = min(SUM_1);
    [Index_6,Index_7] = find(SUM_1==Minimum_1);
    THETA_1 = Index_3(Index_7(1),1); % Solution Set
    
    %--------------------------------  (a_2,w_2)  -------------------------
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
    else
        k6 = k5-1;
    end
    %---------------------------------  Theta_2  --------------------------
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
            S34(j) = abs( S3(b(j), tt(j)) - solution_2(1,k6) * ...
                sin(solution_2(2,k6)*tt(j)+Index_11(k4,1)) );
            SUM_3(k4) = SUM_3(k4) + S34(j);
        end
    end
    Minimum_3 = min(SUM_3);
    [Index_12,Index_13] = find(SUM_3==Minimum_3);
    THETA_2 = Index_11(Index_13(1),1); % Solution Set
    
    y5 = r + solution_1(1)*sin(solution_1(2)*x+THETA_1) + ...
        solution_2(1,k6)*sin(solution_2(2,k6)*x+THETA_2);
    
    Sol_Total(1,24:32) = [solution_1(1), solution_1(2), THETA_1, solution_2(1,k6), ...
        solution_2(2,k6), THETA_2, lambda, a4-1, HasMid];
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Sign Regression  >>>>>>>>>>>>>>>>>>>>>>>>>>
Sig_Reg(2) = sign(a1);

Slope(2,1) = ( b_T(4) - b_T(2) ) / ( tt_T(4) - tt_T(2) ); % 1_3
i3 = 0;
i4 = 0;
Slo_2 = {};
Slo_7 = zeros();
if sign(Slope(2,1)) == Sig_Reg(2)
    Slo_5(1) = 2; % Additional
    for i = 2 : 7
        if i == 2
            continue
        elseif i == 4
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(2,1) * ( tt_T(i) - tt_T(2) ) + b_T(2);
            if y_prim <= b_T(i) && b_T(2) < b_T(3) % min
                Slo_1(1,i3) = 1;
                Slo_6(2,1) = 1001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(2) > b_T(3) % max
                Slo_1(1,i3) = 1;
                Slo_6(2,1) = 1002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(1,i3) = 0;
                Slo_6(2,1) = 1005; % Additional
            end
        end
    end
    Slo_3(2,1) = CONT_R_2(Slo_1(1,:), size(Slo_1,2));
    if Slo_3(2,1) == 4 && Sig_Reg(2) == -1
        i4 = i4 + 1;
        Slo_2(1,1) = {'1_3'};
        Slo_2(1,2) = {'4'};
        Slo_2(1,3) = {'Descending'};
        Slo_2(1,4) = {'1'};
        Slo_2(1,5) = Slo_9;
        Slo_2(1,6) = Slo_10;
        Slo_2(1,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,1));
        Slo_7(i4,2) = 1;
    elseif Slo_3(2,1) == 4 && Sig_Reg(2) == 1
        i4 = i4 + 1;
        Slo_2(1,1) = {'1_3'};
        Slo_2(1,2) = {'4'};
        Slo_2(1,3) = {'Ascending'};
        Slo_2(1,4) = {'2'};
        Slo_2(1,5) = Slo_9;
        Slo_2(1,6) = Slo_10;
        Slo_2(1,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,1));
        Slo_7(i4,2) = 1;
    end
end
%--------------------------------------------------------------------------
Slope(2,2) = ( b_T(5) - b_T(3) ) / ( tt_T(5) - tt_T(3) ); % 2_4
i3 = 0;
if sign(Slope(2,2)) == Sig_Reg(2)
    Slo_5(2) = 3; % Additional
    for i = 3 : 7
        if i == 3
            continue
        elseif i == 5
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(2,2) * ( tt_T(i) - tt_T(3) ) + b_T(3);
            if y_prim >= b_T(i) && b_T(2) < b_T(3) % max
                Slo_1(2,i3) = 1;
                Slo_6(2,2) = 2001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(2) > b_T(3) % min
                Slo_1(2,i3) = 1;
                Slo_6(2,2) = 2002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(2,i3) = 0;
                Slo_6(2,2) = 2005; % Additional
            end
        end
    end
    Slo_3(2,2) = CONT_R_2(Slo_1(2,:), size(Slo_1,2));
    if Slo_3(2,2) == 3 && Sig_Reg(2) == -1
        i4 = i4 + 1;
        Slo_2(2,1) = {'2_4'};
        Slo_2(2,2) = {'3'};
        Slo_2(2,3) = {'Descending'};
        Slo_2(2,4) = {'1'};
        Slo_2(2,5) = Slo_9;
        Slo_2(2,6) = Slo_10;
        Slo_2(2,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,2));
        Slo_7(i4,2) = 2;
    elseif Slo_3(2,2) == 3 && Sig_Reg(2) == 1
        i4 = i4 + 1;
        Slo_2(2,1) = {'2_4'};
        Slo_2(2,2) = {'3'};
        Slo_2(2,3) = {'Ascending'};
        Slo_2(2,4) = {'2'};
        Slo_2(2,5) = Slo_9;
        Slo_2(2,6) = Slo_10;
        Slo_2(2,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,2));
        Slo_7(i4,2) = 2;
    end
end
%--------------------------------------------------------------------------
Slope(2,3) = ( b_T(6) - b_T(4) ) / ( tt_T(6) - tt_T(4) ); % 3_5
i3 = 0;
if sign(Slope(2,3)) == Sig_Reg(2)
    Slo_5(3) = 4; % Additional
    for i = 4 : 7
        if i == 4
            continue
        elseif i == 6
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(2,3) * ( tt_T(i) - tt_T(4) ) + b_T(4);
            if y_prim <= b_T(i) && b_T(2) < b_T(3) % min
                Slo_1(3,i3) = 1;
                Slo_6(2,3) = 3001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(2) > b_T(3) % max
                Slo_1(3,i3) = 1;
                Slo_6(2,3) = 3002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(3,i3) = 0;
                Slo_6(2,3) = 3005; % Additional
            end
        end
    end
    Slo_3(2,3) = CONT_R_2(Slo_1(3,:), size(Slo_1,2));
    if Slo_3(2,3) == 2 && Sig_Reg(2) == -1
        i4 = i4 + 1;
        Slo_2(3,1) = {'3_5'};
        Slo_2(3,2) = {'2'};
        Slo_2(3,3) = {'Descending'};
        Slo_2(3,4) = {'1'};
        Slo_2(3,5) = Slo_9;
        Slo_2(3,6) = Slo_10;
        Slo_2(3,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,3));
        Slo_7(i4,2) = 3;
    elseif Slo_3(2,3) == 2 && Sig_Reg(2) == 1
        i4 = i4 + 1;
        Slo_2(3,1) = {'3_5'};
        Slo_2(3,2) = {'2'};
        Slo_2(3,3) = {'Ascending'};
        Slo_2(3,4) = {'2'};
        Slo_2(3,5) = Slo_9;
        Slo_2(3,6) = Slo_10;
        Slo_2(3,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,3));
        Slo_7(i4,2) = 3;
    end
end
%--------------------------------------------------------------------------
Slope(2,4) = ( b_T(7) - b_T(5) ) / ( tt_T(7) - tt_T(5) ); % 4_6
i3 = 0;
if sign(Slope(2,4)) == Sig_Reg(2)
    Slo_5(4) = 5; % Additional
    for i = 5 : 7
        if i == 5
            continue
        elseif i == 7
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(2,4) * ( tt_T(i) - tt_T(5) ) + b_T(5);
            if y_prim >= b_T(i) && b_T(2) < b_T(3) % max
                Slo_1(4,i3) = 1;
                Slo_6(2,4) = 4001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(2) > b_T(3) % min
                Slo_1(4,i3) = 1;
                Slo_6(2,4) = 4002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(4,i3) = 0;
                Slo_6(2,4) = 4005; % Additional
            end
        end
    end
    Slo_3(2,4) = CONT_R_2(Slo_1(4,:), size(Slo_1,2));
    if Slo_3(2,4) == 1 && Sig_Reg(2) == -1
        i4 = i4 + 1;
        Slo_2(4,1) = {'4_6'};
        Slo_2(4,2) = {'1'};
        Slo_2(4,3) = {'Descending'};
        Slo_2(4,4) = {'1'};
        Slo_2(4,5) = Slo_9;
        Slo_2(4,6) = Slo_10;
        Slo_2(4,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,4));
        Slo_7(i4,2) = 4;
    elseif Slo_3(2,4) == 1 && Sig_Reg(2) == 1
        i4 = i4 + 1;
        Slo_2(4,1) = {'4_6'};
        Slo_2(4,2) = {'1'};
        Slo_2(4,3) = {'Ascending'};
        Slo_2(4,4) = {'2'};
        Slo_2(4,5) = Slo_9;
        Slo_2(4,6) = Slo_10;
        Slo_2(4,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,4));
        Slo_7(i4,2) = 4;
    end
end
%--------------------------------------------------------------------------
Slope(2,5) = ( b_T(6) - b_T(2) ) / ( tt_T(6) - tt_T(2) ); % 1_5
i3 = 0;
if sign(Slope(2,5)) == Sig_Reg(2)
    Slo_5(5) = 2; % Additional
    for i = 2 : 7
        if i == 2
            continue
        elseif i == 7
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(2,5) * ( tt_T(i) - tt_T(2) ) + b_T(2);
            if y_prim <= b_T(i) && b_T(2) < b_T(3) % min
                Slo_1(5,i3) = 1;
                Slo_6(2,5) = 5001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(2) > b_T(3) % max
                Slo_1(5,i3) = 1;
                Slo_6(2,5) = 5002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(5,i3) = 0;
                Slo_6(2,5) = 5005; % Additional
            end
        end
    end
    Slo_3(2,5) = CONT_R_2(Slo_1(5,:), size(Slo_1,2));
    if Slo_3(2,5) == 4 && Sig_Reg(2) == -1
        i4 = i4 + 1;
        Slo_2(5,1) = {'1_5'};
        Slo_2(5,2) = {'4'};
        Slo_2(5,3) = {'Descending'};
        Slo_2(5,4) = {'1'};
        Slo_2(5,5) = Slo_9;
        Slo_2(5,6) = Slo_10;
        Slo_2(5,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,5));
        Slo_7(i4,2) = 5;
    elseif Slo_3(2,5) == 4 && Sig_Reg(2) == 1
        i4 = i4 + 1;
        Slo_2(5,1) = {'1_5'};
        Slo_2(5,2) = {'4'};
        Slo_2(5,3) = {'Ascending'};
        Slo_2(5,4) = {'2'};
        Slo_2(5,5) = Slo_9;
        Slo_2(5,6) = Slo_10;
        Slo_2(5,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,5));
        Slo_7(i4,2) = 5;
    end
end
%--------------------------------------------------------------------------
Slope(2,6) = ( b_T(7) - b_T(3) ) / ( tt_T(7) - tt_T(3) ); % 2_6
i3 = 0;
if sign(Slope(2,6)) == Sig_Reg(2)
    Slo_5(6) = 3; % Additional
    for i = 3 : 7
        if i == 3
            continue
        elseif i == 7
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(2,6) * ( tt_T(i) - tt_T(3) ) + b_T(3);
            if y_prim >= b_T(i) && b_T(2) < b_T(3) % max
                Slo_1(6,i3) = 1;
                Slo_6(2,6) = 6001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(2) > b_T(3) % min
                Slo_1(6,i3) = 1;
                Slo_6(2,6) = 6002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(6,i3) = 0;
                Slo_6(2,6) = 6005; % Additional
            end
        end
    end
    Slo_3(2,6) = CONT_R_2(Slo_1(6,:), size(Slo_1,2));
    if Slo_3(2,6) == 3 && Sig_Reg(2) == -1
        i4 = i4 + 1;
        Slo_2(6,1) = {'2_6'};
        Slo_2(6,2) = {'3'};
        Slo_2(6,3) = {'Descending'};
        Slo_2(6,4) = {'1'};
        Slo_2(6,5) = Slo_9;
        Slo_2(6,6) = Slo_10;
        Slo_2(6,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,6));
        Slo_7(i4,2) = 6;
    elseif Slo_3(2,6) == 3 && Sig_Reg(2) == 1
        i4 = i4 + 1;
        Slo_2(6,1) = {'2_6'};
        Slo_2(6,2) = {'3'};
        Slo_2(6,3) = {'Ascending'};
        Slo_2(6,4) = {'2'};
        Slo_2(6,5) = Slo_9;
        Slo_2(6,6) = Slo_10;
        Slo_2(6,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(2,6));
        Slo_7(i4,2) = 6;
    end
end
%--------------------------------------------------------------------------

if Slo_7 == 0
    Trend_Num_2(1,1) = {'0'};
    Point_Num_2(1,1) = {'0'};
    Trend_2(1,1) = {'No Terend'};
    Trend_Mac_2(1,1) = {'0'};
    Up_Down_2(1,1) = {'No'};
    Up_Down_Mac_2(1,1) = {'0'};
    Min_Max_2(1,1) = {'No'};
else
    Slo_4 = min(Slo_7(:,1));
    ind_slo = find( Slo_7(:,1) == Slo_4 );
    Slo_8 = Slo_7(ind_slo,2);
    
    for i = 1 : length(ind_slo)
        plot(x,Slope(2,Slo_8) * (x-tt_T(Slo_5(Slo_8))) ...
            + b_T(Slo_5(Slo_8)),'b') % Additional
        hold on
    end
    
    Trend_Num_2(1,1) = Slo_2(Slo_8,1);
    Point_Num_2(1,1) = Slo_2(Slo_8,2);
    Trend_2(1,1) = Slo_2(Slo_8,3);
    Trend_Mac_2(1,1) = Slo_2(Slo_8,4);
    Up_Down_2(1,1) = Slo_2(Slo_8,5);
    Up_Down_Mac_2(1,1) = Slo_2(Slo_8,6);
    Min_Max_2(1,1) = Slo_2(Slo_8,7);
end

if Sig_Reg(1) == Sig_Reg(2)
    Trend_Change(1,1) = 0;
else
    Trend_Change(1,1) = 1;
end

t7_t1(1,1) = Sol_Total(1,13) - Sol_Total(1,1) + 1;

plot(x,r,'m')
hold on

% Slo_2
% Slo_7(:,1)
% Slo_7(:,2)

%--------------------------------  Forecast  ------------------------------
First_t = tt_T(7);
NumerationE = 0;
DeltaY1 = a1_Previous * First_t + a0_Previous ...
    + Sol_Total(1,15) * sin(Sol_Total(1,16) * First_t + Sol_Total(1,17)) ...
    - a1 * First_t - a0 ...
    - Sol_Total(1,24) * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26));

while (NumerationE < 3) && ( First_t <= (tt_T(7) + 2 * t7_t1(1,1)) )
    First_t = First_t + 1;
    DeltaY2 = a1_Previous * First_t + a0_Previous ...
        + Sol_Total(1,15) * sin(Sol_Total(1,16) * First_t + Sol_Total(1,17)) ...
        - a1 * First_t - a0 ...
        - Sol_Total(1,24) * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26));
    
    if ( sign(DeltaY1) ~= sign(DeltaY2)) && (First_t <= MM)
        NumerationE = NumerationE + 1;
        NY1(1,NumerationE) = Row(First_t,1);
        NDeltaY(1,NumerationE) = First_t;
        PEncounter(1,NumerationE) = a1 * First_t + a0 + Sol_Total(1,24) ...
            * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26));
    elseif ( sign(DeltaY1) ~= sign(DeltaY2)) && (First_t > MM)
        NumerationE = NumerationE + 1;
        NY1(1,NumerationE) = {'Infinity'};
        NDeltaY(1,NumerationE) = First_t;
        PEncounter(1,NumerationE) = a1 * First_t + a0 + Sol_Total(1,24) ...
            * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26));
    end
    DeltaY1 = DeltaY2;
end

if (NumerationE==0) || (NumerationE==1) || (NumerationE==2) % The # encounters is not enough.
    First_encounter_CT = {'Infinity'};
    First_encounter_DT = -1;
    Second_encounter_CT = {'Infinity'};
    Second_encounter_DT = -1;
    Third_encounter_CT = {'Infinity'};
    Third_encounter_DT = -1;
    
    First_PEncounter = -1;
    Second_PEncounter = -1;
    Third_PEncounter = -1;
    
    XNUD_Encounter = {'NO'};
    Hunt_Encounter = {'NO'};
    
    MaxTyp_PriEncounter = -1;
    MinTyp_PriEncounter = -1;
    x_maxE = -1;
    x_minE = -1;
    
    Max_PercentE = -1;
    Deviation_E = -1;
    
    TypePosition_E = -1;
    TypePos_MaxPer_Dev_E = -1;
    TypePos_E = {'NO'};
elseif NDeltaY(1,1) > MM % Because there is no real data.
    First_encounter_CT = {'Infinity'};
    First_encounter_DT = -1;
    Second_encounter_CT = {'Infinity'};
    Second_encounter_DT = -1;
    Third_encounter_CT = {'Infinity'};
    Third_encounter_DT = -1;
    
    First_PEncounter = -1;
    Second_PEncounter = -1;
    Third_PEncounter = -1;
    
    XNUD_Encounter = {'NO'};
    Hunt_Encounter = {'NO'};
    
    MaxTyp_PriEncounter = -1;
    MinTyp_PriEncounter = -1;
    x_maxE = -1;
    x_minE = -1;
    
    Max_PercentE = -1;
    Deviation_E = -1;
    
    TypePosition_E = -1;
    TypePos_MaxPer_Dev_E = -1;
    TypePos_E = {'NO'};
elseif Spark(7) == 0 % The extreme point is not clear.
    First_encounter_CT = {'Infinity'};
    First_encounter_DT = -1;
    Second_encounter_CT = {'Infinity'};
    Second_encounter_DT = -1;
    Third_encounter_CT = {'Infinity'};
    Third_encounter_DT = -1;
    
    First_PEncounter = -1;
    Second_PEncounter = -1;
    Third_PEncounter = -1;
    
    XNUD_Encounter = {'NO'};
    Hunt_Encounter = {'NO'};
    
    MaxTyp_PriEncounter = -1;
    MinTyp_PriEncounter = -1;
    x_maxE = -1;
    x_minE = -1;
    
    Max_PercentE = -1;
    Deviation_E = -1;
    
    TypePosition_E = -1;
    TypePos_MaxPer_Dev_E = -1;
    TypePos_E = {'NO'};
elseif NDeltaY(1,1) <= Spark(7) % The first encounter is before the first spark.
    First_encounter_CT = {'Infinity'};
    First_encounter_DT = -1;
    Second_encounter_CT = {'Infinity'};
    Second_encounter_DT = -1;
    Third_encounter_CT = {'Infinity'};
    Third_encounter_DT = -1;
    
    First_PEncounter = -1;
    Second_PEncounter = -1;
    Third_PEncounter = -1;
    
    XNUD_Encounter = {'NO'};
    Hunt_Encounter = {'NO'};
    
    MaxTyp_PriEncounter = -1;
    MinTyp_PriEncounter = -1;
    x_maxE = -1;
    x_minE = -1;
    
    Max_PercentE = -1;
    Deviation_E = -1;
    
    TypePosition_E = -1;
    TypePos_MaxPer_Dev_E = -1;
    TypePos_E = {'NO'};
elseif (NDeltaY(1,1) - Spark(7)) < 1 % There is no data between first encounter and spark timing.
    First_encounter_CT = {'Infinity'};
    First_encounter_DT = -1;
    Second_encounter_CT = {'Infinity'};
    Second_encounter_DT = -1;
    Third_encounter_CT = {'Infinity'};
    Third_encounter_DT = -1;
    
    First_PEncounter = -1;
    Second_PEncounter = -1;
    Third_PEncounter = -1;
    
    XNUD_Encounter = {'NO'};
    Hunt_Encounter = {'NO'};
    
    MaxTyp_PriEncounter = -1;
    MinTyp_PriEncounter = -1;
    x_maxE = -1;
    x_minE = -1;
    
    Max_PercentE = -1;
    Deviation_E = -1;
    
    TypePosition_E = -1;
    TypePos_MaxPer_Dev_E = -1;
    TypePos_E = {'NO'};
elseif (NDeltaY(1,1) - tt_T(7) + 1) <= (t7_t1(1,1) / 2)
    First_encounter_CT = NY1(1,1);
    First_encounter_DT = NDeltaY(1,1);
    Second_encounter_CT = NY1(1,2);
    Second_encounter_DT = NDeltaY(1,2);
    Third_encounter_CT = NY1(1,3);
    Third_encounter_DT = NDeltaY(1,3);
    
    First_PEncounter = PEncounter(1,1);
    Second_PEncounter = PEncounter(1,2);
    Third_PEncounter = PEncounter(1,3);
    
    if NDeltaY(1,1) <= MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : NDeltaY(1,1) );
        MaxTyp_PriEncounter = max(MTyp_Pri); % Find max
        MinTyp_PriEncounter = min(MTyp_Pri); % Find min
        x_maxE = find(MTyp_Pri == MaxTyp_PriEncounter) + Spark(7);
        x_minE = find(MTyp_Pri == MinTyp_PriEncounter) + Spark(7);
    elseif NDeltaY(1,1) > MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : MM );
        MaxTyp_PriEncounter = max(MTyp_Pri); % Find max
        MinTyp_PriEncounter = min(MTyp_Pri); % Find min
        x_maxE = find(MTyp_Pri == MaxTyp_PriEncounter) + Spark(7);
        x_minE = find(MTyp_Pri == MinTyp_PriEncounter) + Spark(7);
    end
    
    if (b_T(6) < b_T(7)) && (b_T(7) < PEncounter(1,1)) % b_T(7) is max
        XNUD_Encounter(1,1) = {'XU'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) > PEncounter(1,1)) % b_T(7) is max
        XNUD_Encounter(1,1) = {'XD'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) < PEncounter(1,1)) % b_T(7) is min
        XNUD_Encounter(1,1) = {'NU'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) > PEncounter(1,1)) % b_T(7) is min
        XNUD_Encounter(1,1) = {'ND'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) == PEncounter(1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Encounter(1,1) = {'XE'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) == PEncounter(1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Encounter(1,1) = {'NE'};
    end
    
    if Typ_Pri(Spark(7)) < First_PEncounter % Hunt Up
        Hunt_Encounter = {'Hunt_U'};
        if MaxTyp_PriEncounter <= Typ_Pri(Spark(7))
            Max_PercentE = 0;
        elseif ( Typ_Pri(Spark(7)) < MaxTyp_PriEncounter ) ...
                && ( MaxTyp_PriEncounter < First_PEncounter )
            Max_PercentE = abs( MaxTyp_PriEncounter - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
        elseif MaxTyp_PriEncounter >= First_PEncounter
            Max_PercentE = 100;
        end
        if MinTyp_PriEncounter >= Typ_Pri(Spark(7))
            Deviation_E = 0;
        elseif MinTyp_PriEncounter < Typ_Pri(Spark(7))
            Deviation_E = abs( MinTyp_PriEncounter - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
        end
        %----------------------------------------------
        if (Max_PercentE >= 50) && (Deviation_E >= 100)
            if x_maxE < x_minE
                TypePosition_E = 1; % Max_Percent
            else
                TypePosition_E = 2; % Deviation
            end
            TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
            if Typ_Pri(Spark(7)) < Typ_Pri(NDeltaY(1,1))
                TypePos_E = {'Max_Percent'};
            else
                TypePos_E = {'Deviation'};
            end
        elseif Max_PercentE >= 50
            TypePosition_E = 1; % Max_Percent
            TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
            if Typ_Pri(Spark(7)) < Typ_Pri(NDeltaY(1,1))
                TypePos_E = {'Max_Percent'};
            else
                TypePos_E = {'Deviation'};
            end
        elseif Deviation_E >= 100
            TypePosition_E = 2; % Deviation
            TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
            if Typ_Pri(Spark(7)) < Typ_Pri(NDeltaY(1,1))
                TypePos_E = {'Max_Percent'};
            else
                TypePos_E = {'Deviation'};
            end
        else
            TypePosition_E = 3; % No mode
            if Typ_Pri(Spark(7)) < Typ_Pri(NDeltaY(1,1))
                TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                    abs( First_PEncounter - Typ_Pri(Spark(7)) );
                TypePos_E = {'Max_Percent'};
            else
                TypePos_MaxPer_Dev_E = -abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                    abs( First_PEncounter - Typ_Pri(Spark(7)) );
                TypePos_E = {'Deviation'};
            end
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark(7)) > First_PEncounter % Hunt Down
        Hunt_Encounter = {'Hunt_D'};
        if MinTyp_PriEncounter >= Typ_Pri(Spark(7))
            Max_PercentE = 0;
        elseif ( Typ_Pri(Spark(7)) > MinTyp_PriEncounter ) ...
                && ( MinTyp_PriEncounter > First_PEncounter )
            Max_PercentE = abs( MinTyp_PriEncounter - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
        elseif MinTyp_PriEncounter <= First_PEncounter
            Max_PercentE = 100;
        end
        if MaxTyp_PriEncounter <= Typ_Pri(Spark(7))
            Deviation_E = 0;
        elseif MaxTyp_PriEncounter > Typ_Pri(Spark(7))
            Deviation_E = abs( MaxTyp_PriEncounter - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
        end
        %----------------------------------------------
        if (Max_PercentE >= 50) && (Deviation_E >= 100)
            if x_maxE > x_minE
                TypePosition_E = 1; % Max_Percent
            else
                TypePosition_E = 2; % Deviation
            end
            TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
            if Typ_Pri(Spark(7)) > Typ_Pri(NDeltaY(1,1))
                TypePos_E = {'Max_Percent'};
            else
                TypePos_E = {'Deviation'};
            end
        elseif Max_PercentE >= 50
            TypePosition_E = 1; % Max_Percent
            TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
            if Typ_Pri(Spark(7)) > Typ_Pri(NDeltaY(1,1))
                TypePos_E = {'Max_Percent'};
            else
                TypePos_E = {'Deviation'};
            end
        elseif Deviation_E >= 100
            TypePosition_E = 2; % Deviation
            TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( First_PEncounter - Typ_Pri(Spark(7)) );
            if Typ_Pri(Spark(7)) > Typ_Pri(NDeltaY(1,1))
                TypePos_E = {'Max_Percent'};
            else
                TypePos_E = {'Deviation'};
            end
        else
            TypePosition_E = 3; % No mode
            if Typ_Pri(Spark(7)) > Typ_Pri(NDeltaY(1,1))
                TypePos_MaxPer_Dev_E = abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                    abs( First_PEncounter - Typ_Pri(Spark(7)) );
                TypePos_E = {'Max_Percent'};
            else
                TypePos_MaxPer_Dev_E = -abs( Typ_Pri(NDeltaY(1,1)) - Typ_Pri(Spark(7)) ) * 100 / ...
                    abs( First_PEncounter - Typ_Pri(Spark(7)) );
                TypePos_E = {'Deviation'};
            end
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark(7)) == First_PEncounter % Hunt Equal
        Hunt_Encounter = {'Hunt_E'};
        Max_PercentE = -1;
        Deviation_E = -1;
        
        TypePosition_E = -1;
        TypePos_MaxPer_Dev_E = -1;
        TypePos_E = {'NO'};
    end
else
    First_encounter_CT = {'Infinity'};
    First_encounter_DT = -1;
    Second_encounter_CT = {'Infinity'};
    Second_encounter_DT = -1;
    Third_encounter_CT = {'Infinity'};
    Third_encounter_DT = -1;
    
    First_PEncounter = -1;
    Second_PEncounter = -1;
    Third_PEncounter = -1;
    
    XNUD_Encounter = {'NO'};
    Hunt_Encounter = {'NO'};
    
    MaxTyp_PriEncounter = -1;
    MinTyp_PriEncounter = -1;
    x_maxE = -1;
    x_minE = -1;
    
    Max_PercentE = -1;
    Deviation_E = -1;
    
    TypePosition_E = -1;
    TypePos_MaxPer_Dev_E = -1;
    TypePos_E = {'NO'};
end

% Center of gravity
if NumerationE == 3
    Gravity_tE = ( NDeltaY(1,1) + NDeltaY(1,2) + NDeltaY(1,3) ) / 3;
    Gravity_PE = ( PEncounter(1,1) + PEncounter(1,2) + PEncounter(1,3) ) / 3;
else
    Gravity_tE = -1;
    Gravity_PE = -1;
end

if (NumerationE==0) || (NumerationE==1) || (NumerationE==2) % The # encounters is not enough.
    MaxTyp_PriceEG = -1;
    MinTyp_PriceEG = -1;
    x_maxEG = -1; % Additional
    x_minEG = -1; % Additional
    
    XNUD_EG(1,1) = {'NO'};
    Hunt_EG(1,1) = {'NO'};
    
    Max_PercentEG = -1;
    Deviation_EG = -1;
elseif Gravity_tE > MM % Because there is no real data.
    MaxTyp_PriceEG = -1;
    MinTyp_PriceEG = -1;
    x_maxEG = -1; % Additional
    x_minEG = -1; % Additional
    
    XNUD_EG(1,1) = {'NO'};
    Hunt_EG(1,1) = {'NO'};
    
    Max_PercentEG = -1;
    Deviation_EG = -1;
elseif Spark(7) == 0 % The extreme point is not clear.
    MaxTyp_PriceEG = -1;
    MinTyp_PriceEG = -1;
    x_maxEG = -1; % Additional
    x_minEG = -1; % Additional
    
    XNUD_EG(1,1) = {'NO'};
    Hunt_EG(1,1) = {'NO'};
    
    Max_PercentEG = -1;
    Deviation_EG = -1;
elseif Gravity_tE <= Spark(7) % The first encounter is before the first spark.
    MaxTyp_PriceEG = -1;
    MinTyp_PriceEG = -1;
    x_maxEG = -1; % Additional
    x_minEG = -1; % Additional
    
    XNUD_EG(1,1) = {'NO'};
    Hunt_EG(1,1) = {'NO'};
    
    Max_PercentEG = -1;
    Deviation_EG = -1;
elseif (Gravity_tE - Spark(7)) < 1 % There is no data between first encounter and spark timing.
    MaxTyp_PriceEG = -1;
    MinTyp_PriceEG = -1;
    x_maxEG = -1; % Additional
    x_minEG = -1; % Additional
    
    XNUD_EG(1,1) = {'NO'};
    Hunt_EG(1,1) = {'NO'};
    
    Max_PercentEG = -1;
    Deviation_EG = -1;
elseif (Gravity_tE - tt_T(7) + 1) <= (t7_t1(1,1) / 2)
    if Gravity_tE <= MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : Gravity_tE );
        MaxTyp_PriceEG = max(MTyp_Pri); % Find max
        MinTyp_PriceEG = min(MTyp_Pri); % Find min
        x_maxEG = find(MTyp_Pri == MaxTyp_PriceEG); % Additional
        x_minEG = find(MTyp_Pri == MinTyp_PriceEG); % Additional
    elseif Gravity_tE > MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : MM );
        MaxTyp_PriceEG = max(MTyp_Pri); % Find max
        MinTyp_PriceEG = min(MTyp_Pri); % Find min
        x_maxEG = find(MTyp_Pri == MaxTyp_PriceEG); % Additional
        x_minEG = find(MTyp_Pri == MinTyp_PriceEG); % Additional
    end
    
    if (b_T(6) < b_T(7)) && (b_T(7) < Gravity_PE) % b_T(7) is max
        XNUD_EG(1,1) = {'XU'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) > Gravity_PE) % b_T(7) is max
        XNUD_EG(1,1) = {'XD'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) < Gravity_PE) % b_T(7) is min
        XNUD_EG(1,1) = {'NU'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) > Gravity_PE) % b_T(7) is min
        XNUD_EG(1,1) = {'ND'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) == Gravity_PE) % b_T(7) is equal to PForecast(1,1).
        XNUD_EG(1,1) = {'XE'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) == Gravity_PE) % b_T(7) is equal to PForecast(1,1).
        XNUD_EG(1,1) = {'NE'};
    end
    
    if Typ_Pri(Spark(7)) < Gravity_PE % Hunt Up
        Hunt_EG = {'Hunt_U'};
        if MaxTyp_PriceEG <= Typ_Pri(Spark(7))
            Max_PercentEG = 0;
        elseif ( Typ_Pri(Spark(7)) < MaxTyp_PriceEG ) ...
                && ( MaxTyp_PriceEG < Gravity_PE )
            Max_PercentEG = abs( MaxTyp_PriceEG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PE - Typ_Pri(Spark(7)) );
        elseif MaxTyp_PriceEG >= Gravity_PE
            Max_PercentEG = 100;
        end
        if MinTyp_PriceEG >= Typ_Pri(Spark(7))
            Deviation_EG = 0;
        elseif MinTyp_PriceEG < Typ_Pri(Spark(7))
            Deviation_EG = abs( MinTyp_PriceEG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PE - Typ_Pri(Spark(7)) );
        end
    elseif Typ_Pri(Spark(7)) > Gravity_PE % Hunt Down
        Hunt_EG = {'Hunt_D'};
        if MinTyp_PriceEG >= Typ_Pri(Spark(7))
            Max_PercentEG = 0;
        elseif ( Typ_Pri(Spark(7)) > MinTyp_PriceEG ) ...
                && ( MinTyp_PriceEG > Gravity_PE )
            Max_PercentEG = abs( MinTyp_PriceEG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PE - Typ_Pri(Spark(7)) );
        elseif MinTyp_PriceEG <= Gravity_PE
            Max_PercentEG = 100;
        end
        if MaxTyp_PriceEG <= Typ_Pri(Spark(7))
            Deviation_EG = 0;
        elseif MaxTyp_PriceEG > Typ_Pri(Spark(7))
            Deviation_EG = abs( MaxTyp_PriceEG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PE - Typ_Pri(Spark(7)) );
        end
    elseif Typ_Pri(Spark(7)) == Gravity_PE % Hunt Equal
        Hunt_EG = {'Hunt_E'};
        Max_PercentEG = -1;
        Deviation_EG = -1;
    end
else
    MaxTyp_PriceEG = -1;
    MinTyp_PriceEG = -1;
    x_maxEG = -1; % Additional
    x_minEG = -1; % Additional
    
    XNUD_EG(1,1) = {'NO'};
    Hunt_EG(1,1) = {'NO'};
    
    Max_PercentEG = -1;
    Deviation_EG = -1;
end
%-------------------------------  With noise  -----------------------------
First_t = tt_T(7);
Numeration = 0;
DeltaY3 = a1_Previous * First_t + a0_Previous ...
    + Sol_Total(1,15) * sin(Sol_Total(1,16) * First_t + Sol_Total(1,17)) ...
    + Sol_Total(1,18) * sin(Sol_Total(1,19) * First_t + Sol_Total(1,20)) ...
    - a1 * First_t - a0 ...
    - Sol_Total(1,24) * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26)) ...
    - Sol_Total(1,27) * sin(Sol_Total(1,28) * First_t + Sol_Total(1,29));

while (Numeration < 3) && ( First_t <= (tt_T(7) + 2 * t7_t1(1,1)) )
    First_t = First_t + 1;
    DeltaY4 = a1_Previous * First_t + a0_Previous ...
        + Sol_Total(1,15) * sin(Sol_Total(1,16) * First_t + Sol_Total(1,17)) ...
        + Sol_Total(1,18) * sin(Sol_Total(1,19) * First_t + Sol_Total(1,20)) ...
        - a1 * First_t - a0 ...
        - Sol_Total(1,24) * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26)) ...
        - Sol_Total(1,27) * sin(Sol_Total(1,28) * First_t + Sol_Total(1,29));
    
    if ( sign(DeltaY3) ~= sign(DeltaY4)) && (First_t <= MM)
        Numeration = Numeration + 1;
        NY2(1,Numeration) = Row(First_t,1);
        YForecast(1,Numeration) = First_t;
        PForecast(1,Numeration) = a1 * First_t + a0 + Sol_Total(1,24) ...
            * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26)) + Sol_Total(1,27) ...
            * sin(Sol_Total(1,28) * First_t + Sol_Total(1,29));
        
        if (Spark(7) ~= 0) && (Typ_Pri(Spark(7)) <= PForecast(1,Numeration))
            PForecast_Rev(1,Numeration) = Typ_Pri(Spark(7)) - abs(Typ_Pri(Spark(7)) - PForecast(1,Numeration));
        elseif (Spark(7) ~= 0) && (Typ_Pri(Spark(7)) > PForecast(1,Numeration))
            PForecast_Rev(1,Numeration) = Typ_Pri(Spark(7)) + abs(Typ_Pri(Spark(7)) - PForecast(1,Numeration));
        elseif Spark(7) == 0
            PForecast_Rev(1,Numeration) = 0;
        end
    elseif ( sign(DeltaY3) ~= sign(DeltaY4)) && (First_t > MM)
        Numeration = Numeration + 1;
        NY2(1,Numeration) = {'Infinity'};
        YForecast(1,Numeration) = First_t;
        PForecast(1,Numeration) = a1 * First_t + a0 + Sol_Total(1,24) ...
            * sin(Sol_Total(1,25) * First_t + Sol_Total(1,26)) + Sol_Total(1,27) ...
            * sin(Sol_Total(1,28) * First_t + Sol_Total(1,29));
        
        if (Spark(7) ~= 0) && (Typ_Pri(Spark(7)) <= PForecast(1,Numeration))
            PForecast_Rev(1,Numeration) = Typ_Pri(Spark(7)) - abs(Typ_Pri(Spark(7)) - PForecast(1,Numeration));
        elseif (Spark(7) ~= 0) && (Typ_Pri(Spark(7)) > PForecast(1,Numeration))
            PForecast_Rev(1,Numeration) = Typ_Pri(Spark(7)) + abs(Typ_Pri(Spark(7)) - PForecast(1,Numeration));
        elseif Spark(7) == 0
            PForecast_Rev(1,Numeration) = 0;
        end
    end
    DeltaY3 = DeltaY4;
end

if Spark(7) ~= 0
    Hunt_Price = Typ_Pri(Spark(7));
else
    Hunt_Price = 0;
end

if (Numeration == 0) || (Numeration == 1) || (Numeration == 2) % The # encounters is not enough.
    First_Forecast_CT = {'Infinity'};
    First_Forecast_DT = -1;
    Second_Forecast_CT = {'Infinity'};
    Second_Forecast_DT = -1;
    Third_Forecast_CT = {'Infinity'};
    Third_Forecast_DT = -1;
    
    First_PForecast = -1;
    Second_PForecast = -1;
    Third_PForecast = -1;
    
    XNUD_Forecast = {'NO'};
    Hunt_Forecast = {'NO'};
    SellBuy = {'No'};
    
    MaxTyp_PriForecast = -1;
    MinTyp_PriForecast = -1;
    x_maxF = -1;
    x_minF = -1;
    
%     Max_PercentF = -1;
%     Deviation_F = -1;
    
    TypePosition_F = -1;
    TypePos_MaxPer_Dev_F = -1;
    TypePos_F = {'NO'};
    
    FinalPrice_F = -1;
    PiP = -1;
    NetProfit = -1;
    PerPosition = -1;
    
    TimeOpen = -1;
    PriceOpen = -1;
    TimeClose = -1;
    PriceClose = -1;
    
    New_lineP = -1;
    New_lineP_D = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1; Cut_Number=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
    
elseif YForecast(1,1) > MM % Because there is no real data.
    First_Forecast_CT = {'Infinity'};
    First_Forecast_DT = -1;
    Second_Forecast_CT = {'Infinity'};
    Second_Forecast_DT = -1;
    Third_Forecast_CT = {'Infinity'};
    Third_Forecast_DT = -1;
    
    First_PForecast = -1;
    Second_PForecast = -1;
    Third_PForecast = -1;
    
    XNUD_Forecast = {'NO'};
    Hunt_Forecast = {'NO'};
    SellBuy = {'No'};
    
    MaxTyp_PriForecast = -1;
    MinTyp_PriForecast = -1;
    x_maxF = -1;
    x_minF = -1;
    
%     Max_PercentF = -1;
%     Deviation_F = -1;
    
    TypePosition_F = -1;
    TypePos_MaxPer_Dev_F = -1;
    TypePos_F = {'NO'};
    
    FinalPrice_F = -1;
    PiP = -1;
    NetProfit = -1;
    PerPosition = -1;
    
    TimeOpen = -1;
    PriceOpen = -1;
    TimeClose = -1;
    PriceClose = -1;
    
    New_lineP = -1;
    New_lineP_D = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1; Cut_Number=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
    
elseif Spark(7) == 0 % The extreme point is not clear.
    First_Forecast_CT = {'Infinity'};
    First_Forecast_DT = -1;
    Second_Forecast_CT = {'Infinity'};
    Second_Forecast_DT = -1;
    Third_Forecast_CT = {'Infinity'};
    Third_Forecast_DT = -1;
    
    First_PForecast = -1;
    Second_PForecast = -1;
    Third_PForecast = -1;
    
    XNUD_Forecast = {'NO'};
    Hunt_Forecast = {'NO'};
    SellBuy = {'No'};
    
    MaxTyp_PriForecast = -1;
    MinTyp_PriForecast = -1;
    x_maxF = -1;
    x_minF = -1;
    
%     Max_PercentF = -1;
%     Deviation_F = -1;
    
    TypePosition_F = -1;
    TypePos_MaxPer_Dev_F = -1;
    TypePos_F = {'NO'};
    
    FinalPrice_F = -1;
    PiP = -1;
    NetProfit = -1;
    PerPosition = -1;
    
    TimeOpen = -1;
    PriceOpen = -1;
    TimeClose = -1;
    PriceClose = -1;
    
    New_lineP = -1;
    New_lineP_D = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1; Cut_Number=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
    
elseif YForecast(1,1) <= Spark(7) % The first encounter is before the first spark.
    First_Forecast_CT = {'Infinity'};
    First_Forecast_DT = -1;
    Second_Forecast_CT = {'Infinity'};
    Second_Forecast_DT = -1;
    Third_Forecast_CT = {'Infinity'};
    Third_Forecast_DT = -1;
    
    First_PForecast = -1;
    Second_PForecast = -1;
    Third_PForecast = -1;
    
    XNUD_Forecast = {'NO'};
    Hunt_Forecast = {'NO'};
    SellBuy = {'No'};
    
    MaxTyp_PriForecast = -1;
    MinTyp_PriForecast = -1;
    x_maxF = -1;
    x_minF = -1;
    
%     Max_PercentF = -1;
%     Deviation_F = -1;
    
    TypePosition_F = -1;
    TypePos_MaxPer_Dev_F = -1;
    TypePos_F = {'NO'};
    
    FinalPrice_F = -1;
    PiP = -1;
    NetProfit = -1;
    PerPosition = -1;
    
    TimeOpen = -1;
    PriceOpen = -1;
    TimeClose = -1;
    PriceClose = -1;
    
    New_lineP = -1;
    New_lineP_D = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1; Cut_Number=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
    
elseif (YForecast(1,1) - Spark(7)) < 1 % There is no data between first encounter and spark timing.
    First_Forecast_CT = {'Infinity'};
    First_Forecast_DT = -1;
    Second_Forecast_CT = {'Infinity'};
    Second_Forecast_DT = -1;
    Third_Forecast_CT = {'Infinity'};
    Third_Forecast_DT = -1;
    
    First_PForecast = -1;
    Second_PForecast = -1;
    Third_PForecast = -1;
    
    XNUD_Forecast = {'NO'};
    Hunt_Forecast = {'NO'};
    SellBuy = {'No'};
    
    MaxTyp_PriForecast = -1;
    MinTyp_PriForecast = -1;
    x_maxF = -1;
    x_minF = -1;
    
%     Max_PercentF = -1;
%     Deviation_F = -1;
    
    TypePosition_F = -1;
    TypePos_MaxPer_Dev_F = -1;
    TypePos_F = {'NO'};
    
    FinalPrice_F = -1;
    PiP = -1;
    NetProfit = -1;
    PerPosition = -1;
    
    TimeOpen = -1;
    PriceOpen = -1;
    TimeClose = -1;
    PriceClose = -1;
    
    New_lineP = -1;
    New_lineP_D = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1; Cut_Number=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
    
elseif (YForecast(1,1) - tt_T(7) + 1) <= (t7_t1(1,1) / 2)
    First_Forecast_CT = NY2(1,1);
    First_Forecast_DT = YForecast(1,1);
    Second_Forecast_CT = NY2(1,2);
    Second_Forecast_DT = YForecast(1,2);
    Third_Forecast_CT = NY2(1,3);
    Third_Forecast_DT = YForecast(1,3);
    
    First_PForecast = PForecast(1,1);
    Second_PForecast = PForecast(1,2);
    Third_PForecast = PForecast(1,3);
    
    if YForecast(1,1) <= MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : YForecast(1,1) );
        MaxTyp_PriForecast = max(MTyp_Pri); % Find max
        MinTyp_PriForecast = min(MTyp_Pri); % Find min
        x_maxF = find(MTyp_Pri == MaxTyp_PriForecast) + Spark(7);
        x_minF = find(MTyp_Pri == MinTyp_PriForecast) + Spark(7);
    elseif YForecast(1,1) > MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : MM );
        MaxTyp_PriForecast = max(MTyp_Pri); % Find max
        MinTyp_PriForecast = min(MTyp_Pri); % Find min
        x_maxF = find(MTyp_Pri == MaxTyp_PriForecast) + Spark(7);
        x_minF = find(MTyp_Pri == MinTyp_PriForecast) + Spark(7);
    end
    
    if (b_T(6) < b_T(7)) && (b_T(7) < PForecast(1,1)) % b_T(7) is max
        XNUD_Forecast(1,1) = {'XU'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) > PForecast(1,1)) % b_T(7) is max
        XNUD_Forecast(1,1) = {'XD'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) < PForecast(1,1)) % b_T(7) is min
        XNUD_Forecast(1,1) = {'NU'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) > PForecast(1,1)) % b_T(7) is min
        XNUD_Forecast(1,1) = {'ND'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) == PForecast(1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Forecast(1,1) = {'XE'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) == PForecast(1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Forecast(1,1) = {'NE'};
    end
    
    TimeOpen = Spark(7);
    PriceOpen = Typ_Pri(Spark(7));
    
    if Typ_Pri(Spark(7)) < First_PForecast % Hunt Up
        Hunt_Forecast = {'Hunt_U'};
        SellBuy = {'Buy'};
        TypePosition_F = -1; % No mode
        New_lineP = First_PForecast + 4 * 10^(-(Digit-1));
        New_lineP_D = Typ_Pri(Spark(7)) - abs(First_PForecast - Typ_Pri(Spark(7))) + 4 * 10^(-(Digit-1));
        
        for i = Spark(7) : YForecast(1,1) - 1
            if (((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                    || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))) ...
                    && (((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                    || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D)))
                TypePosition_F = 2; % Deviation
                NSZO = NSZO + 1;
                SequenceZeroOne(NSZO) = 0;
                TimeClose = i;
                PriceClose = New_lineP_D;
                break
            elseif ((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                    || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))
                TypePosition_F = 1; % Max_Percent
                NSZO = NSZO + 1;
                SequenceZeroOne(NSZO) = 1;
                TimeClose = i;
                PriceClose = New_lineP;
                break
            elseif ((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                    || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D))
                TypePosition_F = 2; % Deviation
                NSZO = NSZO + 1;
                SequenceZeroOne(NSZO) = 0;
                TimeClose = i;
                PriceClose = New_lineP_D;
                break
            end
        end
        i = i + 1;
        if TypePosition_F == -1 % No mode
            if Typ_Pri(YForecast(1,1)) >= Typ_Pri(Spark(7)) + 4 * 10^(-(Digit-1))
                TypePosition_F = 3; % No mode & in profit
                TimeClose = YForecast(1,1);
                PriceClose = Typ_Pri(YForecast(1,1));
            else
                while (Typ_Pri(i) < Typ_Pri(Spark(7)) + 4 * 10^(-(Digit-1))) ...
                        && (Typ_Pri(i) > New_lineP_D) && (i < MM)
                    i = i + 1;
                end
                if i > MM
                    TypePosition_F = 100; % No mode & ignore
                    TimeClose = -1;
                    PriceClose = -1;
                elseif Typ_Pri(i) >= Typ_Pri(Spark(7)) + 4 * 10^(-(Digit-1))
                    TypePosition_F = 32; % No mode & Risk Free
                    TimeClose = i;
                    PriceClose = Typ_Pri(Spark(7)) + 4 * 10^(-(Digit-1));
                elseif Typ_Pri(i) <= New_lineP_D
                    TypePosition_F = 33; % No mode & Deviation
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 0;
                    TimeClose = i;
                    PriceClose = New_lineP_D;
                else
                    TypePosition_F = 100; % No mode & ignore
                    TimeClose = -1;
                    PriceClose = -1;
                end
            end
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark(7)) > First_PForecast % Hunt Down
        Hunt_Forecast = {'Hunt_D'};
        SellBuy = {'Sell'};
        TypePosition_F = -1; % No mode
        New_lineP = First_PForecast - 4 * 10^(-(Digit-1));
        New_lineP_D = Typ_Pri(Spark(7)) + abs(First_PForecast - Typ_Pri(Spark(7))) - 4 * 10^(-(Digit-1));
        
        for i = Spark(7) : YForecast(1,1) - 1
            if (((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                    || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))) ...
                    && (((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                    || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D)))
                TypePosition_F = 2; % Deviation
                NSZO = NSZO + 1;
                SequenceZeroOne(NSZO) = 0;
                TimeClose = i;
                PriceClose = New_lineP_D;
                break
            elseif ((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                    || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))
                TypePosition_F = 1; % Max_Percent
                NSZO = NSZO + 1;
                SequenceZeroOne(NSZO) = 1;
                TimeClose = i;
                PriceClose = New_lineP;
                break
            elseif ((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                    || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D))
                TypePosition_F = 2; % Deviation
                NSZO = NSZO + 1;
                SequenceZeroOne(NSZO) = 0;
                TimeClose = i;
                PriceClose = New_lineP_D;
                break
            end
        end
        i = i + 1;
        if TypePosition_F == -1 % No mode
            if Typ_Pri(YForecast(1,1)) <= Typ_Pri(Spark(7)) - 4 * 10^(-(Digit-1))
                TypePosition_F = 3; % No mode & in profit
                TimeClose = YForecast(1,1);
                PriceClose = Typ_Pri(YForecast(1,1));
            else
                while (Typ_Pri(i) > Typ_Pri(Spark(7)) - 4 * 10^(-(Digit-1))) ...
                        && (Typ_Pri(i) < New_lineP_D) && (i < MM)
                    i = i + 1;
                end
                if i > MM
                    TypePosition_F = 100; % No mode & ignore
                    TimeClose = -1;
                    PriceClose = -1;
                elseif Typ_Pri(i) <= Typ_Pri(Spark(7)) - 4 * 10^(-(Digit-1))
                    TypePosition_F = 32; % No mode & Risk Free
                    TimeClose = i;
                    PriceClose = Typ_Pri(Spark(7)) - 4 * 10^(-(Digit-1));
                elseif Typ_Pri(i) >= New_lineP_D
                    TypePosition_F = 33; % No mode & Deviation
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 0;
                    TimeClose = i;
                    PriceClose = New_lineP_D;
                else
                    TypePosition_F = 100; % No mode & ignore
                    TimeClose = -1;
                    PriceClose = -1;
                end
            end
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark(7)) == First_PForecast % Hunt Equal
        Hunt_Forecast = {'Hunt_E'};
        SellBuy = {'No'};
%         Max_PercentF = -1;
%         Deviation_F = -1;
        
        TypePosition_F = -1;
        TypePos_MaxPer_Dev_F = -1;
        TypePos_F = {'NO'};
        
        FinalPrice_F = Typ_Pri(YForecast(1,1));
        PiP = -1;
        NetProfit = -1;
        PerPosition = -1;
        
        TimeClose = -1;
        PriceClose = -1;
        
        New_lineP = -1;
        New_lineP_D = -1;
    end
    %--------------------------  Select cut 1,2,3  ------------------------
    x1 = Sol_Total(1,1) : YForecast(1,3); % First Curves
    x_curve2 = Sol_Total(1,3) : YForecast(1,3); % Second Curves
    
    r1 = a1_Previous * x1 + a0_Previous; % First Curves
    y10 = r1 + Sol_Total(1,15) * sin(Sol_Total(1,16) * x1 + Sol_Total(1,17)) ...
        + Sol_Total(1,18) * sin(Sol_Total(1,19) * x1 + Sol_Total(1,20)); % First Curves
    
    r2 = a1 * x_curve2 + a0; % Second Curves
    y11 = r2 + Sol_Total(1,24) * sin(Sol_Total(1,25) * x_curve2 + Sol_Total(1,26)) ...
        + Sol_Total(1,27) * sin(Sol_Total(1,28) * x_curve2 + Sol_Total(1,29)); % Second Curves
    
%     figure()
%     plot(x1,y10,'b.')
    S_E_1 = Select_Extreme(y10, Sol_Total(1,1), YForecast(1,3)); % First Curves
%     hold on
%     plot(x_curve2,y11,'g.')
    S_E_2 = Select_Extreme(y11, Sol_Total(1,3), YForecast(1,3)); % Second Curves
    
    if (S_E_1(1,1) ~= 0) && (S_E_2(1,1) ~= 0)
        % First Curves
        Ext1_Cut_3_time = S_E_1(size(S_E_1,1),1); % Cut 3
        %Ext1_Cut_3_Price = S_E_1(size(S_E_1,1),2);
        Ext1_Cut_3_MinMax = S_E_1(size(S_E_1,1),3);
        for i = size(S_E_1,1) : -1 : 1 % Cut 2
            if YForecast(1,2) > S_E_1(i,1)
                Ext1_Cut_2_time = S_E_1(i,1);
                %Ext1_Cut_2_Price = S_E_1(i,2);
                Ext1_Cut_2_MinMax = S_E_1(i,3);
                break
            end
        end
        for i = size(S_E_1,1) : -1 : 1 % Cut 1
            if YForecast(1,1) > S_E_1(i,1)
                Ext1_Cut_1_time = S_E_1(i,1);
                %Ext1_Cut_1_Price = S_E_1(i,2);
                Ext1_Cut_1_MinMax = S_E_1(i,3);
                break
            end
        end
        % Second Curves
        Ext2_Cut_3_time = S_E_2(size(S_E_2,1),1); % Cut 3
        %Ext2_Cut_3_Price = S_E_2(size(S_E_2,1),2);
        Ext2_Cut_3_MinMax = S_E_2(size(S_E_2,1),3);
        for i = size(S_E_2,1) : -1 : 1 % Cut 2
            if YForecast(1,2) > S_E_2(i,1)
                Ext2_Cut_2_time = S_E_2(i,1);
                %Ext2_Cut_2_Price = S_E_2(i,2);
                Ext2_Cut_2_MinMax = S_E_2(i,3);
                break
            end
        end
        for i = size(S_E_2,1) : -1 : 1 % Cut 1
            if YForecast(1,1) > S_E_2(i,1)
                Ext2_Cut_1_time = S_E_2(i,1);
                %Ext2_Cut_1_Price = S_E_2(i,2);
                Ext2_Cut_1_MinMax = S_E_2(i,3);
                break
            end
        end
        ii3 = ii3 + 1;
        Trend_Curves_Cut_1(ii3,1) = Result_Trend_Curves(Ext1_Cut_1_time, Ext2_Cut_1_time, ...
            Ext1_Cut_1_MinMax, Ext2_Cut_1_MinMax, YForecast(1,1), S_E_1, S_E_2); % Cut 1
        
        Trend_Curves_Cut_2(ii3,1) = Result_Trend_Curves(Ext1_Cut_2_time, Ext2_Cut_2_time, ...
            Ext1_Cut_2_MinMax, Ext2_Cut_2_MinMax, YForecast(1,2), S_E_1, S_E_2); % Cut 2
        
        Trend_Curves_Cut_3(ii3,1) = Result_Trend_Curves(Ext1_Cut_3_time, Ext2_Cut_3_time, ...
            Ext1_Cut_3_MinMax, Ext2_Cut_3_MinMax, YForecast(1,3), S_E_1, S_E_2); % Cut 3
        
        %-----------------------  Three reward to risk  -------------------
        if (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==0)
            % RR=1, TP=100, SL=100
            XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
            TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
            PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
            Cut_Number = -1;
            %  Reverse
            XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
            TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
            PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
            
            % RR=1, TP=50, SL=50
            XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
            TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
            PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
            %  Reverse
            XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
            TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
            PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
            
            % RR=3, TP=150, SL=50
            XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
            TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
            PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
            %  Reverse
            XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
            TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
            PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==0) % Cut 1
            if abs(PForecast(1,1) - Typ_Pri(Spark(7))) >= Dist_HC
                Cut_Number = 1;
                % R/R=1, TP=100, SL=100
                [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                    PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                    NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                
                if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                    PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                    NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                
                if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                end
                
                % R/R=1, TP=50, SL=50
                [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                    PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                    NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                
                if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                    PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                    NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                
                if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                end
                
                % R/R=3, TP=150, SL=50
                [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                    PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                    NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                
                if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                    PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                    NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                
                if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                end
            else
                Cut_Number = -1;
                % R/R=1, TP=100, SL=100
                XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                %  Reverse
                XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                
                % R/R=1, TP=50, SL=50
                XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                %  Reverse
                XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                
                % R/R=3, TP=150, SL=50
                XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                %  Reverse
                XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==0) % Cut 2
            if abs(PForecast(1,2) - Typ_Pri(Spark(7))) >= Dist_HC
                Cut_Number = 2;
                % R/R=1, TP=100, SL=100
                [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                    PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                    NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                
                if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                    PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                    NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                
                if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                end
                
                % R/R=1, TP=50, SL=50
                [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                    PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                    NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                
                if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                    PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                    NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                
                if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                end
                
                % R/R=3, TP=150, SL=50
                [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                    PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                    NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                
                if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                    PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                    NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                
                if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                end
            else
                Cut_Number = -1;
                % R/R=1, TP=100, SL=100
                XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                %  Reverse
                XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                
                % R/R=1, TP=50, SL=50
                XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                %  Reverse
                XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                
                % R/R=3, TP=150, SL=50
                XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                %  Reverse
                XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==1) % Cut 3
            if abs(PForecast(1,3) - Typ_Pri(Spark(7))) >= Dist_HC
                Cut_Number = 3;
                % R/R=1, TP=100, SL=100
                [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                    PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                    NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                
                if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                    PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                    NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                
                if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                end
                
                % R/R=1, TP=50, SL=50
                [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                    PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                    NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                
                if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                    PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                    NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                
                if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                end
                
                % R/R=3, TP=150, SL=50
                [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                    PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                    NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                
                if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                    PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                    = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                    NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                
                if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                end
            else
                Cut_Number = -1;
                % R/R=1, TP=100, SL=100
                XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                %  Reverse
                XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                
                % R/R=1, TP=50, SL=50
                XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                %  Reverse
                XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                
                % R/R=3, TP=150, SL=50
                XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                %  Reverse
                XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==0)
            if abs(Typ_Pri(Spark(7)) - PForecast(1,1)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,2)) % Cut 1
                if abs(PForecast(1,1) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 1;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            elseif abs(Typ_Pri(Spark(7)) - PForecast(1,1)) > abs(Typ_Pri(Spark(7)) - PForecast(1,2)) % Cut 2
                if abs(PForecast(1,2) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 2;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==1)
            if abs(Typ_Pri(Spark(7)) - PForecast(1,1)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,3)) % Cut 1
                if abs(PForecast(1,1) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 1;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            elseif abs(Typ_Pri(Spark(7)) - PForecast(1,1)) > abs(Typ_Pri(Spark(7)) - PForecast(1,3)) % Cut 3
                if abs(PForecast(1,3) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 3;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==1)
            if abs(Typ_Pri(Spark(7)) - PForecast(1,2)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,3)) % Cut 2
                if abs(PForecast(1,2) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 2;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            elseif abs(Typ_Pri(Spark(7)) - PForecast(1,2)) > abs(Typ_Pri(Spark(7)) - PForecast(1,3)) % Cut 3
                if abs(PForecast(1,3) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 3;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % RR=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % RR=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==1)
            if ( abs(Typ_Pri(Spark(7)) - PForecast(1,1)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,2)) ) && ...
                    ( abs(Typ_Pri(Spark(7)) - PForecast(1,1)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,3)) ) % Cut 1
                if abs(PForecast(1,1) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 1;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast(1,1), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,1), PForecast_Rev(1,1), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            elseif ( abs(Typ_Pri(Spark(7)) - PForecast(1,2)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,1)) ) && ...
                    ( abs(Typ_Pri(Spark(7)) - PForecast(1,2)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,3)) ) % Cut 2
                if abs(PForecast(1,2) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 2;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast(1,2), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,2), PForecast_Rev(1,2), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            elseif ( abs(Typ_Pri(Spark(7)) - PForecast(1,3)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,1)) ) && ...
                    ( abs(Typ_Pri(Spark(7)) - PForecast(1,3)) <= abs(Typ_Pri(Spark(7)) - PForecast(1,2)) ) % Cut 3
                if abs(PForecast(1,3) - Typ_Pri(Spark(7))) >= Dist_HC
                    Cut_Number = 3;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
                        PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, NSZO_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1 == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
                        PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, NSZO_Rev_R1] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R1, b_T(6), b_T(7), Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1 == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
                        PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, NSZO_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2 == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
                        PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, NSZO_Rev_R2] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R2, b_T(6), b_T(7), Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2 == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
                        PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, NSZO_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast(1,3), Spark(7), ...
                        NSZO_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3 == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
                        PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3, NSZO_Rev_R3] ...
                        = Select_Profit_1_S(Typ_Pri, YForecast(1,3), PForecast_Rev(1,3), Spark(7), ...
                        NSZO_Rev_R3, b_T(6), b_T(7), Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3 == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
                    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
                    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
                    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
                    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
                    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
                    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
                    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
                    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
                    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
                    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
                    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
                    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
                end
            end
        end
    else
        Cut_Number = -1;
        % R/R=1, TP=100, SL=100
        XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
        TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
        PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
        %  Reverse
        XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
        TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
        PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
        
        % R/R=1, TP=50, SL=50
        XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
        TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
        PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
        %  Reverse
        XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
        TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
        PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
        
        % R/R=3, TP=150, SL=50
        XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
        TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
        PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
        %  Reverse
        XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
        TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
        PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
    end
else
    First_Forecast_CT = {'Infinity'};
    First_Forecast_DT = -1;
    Second_Forecast_CT = {'Infinity'};
    Second_Forecast_DT = -1;
    Third_Forecast_CT = {'Infinity'};
    Third_Forecast_DT = -1;
    
    First_PForecast = -1;
    Second_PForecast = -1;
    Third_PForecast = -1;
    
    XNUD_Forecast = {'NO'};
    Hunt_Forecast = {'NO'};
    SellBuy = {'No'};
    
    MaxTyp_PriForecast = -1;
    MinTyp_PriForecast = -1;
    x_maxF = -1;
    x_minF = -1;
    
%     Max_PercentF = -1;
%     Deviation_F = -1;
    
    TypePosition_F = -1;
    TypePos_MaxPer_Dev_F = -1;
    TypePos_F = {'NO'};
    
    FinalPrice_F = -1;
    PiP = -1;
    NetProfit = -1;
    PerPosition = -1;
    
    TimeOpen = -1;
    PriceOpen = -1;
    TimeClose = -1;
    PriceClose = -1;
    
    New_lineP = -1;
    New_lineP_D = -1;
    
    Cut_Number=-1;
    % R/R=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1={'No'}; Hunt_Forecast_Sel_R1={'No'}; SellBuy_Sel_R1={'No'}; 
    TimeOpen_Sel_R1=-1; PriceOpen_Sel_R1=-1; TimeClose_Sel_R1=-1; 
    PriceClose_Sel_R1=-1; TypePosition_Sel_R1=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1={'No'}; Hunt_Forecast_Rev_R1={'No'}; SellBuy_Rev_R1={'No'};
    TimeOpen_Rev_R1=-1; PriceOpen_Rev_R1=-1; TimeClose_Rev_R1=-1; 
    PriceClose_Rev_R1=-1; TypePosition_Rev_R1=-1;
    
    % R/R=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2={'No'}; Hunt_Forecast_Sel_R2={'No'}; SellBuy_Sel_R2={'No'}; 
    TimeOpen_Sel_R2=-1; PriceOpen_Sel_R2=-1; TimeClose_Sel_R2=-1; 
    PriceClose_Sel_R2=-1; TypePosition_Sel_R2=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2={'No'}; Hunt_Forecast_Rev_R2={'No'}; SellBuy_Rev_R2={'No'};
    TimeOpen_Rev_R2=-1; PriceOpen_Rev_R2=-1; TimeClose_Rev_R2=-1; 
    PriceClose_Rev_R2=-1; TypePosition_Rev_R2=-1;
    
    % R/R=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3={'No'}; Hunt_Forecast_Sel_R3={'No'}; SellBuy_Sel_R3={'No'}; 
    TimeOpen_Sel_R3=-1; PriceOpen_Sel_R3=-1; TimeClose_Sel_R3=-1; 
    PriceClose_Sel_R3=-1; TypePosition_Sel_R3=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3={'No'}; Hunt_Forecast_Rev_R3={'No'}; SellBuy_Rev_R3={'No'};
    TimeOpen_Rev_R3=-1; PriceOpen_Rev_R3=-1; TimeClose_Rev_R3=-1; 
    PriceClose_Rev_R3=-1; TypePosition_Rev_R3=-1;
end

% Center of gravity
if Numeration == 3
    Gravity_tF = ( YForecast(1,1) + YForecast(1,2) + YForecast(1,3) ) / 3;
    Gravity_PF = ( PForecast(1,1) + PForecast(1,2) + PForecast(1,3) ) / 3;
else
    Gravity_tF = -1;
    Gravity_PF = -1;
end

if (Numeration ==0) || (Numeration ==1) || (Numeration ==2) % The # encounters is not enough.
    MaxTyp_PriceFG = -1;
    MinTyp_PriceFG = -1;
    x_maxFG = -1; % Additional
    x_minFG = -1; % Additional
    
    XNUD_FG(1,1) = {'NO'};
    Hunt_FG(1,1) = {'NO'};
    
    Max_PercentFG = -1;
    Deviation_FG = -1;
elseif Gravity_tF > MM % Because there is no real data.
    MaxTyp_PriceFG = -1;
    MinTyp_PriceFG = -1;
    x_maxFG = -1; % Additional
    x_minFG = -1; % Additional
    
    XNUD_FG(1,1) = {'NO'};
    Hunt_FG(1,1) = {'NO'};
    
    Max_PercentFG = -1;
    Deviation_FG = -1;
elseif Spark(7) == 0 % The extreme point is not clear.
    MaxTyp_PriceFG = -1;
    MinTyp_PriceFG = -1;
    x_maxFG = -1; % Additional
    x_minFG = -1; % Additional
    
    XNUD_FG(1,1) = {'NO'};
    Hunt_FG(1,1) = {'NO'};
    
    Max_PercentFG = -1;
    Deviation_FG = -1;
elseif Gravity_tF <= Spark(7) % The first encounter is before the first spark.
    MaxTyp_PriceFG = -1;
    MinTyp_PriceFG = -1;
    x_maxFG = -1; % Additional
    x_minFG = -1; % Additional
    
    XNUD_FG(1,1) = {'NO'};
    Hunt_FG(1,1) = {'NO'};
    
    Max_PercentFG = -1;
    Deviation_FG = -1;
elseif (Gravity_tF - Spark(7)) < 1 % There is no data between first encounter and spark timing.
    MaxTyp_PriceFG = -1;
    MinTyp_PriceFG = -1;
    x_maxFG = -1; % Additional
    x_minFG = -1; % Additional
    
    XNUD_FG(1,1) = {'NO'};
    Hunt_FG(1,1) = {'NO'};
    
    Max_PercentFG = -1;
    Deviation_FG = -1;
elseif (Gravity_tF - tt_T(7) + 1) <= (t7_t1(1,1) / 2)
    if Gravity_tF <= MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : Gravity_tF );
        MaxTyp_PriceFG = max(MTyp_Pri); % Find max
        MinTyp_PriceFG = min(MTyp_Pri); % Find min
        x_maxFG = find(MTyp_Pri == MaxTyp_PriceFG); % Additional
        x_minFG = find(MTyp_Pri == MinTyp_PriceFG); % Additional
    elseif Gravity_tF > MM
        MTyp_Pri = Typ_Pri( Spark(7)+1 : MM );
        MaxTyp_PriceFG = max(MTyp_Pri);
        MinTyp_PriceFG = min(MTyp_Pri);
        x_maxFG = find(MTyp_Pri == MaxTyp_PriceFG); % Additional
        x_minFG = find(MTyp_Pri == MinTyp_PriceFG); % Additional
    end
    
    if (b_T(6) < b_T(7)) && (b_T(7) < Gravity_PF) % b_T(7) is max
        XNUD_FG(1,1) = {'XU'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) > Gravity_PF) % b_T(7) is max
        XNUD_FG(1,1) = {'XD'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) < Gravity_PF) % b_T(7) is min
        XNUD_FG(1,1) = {'NU'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) > Gravity_PF) % b_T(7) is min
        XNUD_FG(1,1) = {'ND'};
    elseif (b_T(6) < b_T(7)) && (b_T(7) == Gravity_PF) % b_T(7) is equal to PForecast(1,1).
        XNUD_FG(1,1) = {'XE'};
    elseif (b_T(6) > b_T(7)) && (b_T(7) == Gravity_PF) % b_T(7) is equal to PForecast(1,1).
        XNUD_FG(1,1) = {'NE'};
    end
    
    if Typ_Pri(Spark(7)) < Gravity_PF % Hunt Up
        Hunt_FG = {'Hunt_U'};
        if MaxTyp_PriceFG <= Typ_Pri(Spark(7))
            Max_PercentFG(1,1) = 0;
        elseif ( Typ_Pri(Spark(7)) < MaxTyp_PriceFG ) ...
                && ( MaxTyp_PriceFG < Gravity_PF )
            Max_PercentFG(1,1) = abs( MaxTyp_PriceFG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PF - Typ_Pri(Spark(7)) );
        elseif MaxTyp_PriceFG >= Gravity_PF
            Max_PercentFG(1,1) = 100;
        end
        if MinTyp_PriceFG >= Typ_Pri(Spark(7))
            Deviation_FG(1,1) = 0;
        elseif MinTyp_PriceFG < Typ_Pri(Spark(7))
            Deviation_FG(1,1) = abs( MinTyp_PriceFG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PF - Typ_Pri(Spark(7)) );
        end
    elseif Typ_Pri(Spark(7)) > Gravity_PF % Hunt Down
        Hunt_FG = {'Hunt_D'};
        if MinTyp_PriceFG >= Typ_Pri(Spark(7))
            Max_PercentFG(1,1) = 0;
        elseif ( Typ_Pri(Spark(7)) > MinTyp_PriceFG ) ...
                && ( MinTyp_PriceFG(1,1) > Gravity_PF )
            Max_PercentFG(1,1) = abs( MinTyp_PriceFG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PF - Typ_Pri(Spark(7)) );
        elseif MinTyp_PriceFG <= Gravity_PF
            Max_PercentFG(1,1) = 100;
        end
        if MaxTyp_PriceFG <= Typ_Pri(Spark(7))
            Deviation_FG(1,1) = 0;
        elseif MaxTyp_PriceFG > Typ_Pri(Spark(7))
            Deviation_FG(1,1) = abs( MaxTyp_PriceFG - Typ_Pri(Spark(7)) ) * 100 / ...
                abs( Gravity_PF - Typ_Pri(Spark(7)) );
        end
    elseif Typ_Pri(Spark(7)) == Gravity_PF % Hunt Equal
        Hunt_FG = {'Hunt_E'};
        Max_PercentFG = -1;
        Deviation_FG = -1;
    end
else
    MaxTyp_PriceFG = -1;
    MinTyp_PriceFG = -1;
    x_maxFG = -1; % Additional
    x_minFG = -1; % Additional
    
    XNUD_FG(1,1) = {'NO'};
    Hunt_FG(1,1) = {'NO'};
    
    Max_PercentFG = -1;
    Deviation_FG = -1;
end
%---------------------------------  No noise  -----------------------------
% %x1 = 1 : Third_encounter_DT(1);
% if NumerationE == 3
% x_curve1 = Sol_Total(1,1) : NDeltaY(1,3);
% x_curve2 = Sol_Total(1,3) : NDeltaY(1,3);
% 
% r1 = a1_Previous * x_curve1 + a0_Previous;
% y6 = r1 + Sol_Total(1,15) * sin(Sol_Total(1,16) * x_curve1 + Sol_Total(1,17));
% 
% r2 = a1 * x_curve2 + a0;
% y7 = r2 + Sol_Total(1,24) * sin(Sol_Total(1,25) * x_curve2 + Sol_Total(1,26));
% if NDeltaY(1,3) > MM
%     figure()
%     p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
%     hold on
%     p2 = plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
%     hold on
% elseif Second_encounter_DT(1,1) > MM
%     figure()
%     p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
%     hold on
%     p2 = plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
%     hold on
% elseif First_encounter_DT(1,1) > MM
%     figure()
%     p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
%     hold on
%     p2 = plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
%     hold on
% else
%     figure()
%     %p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
%     p1 = plot(x_curve1,Typ_Pri(Sol_Total(1,1):NDeltaY(1,3)),'k.');
%     hold on
%     bi = 7;
%     while ( bi < length(tt_T) ) && ( tt_T(bi) < NDeltaY(1,3) )
%         bi = bi + 1;
%     end
%     p2 = plot(tt_T(1:bi),b_T(1:bi),'mo','linewidth',2);
%     hold on
% end
% p3 = plot(x_curve1,y6,'b');
% hold on
% p4 = plot(x_curve2,y7,'g');
% title('Systems 1 and 2: r+a_1*sin(w_1t+\theta_1)')
% %legend('Real Data','Six Points','First System','Second System')
% %axis tight
% Y_1 = ylim;
% 
% for i = 1 : 3 % Additional
%     hold on
%     X_1 = [NDeltaY(1,i), NDeltaY(1,i)];
%     p5 = line(X_1,Y_1,'Color',[.5 .3 1.0]);%[.6 .1 .9]
%     hold on
%     p6 = plot(NDeltaY(1,i),PEncounter(1,i),'ch','linewidth',2); % Additional
% end
% hold on
% p7 = plot(x_maxE,MaxTyp_PriEncounter,'gp','linewidth',1.5); % Additional
% hold on
% plot(x_minE,MinTyp_PriEncounter,'gp','linewidth',1.5); % Additional
% hold on
% X_1 = [Spark(7), Spark(7)];
% p8 = line(X_1,Y_1,'Color',[.3 .9 .9]); % Additional
% hold on
% p9 = plot(Gravity_tE,Gravity_PE,'r>','linewidth',1.5); %p9.Marker = '>';
% hold on
% X_1 = [Gravity_tE, Gravity_tE];
% p10 = line(X_1,Y_1,'Color',[1.0 .2 1.0]);%[.9 .1 .9]
% p10.LineStyle = '--';
% hold on
% p11 = plot(x_maxEG+Spark(7),MaxTyp_PriceEG,'g^','linewidth',1.5); % Additional
% hold on
% plot(x_minEG+Spark(7),MinTyp_PriceEG,'g^','linewidth',1.5); % Additional
% legend([p1,p2,p3,p4,p5,p6,p7(1,1),p8,p9,p10,p11(1,1)],'Real Data','Six Points','First System',...
%     'Second System','Encounter of Two System','Encounter of Two System',...
%     'Min & Max Real Data','Hunt','Center of gravity','Center of gravity',...
%     'Min & Max Real Data (gravity)')
% hold on
% plot(Spark(7),Typ_Pri(Spark(7)),'rh','linewidth',2) % Additional
% 
% end

%---------------------------  Life Time (no noise)  -----------------------
%x_curve1 = Sol_Total(1,11) : NDeltaY(1,3);
x_curve2 = Spark(7) + 1 : (tt_T(7) + 2 * t7_t1(1,1));

y16 =@(t_1) a1_Previous * t_1 + a0_Previous ...
    + Sol_Total(1,15) * sin(Sol_Total(1,16) * t_1 + Sol_Total(1,17));

y17 =@(t_1) a1 * t_1 + a0 ...
    + Sol_Total(1,24) * sin(Sol_Total(1,25) * t_1 + Sol_Total(1,26));

dis1(1) = abs( y16(Sol_Total(1,1)) - Sol_Total(1,2) );
dis1(2) = abs( y16(Sol_Total(1,3)) - Sol_Total(1,4) );
dis1(3) = abs( y16(Sol_Total(1,5)) - Sol_Total(1,6) );
dis1(4) = abs( y16(Sol_Total(1,7)) - Sol_Total(1,8) );
dis1(5) = abs( y16(Sol_Total(1,9)) - Sol_Total(1,10) );
dis1(6) = abs( y16(Sol_Total(1,11)) - Sol_Total(1,12) );
max_dis1 = max(dis1);
average1 =  ( dis1(1) + dis1(2) + dis1(3) + dis1(4) + dis1(5) + dis1(6) ) / 6;

dis2(1) = abs( y17(Sol_Total(1,3)) - Sol_Total(1,4) );
dis2(2) = abs( y17(Sol_Total(1,5)) - Sol_Total(1,6) );
dis2(3) = abs( y17(Sol_Total(1,7)) - Sol_Total(1,8) );
dis2(4) = abs( y17(Sol_Total(1,9)) - Sol_Total(1,10) );
dis2(5) = abs( y17(Sol_Total(1,11)) - Sol_Total(1,12) );
dis2(6) = abs( y17(Sol_Total(1,13)) - Sol_Total(1,14) );
max_dis2 = max(dis2);
average2 =  ( dis2(1) + dis2(2) + dis2(3) + dis2(4) + dis2(5) + dis2(6) ) / 6;

y20 = y16(x_curve2) + (max_dis1 + average1) / 2;
y21 = y16(x_curve2) - (max_dis1 + average1) / 2;

y22 = y17(x_curve2) + (max_dis2 + average2) / 2;
y23 = y17(x_curve2) - (max_dis2 + average2) / 2;

LifeTime1 = -1; LifeTime2 = -1;
if (tt_T(7) + 2 * t7_t1(1,1)) <= MM
    for i = Spark(7) + 1 : (tt_T(7) + 2 * t7_t1(1,1)) - 1
        if ( Typ_Pri(i+1) >= y20(i-Spark(7)+1) ) && ( Typ_Pri(i) < y20(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y20(i-Spark(7)+1) ) && ( Typ_Pri(i) > y20(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y21(i-Spark(7)+1) ) && ( Typ_Pri(i) > y21(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y21(i-Spark(7)+1) ) && ( Typ_Pri(i) < y21(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        end
    end
else
    for i = Spark(7) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y20(i-Spark(7)+1) ) && ( Typ_Pri(i) < y20(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y20(i-Spark(7)+1) ) && ( Typ_Pri(i) > y20(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y21(i-Spark(7)+1) ) && ( Typ_Pri(i) > y21(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y21(i-Spark(7)+1) ) && ( Typ_Pri(i) < y21(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        end
    end
end
if (tt_T(7) + 2 * t7_t1(1,1)) <= MM
    for i = Spark(7) + 1 : (tt_T(7) + 2 * t7_t1(1,1)) - 1
        if ( Typ_Pri(i+1) >= y22(i-Spark(7)+1) ) && ( Typ_Pri(i) < y22(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y22(i-Spark(7)+1) ) && ( Typ_Pri(i) > y22(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y23(i-Spark(7)+1) ) && ( Typ_Pri(i) > y23(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y23(i-Spark(7)+1) ) && ( Typ_Pri(i) < y23(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        end
    end
else
    for i = Spark(7) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y22(i-Spark(7)+1) ) && ( Typ_Pri(i) < y22(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y22(i-Spark(7)+1) ) && ( Typ_Pri(i) > y22(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y23(i-Spark(7)+1) ) && ( Typ_Pri(i) > y23(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y23(i-Spark(7)+1) ) && ( Typ_Pri(i) < y23(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        end
    end
end
if (LifeTime1 == -1) && (LifeTime2 == -1)
    LifeTimeA = -1;
elseif (LifeTime1 ~= -1) && (LifeTime2 == -1)
    LifeTimeA = LifeTime1; % Life time for no nise
elseif (LifeTime1 == -1) && (LifeTime2 ~= -1)
    LifeTimeA = LifeTime2; % Life time for no nise
elseif (LifeTime1 ~= -1) && (LifeTime2 ~= -1)
    LifeTimeA = max(LifeTime1,LifeTime2); % Life time for no nise
else
    LifeTimeA = -1;
end

% hold on
% plot(x_curve2,y20,'b-.',x_curve2,y21,'b-.');
% hold on
% plot(x_curve2,y22,'g-.',x_curve2,y23,'g-.');
% hold on
% if LifeTime1 ~= -1
%     p12 = line(LifeTime1,Typ_Pri(LifeTime1),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p12.Marker = 's';
% end
% hold on
% if LifeTime2 ~= -1
%     p13 = line(LifeTime2,Typ_Pri(LifeTime2),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p13.Marker = 's';
% end

%--------------------------------  With noise  ----------------------------
x1 = 1 : Third_Forecast_DT(1);
if (Numeration == 3) && ( (TypePosition_Sel_R2 == 1) || ...
        (TypePosition_Sel_R2 == 2) || (TypePosition_Sel_R2 == 33) )
x_curve1 = Sol_Total(1,1) : YForecast(1,3);
x_curve2 = Sol_Total(1,3) : YForecast(1,3);

r1 = a1_Previous * x_curve1 + a0_Previous;
y10 = r1 + Sol_Total(1,15) * sin(Sol_Total(1,16) * x_curve1 + Sol_Total(1,17)) ...
    + Sol_Total(1,18) * sin(Sol_Total(1,19) * x_curve1 + Sol_Total(1,20));

r2 = a1 * x_curve2 + a0;
y11 = r2 + Sol_Total(1,24) * sin(Sol_Total(1,25) * x_curve2 + Sol_Total(1,26)) ...
    + Sol_Total(1,27) * sin(Sol_Total(1,28) * x_curve2 + Sol_Total(1,29));

if YForecast(1,3) > MM
    figure()
    p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
    hold on
    p2 = plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
    hold on
elseif Second_Forecast_DT(1,1) > MM
    figure()
    p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
    hold on
    p2 = plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
    hold on
elseif First_Forecast_DT(1,1) > MM
    figure()
    p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
    hold on
    p2 = plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
    hold on
else
    figure()
    %p1 = plot(A(1:MM,8), Typ_Pri(1:MM),'k.');
    p1 = plot(x_curve1,Typ_Pri(Sol_Total(1,1):YForecast(1,3)),'k.');
    hold on
    bi = 7;
    while ( bi < length(tt_T) ) && ( tt_T(bi) < YForecast(1,3) )
        bi = bi + 1;
    end
    p2 = plot(tt_T(1:bi),b_T(1:bi),'mo','linewidth',2);
    hold on
end
p3 = plot(x_curve1,y10,'b');
hold on
p4 = plot(x_curve2,y11,'g');
title('Systems 1 and 2: r+a_1*sin(w_1t+\theta_1)+a_2*sin(w_2t+\theta_2)')
% legend('Real Data','Six Points','First System','Second System')
%axis tight
Y_1 = ylim;

for i = 1 : 3 % Additional
    hold on
    X_1 = [YForecast(1,i), YForecast(1,i)];
    p5 = line(X_1,Y_1,'Color',[.5 .3 1.0]);%[.6 .1 .9]
    hold on
    p6 = plot(YForecast(1,i),PForecast(1,i),'ch','linewidth',2); % Additional
end
hold on
p7 = plot(x_maxF,MaxTyp_PriForecast,'gp','linewidth',1.5); % Additional
hold on
plot(x_minF,MinTyp_PriForecast,'gp','linewidth',1.5); % Additional
hold on
X_1 = [Spark(7), Spark(7)];
p8 = line(X_1,Y_1,'Color',[.3 .9 .9]); % Additional
% hold on
% p9 = plot(Gravity_tF,Gravity_PF,'r>','linewidth',1.5);%p9.Marker = '>';
% hold on
% X_1 = [Gravity_tF, Gravity_tF];
% p10 = line(X_1,Y_1,'Color',[1.0 .2 1.0]);%[.9 .1 .9]
% p10.LineStyle = '--';
% hold on
% p11 = plot(x_maxFG+Spark(7),MaxTyp_PriceFG,'g^','linewidth',1.5); % Additional
hold on
plot(x_minFG+Spark(7),MinTyp_PriceFG,'g^','linewidth',1.5); % Additional
legend([p1,p2,p3,p4,p5,p6,p7(1,1),p8],'Real Data','Six Points','First System',...
    'Second System','Encounter of Two System','Encounter of Two System',...
    'Min & Max Real Data','Hunt')
% legend([p1,p2,p3,p4,p5,p6,p7(1,1),p8,p9,p10,p11(1,1)],'Real Data','Six Points','First System',...
%     'Second System','Encounter of Two System','Encounter of Two System',...
%     'Min & Max Real Data','Hunt','Center of gravity','Center of gravity',...
%     'Min & Max Real Data (gravity)')
hold on
plot(Spark(7),Typ_Pri(Spark(7)),'rh','linewidth',2) % Additional

if (Spark(7) ~= 0) && (TypePosition_Sel_R2 ~= -1) && (Typ_Pri(Spark(7)) < PForecast(1,Cut_Number))
    X_2 = xlim;
    New_lineP = (PForecast(1,Cut_Number) + Typ_Pri(Spark(7))) / 2  + 4 * 10^(-(Digit-1)); % TP = 50
    New_lineP_D = Typ_Pri(Spark(7)) - abs(PForecast(1,Cut_Number) - Typ_Pri(Spark(7))) / 2 + 4 * 10^(-(Digit-1)); % SL = 50
    
    Y_2 = [New_lineP, New_lineP];
    Y_3 = [New_lineP_D, New_lineP_D];
    Y_4 = [Hunt_Price + 4 * 10^(-(Digit-1)), Hunt_Price + 4 * 10^(-(Digit-1))];
    hold on
    line(X_2,Y_2,'Color',[.5 .3 1.0]); % Hunt Up
    hold on
    line(X_2,Y_3,'Color',[1.0 .2 1.0]);
    hold on
    line(X_2,Y_4,'Color',[.3 .9 .9]);
    hold on
    plot(TimeClose_Sel_R2,PriceClose_Sel_R2,'b^','linewidth',1.5); % Additional
elseif (Spark(7) ~= 0) && (TypePosition_Sel_R2 ~= -1) && (Typ_Pri(Spark(7)) > PForecast(1,Cut_Number))
    X_2 = xlim;
    New_lineP = (PForecast(1,Cut_Number) + Typ_Pri(Spark(7))) / 2  - 4 * 10^(-(Digit-1)); % TP = 50
    New_lineP_D = Typ_Pri(Spark(7)) + abs(PForecast(1,Cut_Number) - Typ_Pri(Spark(7))) / 2 - 4 * 10^(-(Digit-1)); % SL = 50
    
    Y_2 = [New_lineP, New_lineP];
    Y_3 = [New_lineP_D, New_lineP_D];
    Y_4 = [Hunt_Price - 4 * 10^(-(Digit-1)), Hunt_Price - 4 * 10^(-(Digit-1))];
    hold on
    line(X_2,Y_2,'Color',[.5 .3 1.0]); % Hunt Down
    hold on
    line(X_2,Y_3,'Color',[1.0 .2 1.0]);
    hold on
    line(X_2,Y_4,'Color',[.3 .9 .9]);
    hold on
    plot(TimeClose_Sel_R2,PriceClose_Sel_R2,'b^','linewidth',1.5); % Additional
elseif (TypePosition_Sel_R2 == 100)
    X_2 = xlim;
    New_lineP = (PForecast(1,Cut_Number) + Typ_Pri(Spark(7))) / 2  - 4 * 10^(-(Digit-1)); % TP = 50
    New_lineP_D = Typ_Pri(Spark(7)) + abs(PForecast(1,Cut_Number) - Typ_Pri(Spark(7))) / 2 - 4 * 10^(-(Digit-1)); % SL = 50
    
    Y_2 = [New_lineP, New_lineP];
    Y_3 = [New_lineP_D, New_lineP_D];
    Y_4 = [Hunt_Price - 4 * 10^(-(Digit-1)), Hunt_Price - 4 * 10^(-(Digit-1))];
    hold on
    line(X_2,Y_2,'Color',[.5 .3 1.0]); % Hunt Down
    hold on
    line(X_2,Y_3,'Color',[1.0 .2 1.0]);
    hold on
    line(X_2,Y_4,'Color',[.3 .9 .9]);
    hold on
    plot(TimeClose_Sel_R2,PriceClose_Sel_R2,'b^','linewidth',1.5); % Additional
end

end

%---------------------------  Life time (with noise)  -----------------------
%x_curve1 = Sol_Total(1,11) : YForecast(1,3);
x_curve2 = Spark(7) + 1 : (tt_T(7) + 2 * t7_t1(1,1));

y18 =@(t_1) a1_Previous * t_1 + a0_Previous ...
    + Sol_Total(1,15) * sin(Sol_Total(1,16) * t_1 + Sol_Total(1,17)) ...
    + Sol_Total(1,18) * sin(Sol_Total(1,19) * t_1 + Sol_Total(1,20));

y19 =@(t_1) a1 * t_1 + a0 ...
    + Sol_Total(1,24) * sin(Sol_Total(1,25) * t_1 + Sol_Total(1,26))...
    + Sol_Total(1,27) * sin(Sol_Total(1,28) * t_1 + Sol_Total(1,29));

dis1(1) = abs( y18(Sol_Total(1,1)) - Sol_Total(1,2) );
dis1(2) = abs( y18(Sol_Total(1,3)) - Sol_Total(1,4) );
dis1(3) = abs( y18(Sol_Total(1,5)) - Sol_Total(1,6) );
dis1(4) = abs( y18(Sol_Total(1,7)) - Sol_Total(1,8) );
dis1(5) = abs( y18(Sol_Total(1,9)) - Sol_Total(1,10) );
dis1(6) = abs( y18(Sol_Total(1,11)) - Sol_Total(1,12) );
max_dis1 = max(dis1);
average1 =  ( dis1(1) + dis1(2) + dis1(3) + dis1(4) + dis1(5) + dis1(6) ) / 6;

dis2(1) = abs( y19(Sol_Total(1,3)) - Sol_Total(1,4) );
dis2(2) = abs( y19(Sol_Total(1,5)) - Sol_Total(1,6) );
dis2(3) = abs( y19(Sol_Total(1,7)) - Sol_Total(1,8) );
dis2(4) = abs( y19(Sol_Total(1,9)) - Sol_Total(1,10) );
dis2(5) = abs( y19(Sol_Total(1,11)) - Sol_Total(1,12) );
dis2(6) = abs( y19(Sol_Total(1,13)) - Sol_Total(1,14) );
max_dis2 = max(dis2);
average2 =  ( dis2(1) + dis2(2) + dis2(3) + dis2(4) + dis2(5) + dis2(6) ) / 6;

y24 = y18(x_curve2) + (max_dis1 + average1) / 2;
y25 = y18(x_curve2) - (max_dis1 + average1) / 2;

y26 = y19(x_curve2) + (max_dis2 + average2) / 2;
y27 = y19(x_curve2) - (max_dis2 + average2) / 2;

LifeTime1 = -1; LifeTime2 = -1;
if (tt_T(7) + 2 * t7_t1(1,1)) <= MM
    for i = Spark(7) + 1 : (tt_T(7) + 2 * t7_t1(1,1)) - 1
        if ( Typ_Pri(i+1) >= y24(i-Spark(7)+1) ) && ( Typ_Pri(i) < y24(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y24(i-Spark(7)+1) ) && ( Typ_Pri(i) > y24(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y25(i-Spark(7)+1) ) && ( Typ_Pri(i) > y25(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y25(i-Spark(7)+1) ) && ( Typ_Pri(i) < y25(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        end
    end
else
    for i = Spark(7) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y24(i-Spark(7)+1) ) && ( Typ_Pri(i) < y24(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y24(i-Spark(7)+1) ) && ( Typ_Pri(i) > y24(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y25(i-Spark(7)+1) ) && ( Typ_Pri(i) > y25(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y25(i-Spark(7)+1) ) && ( Typ_Pri(i) < y25(i-Spark(7)) )
            LifeTime1 = i+1;
            break
        end
    end
end
if (tt_T(7) + 2 * t7_t1(1,1)) <= MM
    for i = Spark(7) + 1 : (tt_T(7) + 2 * t7_t1(1,1)) - 1
        if ( Typ_Pri(i+1) >= y26(i-Spark(7)+1) ) && ( Typ_Pri(i) < y26(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y26(i-Spark(7)+1) ) && ( Typ_Pri(i) > y26(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y27(i-Spark(7)+1) ) && ( Typ_Pri(i) > y27(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y27(i-Spark(7)+1) ) && ( Typ_Pri(i) < y27(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        end
    end
else
    for i = Spark(7) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y26(i-Spark(7)+1) ) && ( Typ_Pri(i) < y26(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y26(i-Spark(7)+1) ) && ( Typ_Pri(i) > y26(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y27(i-Spark(7)+1) ) && ( Typ_Pri(i) > y27(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y27(i-Spark(7)+1) ) && ( Typ_Pri(i) < y27(i-Spark(7)) )
            LifeTime2 = i+1;
            break
        end
    end
end
if (LifeTime1 == -1) && (LifeTime2 == -1)
    LifeTimeB = -1;
elseif (LifeTime1 ~= -1) && (LifeTime2 == -1)
    LifeTimeB = LifeTime1; % Life time for with nise
elseif (LifeTime1 == -1) && (LifeTime2 ~= -1)
    LifeTimeB = LifeTime2; % Life time for with nise
elseif (LifeTime1 ~= -1) && (LifeTime2 ~= -1)
    LifeTimeB = max(LifeTime1,LifeTime2); % Life time for with nise
else
    LifeTimeB = -1;
end

% hold on
% plot(x_curve2,y24,'b-.',x_curve2,y25,'b-.');
% hold on
% plot(x_curve2,y26,'g-.',x_curve2,y27,'g-.');
% hold on
% if LifeTime1 ~= -1
%     p12 = line(LifeTime1,Typ_Pri(LifeTime1),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p12.Marker = 's';
% end
% hold on
% if LifeTime2 ~= -1
%     p13 = line(LifeTime2,Typ_Pri(LifeTime2),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p13.Marker = 's';
% end

%-------------------------------  Additional  -----------------------------
% x_curve3 = Sol_Total(1,1) : 2*( tt_T(7) + t7_t1(1,1) );
% x_curve4 = Sol_Total(1,3) : 2*( tt_T(7) + t7_t1(1,1) );
% 
% y12 = a1_Previous * x_curve3 + a0_Previous ...
%     + Sol_Total(1,15) * sin(Sol_Total(1,16) * x_curve3 + Sol_Total(1,17));
% 
% y13 = a1 * x_curve4 + a0 ...
%     + Sol_Total(1,24) * sin(Sol_Total(1,25) * x_curve4 + Sol_Total(1,26));
% 
% if 2*( tt_T(7) + t7_t1(1,1) ) > MM
%     figure()
%     plot(A(Sol_Total(1,1):MM,8), Typ_Pri(Sol_Total(1,1):MM),'k.');
%     hold on
%     plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
%     hold on
% else
%     figure()
%     plot(Sol_Total(1,1):2*( tt_T(7) + t7_t1(1,1) ),Typ_Pri(Sol_Total(1,1):2*( tt_T(7) + t7_t1(1,1) )),'k.');
%     hold on
%     bi = 7;
%     while ( bi < length(tt_T) ) && ( tt_T(bi) < 2*( tt_T(7) + t7_t1(1,1) ))
%         bi = bi + 1;
%     end
%     plot(tt_T(1:bi),b_T(1:bi),'mo','linewidth',2);
%     hold on
% end
% plot(x_curve3,y12,'b',x_curve4,y13,'g');
% 
% % With noise
% y14 = a1_Previous * x_curve3 + a0_Previous ...
%     + Sol_Total(1,15) * sin(Sol_Total(1,16) * x_curve3 + Sol_Total(1,17)) ...
%     + Sol_Total(1,18) * sin(Sol_Total(1,19) * x_curve3 + Sol_Total(1,20));
% 
% y15 = a1 * x_curve4 + a0 ...
%     + Sol_Total(1,24) * sin(Sol_Total(1,25) * x_curve4 + Sol_Total(1,26)) ...
%     + Sol_Total(1,27) * sin(Sol_Total(1,28) * x_curve4 + Sol_Total(1,29));
% 
% if 2*( tt_T(7) + t7_t1(1,1) ) > MM
%     figure()
%     plot(A(Sol_Total(1,1):MM,8), Typ_Pri(Sol_Total(1,1):MM),'k.');
%     hold on
%     plot(tt_T(1:length(tt_T)),b_T(1:length(b_T)),'mo','linewidth',2);
%     hold on
% else
%     figure()
%     plot(Sol_Total(1,1):2*( tt_T(7) + t7_t1(1,1) ),Typ_Pri(Sol_Total(1,1):2*( tt_T(7) + t7_t1(1,1) )),'k.');
%     hold on
%     bi = 7;
%     while ( bi < length(tt_T) ) && ( tt_T(bi) < 2*( tt_T(7) + t7_t1(1,1) ))
%         bi = bi + 1;
%     end
%     plot(tt_T(1:bi),b_T(1:bi),'mo','linewidth',2);
%     hold on
% end
% plot(x_curve3,y14,'b',x_curve4,y15,'g');

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%---------------------------  ii = 3 : SIZE - 6 + 1  ----------------------
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for ii = 3 : SIZE - m_1 + 1
    for i = 1 : m_1 + 1
        Sol_Total(ii-1,i*2-1) = tt_T(ii+i-2); % Discrete time
        Sol_Total(ii-1,i*2) = b_T(ii+i-2); % Price
        Data_DT2(ii-1,i) = TIME(ii+i-2);
        Time_Spark(ii-1,i) = Spark(ii+i-2);
        Time_Delay(ii-1,i) = Delay(ii+i-2);
    end

Sol_Total(ii-1,15:23) = Sol_Total(ii-2,24:32);

tt = tt_T(ii:ii+5);
b = b_T(ii:ii+5);

x = A(tt(1) : tt(6),8); % Discrete time
y = Typ_Pri(tt_T(ii):tt_T(ii+5)); % Typical Price

a1_Previous = a1;
a0_Previous = a0;
[a1,a0] = regression(x,y);
%--------------------------------  (a_1,w_1)  -----------------------------
ttt(1:m_1-1) = ( tt(1:m_1-1) + tt(2:m_1) ) / 2;
bb(1:m_1-1) = ( b(1:m_1-1) + b(2:m_1) ) / 2;

S1=@(y, t, m, n) y-m*t-n;
S2=@(t, w) sin(w * t);
S1_1=@(y, t) S1(y, t, a1, a0);

r = a1 * x + a0;

lambda = 2;
[solution_1,Maximum] = SELECT_5(m_1, M, lambda, ttt, bb, x, r, S1_1, S2);
while Maximum > 2
    lambda = lambda - 0.1;
    [solution_1,Maximum] = SELECT_5(m_1, M, lambda, ttt, bb, x, r,S1_1, S2);
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if solution_1(2) < 10^(-4)
    HasMid = 0;
    lambda = 2;
    [y5, s_1_1, s_1_2, THETA_1, s_2_1, s_2_2, THETA_2, lambda, a4] = ...
        PLOT_APS_6(m_1, M, lambda, tt, b, S1, S2, S1_1, x, a1, a0, r);
    Sol_Total(ii-1,24:32) = [s_1_1, s_1_2, THETA_1, s_2_1, s_2_2, THETA_2, lambda, a4-1, HasMid];
else
    HasMid = 1;
    %---------------------------------  Theta_1  --------------------------
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
            S32(j) = abs( S1(b(j), tt(j), a1, a0) - solution_1(1) * ...
                sin(solution_1(2)*tt(j)+Index_3(k2,1)) );
            SUM_1(k2) = SUM_1(k2) + S32(j);
        end
    end
    Minimum_1 = min(SUM_1);
    [Index_6,Index_7] = find(SUM_1==Minimum_1);
    THETA_1 = Index_3(Index_7(1),1); % Solution Set
    
    %--------------------------------  (a_2,w_2)  -------------------------
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
    else
        k6 = k5-1;
    end
    %---------------------------------  Theta_2  --------------------------
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
            S34(j) = abs( S3(b(j), tt(j)) - solution_2(1,k6) * ...
                sin(solution_2(2,k6)*tt(j)+Index_11(k4,1)) );
            SUM_3(k4) = SUM_3(k4) + S34(j);
        end
    end
    Minimum_3 = min(SUM_3);
    [Index_12,Index_13] = find(SUM_3==Minimum_3);
    THETA_2 = Index_11(Index_13(1),1); % Solution Set
    
    y5 = r + solution_1(1)*sin(solution_1(2)*x+THETA_1) + ...
        solution_2(1,k6)*sin(solution_2(2,k6)*x+THETA_2);
    
    Sol_Total(ii-1,24:32) = [solution_1(1), solution_1(2), THETA_1, solution_2(1,k6), ...
        solution_2(2,k6), THETA_2, lambda, a4-1, HasMid];
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Sign Regression  >>>>>>>>>>>>>>>>>>>>>>>>>>
Sig_Reg(ii) = sign(a1);

Slope(ii,1) = ( b_T(ii+2) - b_T(ii) ) / ( tt_T(ii+2) - tt_T(ii) ); % 1_3
i3 = 0;
i4 = 0;
Slo_2 = {};
Slo_7 = zeros();
if sign(Slope(ii,1)) == Sig_Reg(ii)
    Slo_5(1) = ii; % Additional
    for i = ii : ii + m_1 - 1
        if i == ii
            continue
        elseif i == ii+2
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(ii,1) * ( tt_T(i) - tt_T(ii) ) + b_T(ii);
            if y_prim <= b_T(i) && b_T(ii) < b_T(ii+1) % min
                Slo_1(1,i3) = 1;
                Slo_6(ii,1) = 1001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(ii) > b_T(ii+1) % max
                Slo_1(1,i3) = 1;
                Slo_6(ii,1) = 1002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(1,i3) = 0;
                Slo_6(ii,1) = 1005; % Additional
            end
        end
    end
    Slo_3(ii,1) = CONT_R_2(Slo_1(1,:), size(Slo_1,2));
    if Slo_3(ii,1) == 4 && Sig_Reg(ii) == -1
        i4 = i4 + 1;
        Slo_2(1,1) = {'1_3'};
        Slo_2(1,2) = {'4'};
        Slo_2(1,3) = {'Descending'};
        Slo_2(1,4) = {'1'};
        Slo_2(1,5) = Slo_9;
        Slo_2(1,6) = Slo_10;
        Slo_2(1,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,1));
        Slo_7(i4,2) = 1;
    elseif Slo_3(ii,1) == 4 && Sig_Reg(ii) == 1
        i4 = i4 + 1;
        Slo_2(1,1) = {'1_3'};
        Slo_2(1,2) = {'4'};
        Slo_2(1,3) = {'Ascending'};
        Slo_2(1,4) = {'2'};
        Slo_2(1,5) = Slo_9;
        Slo_2(1,6) = Slo_10;
        Slo_2(1,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,1));
        Slo_7(i4,2) = 1;
    end
end
%--------------------------------------------------------------------------
Slope(ii,2) = ( b_T(ii+3) - b_T(ii+1) ) / ( tt_T(ii+3) - tt_T(ii+1) ); % 2_4
i3 = 0;
if sign(Slope(ii,2)) == Sig_Reg(ii)
    Slo_5(2) = ii+1; % Additional
    for i = ii+1 : ii + m_1 - 1
        if i == ii+1
            continue
        elseif i == ii+3
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(ii,2) * ( tt_T(i) - tt_T(ii+1) ) + b_T(ii+1);
            if y_prim >= b_T(i) && b_T(ii) < b_T(ii+1) % max
                Slo_1(2,i3) = 1;
                Slo_6(ii,2) = 2001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(ii) > b_T(ii+1) % min
                Slo_1(2,i3) = 1;
                Slo_6(ii,2) = 2002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(2,i3) = 0;
                Slo_6(ii,2) = 2005; % Additional
            end
        end
    end
    Slo_3(ii,2) = CONT_R_2(Slo_1(2,:), size(Slo_1,2));
    if Slo_3(ii,2) == 3 && Sig_Reg(ii) == -1
        i4 = i4 + 1;
        Slo_2(2,1) = {'2_4'};
        Slo_2(2,2) = {'3'};
        Slo_2(2,3) = {'Descending'};
        Slo_2(2,4) = {'1'};
        Slo_2(2,5) = Slo_9;
        Slo_2(2,6) = Slo_10;
        Slo_2(2,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,2));
        Slo_7(i4,2) = 2;
    elseif Slo_3(ii,2) == 3 && Sig_Reg(ii) == 1
        i4 = i4 + 1;
        Slo_2(2,1) = {'2_4'};
        Slo_2(2,2) = {'3'};
        Slo_2(2,3) = {'Ascending'};
        Slo_2(2,4) = {'2'};
        Slo_2(2,5) = Slo_9;
        Slo_2(2,6) = Slo_10;
        Slo_2(2,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,2));
        Slo_7(i4,2) = 2;
    end
end
%--------------------------------------------------------------------------
Slope(ii,3) = ( b_T(ii+4) - b_T(ii+2) ) / ( tt_T(ii+4) - tt_T(ii+2) ); % 3_5
i3 = 0;
if sign(Slope(ii,3)) == Sig_Reg(ii)
    Slo_5(3) = ii+2; % Additional
    for i = ii+2 : ii + m_1 - 1
        if i == ii+2
            continue
        elseif i == ii+4
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(ii,3) * ( tt_T(i) - tt_T(ii+2) ) + b_T(ii+2);
            if y_prim <= b_T(i) && b_T(ii) < b_T(ii+1) % min
                Slo_1(3,i3) = 1;
                Slo_6(ii,3) = 3001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(ii) > b_T(ii+1) % max
                Slo_1(3,i3) = 1;
                Slo_6(ii,3) = 3002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(3,i3) = 0;
                Slo_6(ii,3) = 3005; % Additional
            end
        end
    end
    Slo_3(ii,3) = CONT_R_2(Slo_1(3,:), size(Slo_1,2));
    if Slo_3(ii,3) == 2 && Sig_Reg(ii) == -1
        i4 = i4 + 1;
        Slo_2(3,1) = {'3_5'};
        Slo_2(3,2) = {'2'};
        Slo_2(3,3) = {'Descending'};
        Slo_2(3,4) = {'1'};
        Slo_2(3,5) = Slo_9;
        Slo_2(3,6) = Slo_10;
        Slo_2(3,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,3));
        Slo_7(i4,2) = 3;
    elseif Slo_3(ii,3) == 2 && Sig_Reg(ii) == 1
        i4 = i4 + 1;
        Slo_2(3,1) = {'3_5'};
        Slo_2(3,2) = {'2'};
        Slo_2(3,3) = {'Ascending'};
        Slo_2(3,4) = {'2'};
        Slo_2(3,5) = Slo_9;
        Slo_2(3,6) = Slo_10;
        Slo_2(3,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,3));
        Slo_7(i4,2) = 3;
    end
end
%--------------------------------------------------------------------------
Slope(ii,4) = ( b_T(ii+5) - b_T(ii+3) ) / ( tt_T(ii+5) - tt_T(ii+3) ); % 4_6
i3 = 0;
if sign(Slope(ii,4)) == Sig_Reg(ii)
    Slo_5(4) = ii+3; % Additional
    for i = ii+3 : ii + m_1 - 1
        if i == ii+3
            continue
        elseif i == ii+5
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(ii,4) * ( tt_T(i) - tt_T(ii+3) ) + b_T(ii+3);
            if y_prim >= b_T(i) && b_T(ii) < b_T(ii+1) % max
                Slo_1(4,i3) = 1;
                Slo_6(ii,4) = 4001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(ii) > b_T(ii+1) % min
                Slo_1(4,i3) = 1;
                Slo_6(ii,4) = 4002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(4,i3) = 0;
                Slo_6(ii,4) = 4005; % Additional
            end
        end
    end
    Slo_3(ii,4) = CONT_R_2(Slo_1(4,:), size(Slo_1,2));
    if Slo_3(ii,4) == 1 && Sig_Reg(ii) == -1
        i4 = i4 + 1;
        Slo_2(4,1) = {'4_6'};
        Slo_2(4,2) = {'1'};
        Slo_2(4,3) = {'Descending'};
        Slo_2(4,4) = {'1'};
        Slo_2(4,5) = Slo_9;
        Slo_2(4,6) = Slo_10;
        Slo_2(4,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,4));
        Slo_7(i4,2) = 4;
    elseif Slo_3(ii,4) == 1 && Sig_Reg(ii) == 1
        i4 = i4 + 1;
        Slo_2(4,1) = {'4_6'};
        Slo_2(4,2) = {'1'};
        Slo_2(4,3) = {'Ascending'};
        Slo_2(4,4) = {'2'};
        Slo_2(4,5) = Slo_9;
        Slo_2(4,6) = Slo_10;
        Slo_2(4,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,4));
        Slo_7(i4,2) = 4;
    end
end
%--------------------------------------------------------------------------
Slope(ii,5) = ( b_T(ii+4) - b_T(ii) ) / ( tt_T(ii+4) - tt_T(ii) ); % 1_5
i3 = 0;
if sign(Slope(ii,5)) == Sig_Reg(ii)
    Slo_5(5) = ii; % Additional
    for i = ii : ii + m_1 - 1
        if i == ii
            continue
        elseif i == ii+4
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(ii,5) * ( tt_T(i) - tt_T(ii) ) + b_T(ii);
            if y_prim <= b_T(i) && b_T(ii) < b_T(ii+1) % min
                Slo_1(5,i3) = 1;
                Slo_6(ii,5) = 5001; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            elseif y_prim >= b_T(i) && b_T(ii) > b_T(ii+1) % max
                Slo_1(5,i3) = 1;
                Slo_6(ii,5) = 5002; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            else
                Slo_1(5,i3) = 0;
                Slo_6(ii,5) = 5005; % Additional
            end
        end
    end
    Slo_3(ii,5) = CONT_R_2(Slo_1(5,:), size(Slo_1,2));
    if Slo_3(ii,5) == 4 && Sig_Reg(ii) == -1
        i4 = i4 + 1;
        Slo_2(5,1) = {'1_5'};
        Slo_2(5,2) = {'4'};
        Slo_2(5,3) = {'Descending'};
        Slo_2(5,4) = {'1'};
        Slo_2(5,5) = Slo_9;
        Slo_2(5,6) = Slo_10;
        Slo_2(5,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,5));
        Slo_7(i4,2) = 5;
    elseif Slo_3(ii,5) == 4 && Sig_Reg(ii) == 1
        i4 = i4 + 1;
        Slo_2(5,1) = {'1_5'};
        Slo_2(5,2) = {'4'};
        Slo_2(5,3) = {'Ascending'};
        Slo_2(5,4) = {'2'};
        Slo_2(5,5) = Slo_9;
        Slo_2(5,6) = Slo_10;
        Slo_2(5,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,5));
        Slo_7(i4,2) = 5;
    end
end
%--------------------------------------------------------------------------
Slope(ii,6) = ( b_T(ii+5) - b_T(ii+1) ) / ( tt_T(ii+5) - tt_T(ii+1) ); % 2_6
i3 = 0;
if sign(Slope(ii,6)) == Sig_Reg(ii)
    Slo_5(6) = ii+1; % Additional
    for i = ii+1 : ii + m_1 - 1
        if i == ii+1
            continue
        elseif i == ii+5
            continue
        else
            i3 = i3 + 1;
            y_prim = Slope(ii,6) * ( tt_T(i) - tt_T(ii+1) ) + b_T(ii+1);
            if y_prim >= b_T(i) && b_T(ii) < b_T(ii+1) % max
                Slo_1(6,i3) = 1;
                Slo_6(ii,6) = 6001; % Additional
                Slo_9 = {'Down'};
                Slo_10 = {'1'};
                Slo_11 = {'Max'};
            elseif y_prim <= b_T(i) && b_T(ii) > b_T(ii+1) % min
                Slo_1(6,i3) = 1;
                Slo_6(ii,6) = 6002; % Additional
                Slo_9 = {'Up'};
                Slo_10 = {'2'};
                Slo_11 = {'Min'};
            else
                Slo_1(6,i3) = 0;
                Slo_6(ii,6) = 6005; % Additional
            end
        end
    end
    Slo_3(ii,6) = CONT_R_2(Slo_1(6,:), size(Slo_1,2));
    if Slo_3(ii,6) == 3 && Sig_Reg(ii) == -1
        i4 = i4 + 1;
        Slo_2(6,1) = {'2_6'};
        Slo_2(6,2) = {'3'};
        Slo_2(6,3) = {'Descending'};
        Slo_2(6,4) = {'1'};
        Slo_2(6,5) = Slo_9;
        Slo_2(6,6) = Slo_10;
        Slo_2(6,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,6));
        Slo_7(i4,2) = 6;
    elseif Slo_3(ii,6) == 3 && Sig_Reg(ii) == 1
        i4 = i4 + 1;
        Slo_2(6,1) = {'2_6'};
        Slo_2(6,2) = {'3'};
        Slo_2(6,3) = {'Ascending'};
        Slo_2(6,4) = {'2'};
        Slo_2(6,5) = Slo_9;
        Slo_2(6,6) = Slo_10;
        Slo_2(6,7) = Slo_11;
        Slo_7(i4,1) = abs(a1 - Slope(ii,6));
        Slo_7(i4,2) = 6;
    end
end
%--------------------------------------------------------------------------

if Slo_7 == 0
    Trend_Num_1(ii-1,1) = Trend_Num_2(ii-2,1);
    Point_Num_1(ii-1,1) = Point_Num_2(ii-2,1);
    Trend_1(ii-1,1) = Trend_2(ii-2,1);
    Trend_Mac_1(ii-1,1) = Trend_Mac_2(ii-2,1);
    Up_Down_1(ii-1,1) = Up_Down_2(ii-2,1);
    Up_Down_Mac_1(ii-1,1) = Up_Down_Mac_2(ii-2,1);
    Min_Max_1(ii-1,1) = Min_Max_2(ii-2,1);
    
    Trend_Num_2(ii-1,1) = {'0'};
    Point_Num_2(ii-1,1) = {'0'};
    Trend_2(ii-1,1) = {'No Terend'};
    Trend_Mac_2(ii-1,1) = {'0'};
    Up_Down_2(ii-1,1) = {'No'};
    Up_Down_Mac_2(ii-1,1) = {'0'};
    Min_Max_2(ii-1,1) = {'No'};
else
    Slo_4 = min(Slo_7(:,1));
    ind_slo = find( Slo_7(:,1) == Slo_4 );
    Slo_8 = Slo_7(ind_slo,2);
    
%     for i = 1 : length(ind_slo)
%         plot(x,Slope(ii,Slo_8) * (x-tt_T(Slo_5(Slo_8))) + b_T(Slo_5(Slo_8)),'b') % Additional
%         hold on
%     end
    
    Trend_Num_1(ii-1,1) = Trend_Num_2(ii-2,1);
    Point_Num_1(ii-1,1) = Point_Num_2(ii-2,1);
    Trend_1(ii-1,1) = Trend_2(ii-2,1);
    Trend_Mac_1(ii-1,1) = Trend_Mac_2(ii-2,1);
    Up_Down_1(ii-1,1) = Up_Down_2(ii-2,1);
    Up_Down_Mac_1(ii-1,1) = Up_Down_Mac_2(ii-2,1);
    Min_Max_1(ii-1,1) = Min_Max_2(ii-2,1);
    
    Trend_Num_2(ii-1,1) = Slo_2(Slo_8,1);
    Point_Num_2(ii-1,1) = Slo_2(Slo_8,2);
    Trend_2(ii-1,1) = Slo_2(Slo_8,3);
    Trend_Mac_2(ii-1,1) = Slo_2(Slo_8,4);
    Up_Down_2(ii-1,1) = Slo_2(Slo_8,5);
    Up_Down_Mac_2(ii-1,1) = Slo_2(Slo_8,6);
    Min_Max_2(ii-1,1) = Slo_2(Slo_8,7);
end
if Sig_Reg(ii-1) == Sig_Reg(ii)
    Trend_Change(ii-1,1) = 0;
else
    Trend_Change(ii-1,1) = 1;
end

t7_t1(ii-1,1) = Sol_Total(ii-1,13) - Sol_Total(ii-1,1) + 1;

% plot(x,r)
% hold on

% Slo_2
% Slo_7(:,1)
% Slo_7(:,2)

%--------------------------------  Forecast  ------------------------------
First_t = tt_T(ii+5);
NumerationE = 0;
DeltaY1 = a1_Previous * First_t + a0_Previous ...
    + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * First_t + Sol_Total(ii-1,17)) ...
    - a1 * First_t - a0 ...
    - Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26));

while (NumerationE < 3) && ( First_t <= (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) )
    First_t = First_t + 1;
    DeltaY2 = a1_Previous * First_t + a0_Previous ...
        + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * First_t + Sol_Total(ii-1,17)) ...
        - a1 * First_t - a0 ...
        - Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26));
    
    if ( sign(DeltaY1) ~= sign(DeltaY2)) && (First_t <= MM)
        NumerationE = NumerationE + 1;
        NY1(ii-1,NumerationE) = Row(First_t,1);
        NDeltaY(ii-1,NumerationE) = First_t;
        PEncounter(ii-1,NumerationE) = a1 * First_t + a0 + Sol_Total(ii-1,24) ...
            * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26));
    elseif ( sign(DeltaY1) ~= sign(DeltaY2)) && (First_t > MM)
        NumerationE = NumerationE + 1;
        NY1(ii-1,NumerationE) = {'Infinity'};
        NDeltaY(ii-1,NumerationE) = First_t;
        PEncounter(ii-1,NumerationE) = a1 * First_t + a0 + Sol_Total(ii-1,24) ...
            * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26));
    end
    DeltaY1 = DeltaY2;
end

if (NumerationE==0) || (NumerationE==1) || (NumerationE==2) % The # encounters is not enough.
    First_encounter_CT(ii-1,1) = {'Infinity'};
    First_encounter_DT(ii-1,1) = -1;
    Second_encounter_CT(ii-1,1) = {'Infinity'};
    Second_encounter_DT(ii-1,1) = -1;
    Third_encounter_CT(ii-1,1) = {'Infinity'};
    Third_encounter_DT(ii-1,1) = -1;
    
    First_PEncounter(ii-1,1) = -1;
    Second_PEncounter(ii-1,1) = -1;
    Third_PEncounter(ii-1,1) = -1;
    
    XNUD_Encounter(ii-1,1) = {'NO'};
    Hunt_Encounter(ii-1,1) = {'NO'};
    
    MaxTyp_PriEncounter(ii-1,1) = -1;
    MinTyp_PriEncounter(ii-1,1) = -1;
    
    Max_PercentE(ii-1,1) = -1;
    Deviation_E(ii-1,1) = -1;
    
    TypePosition_E(ii-1,1) = -1;
    TypePos_MaxPer_Dev_E(ii-1,1) = -1;
    TypePos_E(ii-1,1) = {'NO'};
elseif NDeltaY(ii-1,1) > MM % Because there is no real data.
    First_encounter_CT(ii-1,1) = {'Infinity'};
    First_encounter_DT(ii-1,1) = -1;
    Second_encounter_CT(ii-1,1) = {'Infinity'};
    Second_encounter_DT(ii-1,1) = -1;
    Third_encounter_CT(ii-1,1) = {'Infinity'};
    Third_encounter_DT(ii-1,1) = -1;
    
    First_PEncounter(ii-1,1) = -1;
    Second_PEncounter(ii-1,1) = -1;
    Third_PEncounter(ii-1,1) = -1;
    
    XNUD_Encounter(ii-1,1) = {'NO'};
    Hunt_Encounter(ii-1,1) = {'NO'};
    
    MaxTyp_PriEncounter(ii-1,1) = -1;
    MinTyp_PriEncounter(ii-1,1) = -1;
    
    Max_PercentE(ii-1,1) = -1;
    Deviation_E(ii-1,1) = -1;
    
    TypePosition_E(ii-1,1) = -1;
    TypePos_MaxPer_Dev_E(ii-1,1) = -1;
    TypePos_E(ii-1,1) = {'NO'};
elseif Spark(ii+5) == 0 % The extreme point is not clear.
    First_encounter_CT(ii-1,1) = {'Infinity'};
    First_encounter_DT(ii-1,1) = -1;
    Second_encounter_CT(ii-1,1) = {'Infinity'};
    Second_encounter_DT(ii-1,1) = -1;
    Third_encounter_CT(ii-1,1) = {'Infinity'};
    Third_encounter_DT(ii-1,1) = -1;
    
    First_PEncounter(ii-1,1) = -1;
    Second_PEncounter(ii-1,1) = -1;
    Third_PEncounter(ii-1,1) = -1;
    
    XNUD_Encounter(ii-1,1) = {'NO'};
    Hunt_Encounter(ii-1,1) = {'NO'};
    
    MaxTyp_PriEncounter(ii-1,1) = -1;
    MinTyp_PriEncounter(ii-1,1) = -1;
    
    Max_PercentE(ii-1,1) = -1;
    Deviation_E(ii-1,1) = -1;
    
    TypePosition_E(ii-1,1) = -1;
    TypePos_MaxPer_Dev_E(ii-1,1) = -1;
    TypePos_E(ii-1,1) = {'NO'};
elseif NDeltaY(ii-1,1) <= Spark(ii+5) % The first encounter is before the first spark.
    First_encounter_CT(ii-1,1) = {'Infinity'};
    First_encounter_DT(ii-1,1) = -1;
    Second_encounter_CT(ii-1,1) = {'Infinity'};
    Second_encounter_DT(ii-1,1) = -1;
    Third_encounter_CT(ii-1,1) = {'Infinity'};
    Third_encounter_DT(ii-1,1) = -1;
    
    First_PEncounter(ii-1,1) = -1;
    Second_PEncounter(ii-1,1) = -1;
    Third_PEncounter(ii-1,1) = -1;
    
    XNUD_Encounter(ii-1,1) = {'NO'};
    Hunt_Encounter(ii-1,1) = {'NO'};
    
    MaxTyp_PriEncounter(ii-1,1) = -1;
    MinTyp_PriEncounter(ii-1,1) = -1;
    
    Max_PercentE(ii-1,1) = -1;
    Deviation_E(ii-1,1) = -1;
    
    TypePosition_E(ii-1,1) = -1;
    TypePos_MaxPer_Dev_E(ii-1,1) = -1;
    TypePos_E(ii-1,1) = {'NO'};
elseif (NDeltaY(ii-1,1) - Spark(ii+5)) < 1 % There is no data between first encounter and spark timing.
    First_encounter_CT(ii-1,1) = {'Infinity'};
    First_encounter_DT(ii-1,1) = -1;
    Second_encounter_CT(ii-1,1) = {'Infinity'};
    Second_encounter_DT(ii-1,1) = -1;
    Third_encounter_CT(ii-1,1) = {'Infinity'};
    Third_encounter_DT(ii-1,1) = -1;
    
    First_PEncounter(ii-1,1) = -1;
    Second_PEncounter(ii-1,1) = -1;
    Third_PEncounter(ii-1,1) = -1;
    
    XNUD_Encounter(ii-1,1) = {'NO'};
    Hunt_Encounter(ii-1,1) = {'NO'};
    
    MaxTyp_PriEncounter(ii-1,1) = -1;
    MinTyp_PriEncounter(ii-1,1) = -1;
    
    Max_PercentE(ii-1,1) = -1;
    Deviation_E(ii-1,1) = -1;
    
    TypePosition_E(ii-1,1) = -1;
    TypePos_MaxPer_Dev_E(ii-1,1) = -1;
    TypePos_E(ii-1,1) = {'NO'};
elseif (NDeltaY(ii-1,1) - tt_T(ii+5) + 1) <= (t7_t1(ii-1,1) / 2)
    First_encounter_CT(ii-1,1) = NY1(ii-1,1);
    First_encounter_DT(ii-1,1) = NDeltaY(ii-1,1);
    Second_encounter_CT(ii-1,1) = NY1(ii-1,2);
    Second_encounter_DT(ii-1,1) = NDeltaY(ii-1,2);
    Third_encounter_CT(ii-1,1) = NY1(ii-1,3);
    Third_encounter_DT(ii-1,1) = NDeltaY(ii-1,3);
    
    First_PEncounter(ii-1,1) = PEncounter(ii-1,1);
    Second_PEncounter(ii-1,1) = PEncounter(ii-1,2);
    Third_PEncounter(ii-1,1) = PEncounter(ii-1,3);
    
    if NDeltaY(ii-1,1) <= MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : NDeltaY(ii-1,1) );
        MaxTyp_PriEncounter(ii-1,1) = max(MTyp_Pri); % Find max
        MinTyp_PriEncounter(ii-1,1) = min(MTyp_Pri); % Find min
        x_maxE = find(MTyp_Pri == MaxTyp_PriEncounter(ii-1,1)) + Spark(ii+5);
        x_minE = find(MTyp_Pri == MinTyp_PriEncounter(ii-1,1)) + Spark(ii+5);
    elseif NDeltaY(ii-1,1) > MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : MM );
        MaxTyp_PriEncounter(ii-1,1) = max(MTyp_Pri); % Find max
        MinTyp_PriEncounter(ii-1,1) = min(MTyp_Pri); % Find min
        x_maxE = find(MTyp_Pri == MaxTyp_PriEncounter(ii-1,1)) + Spark(ii+5);
        x_minE = find(MTyp_Pri == MinTyp_PriEncounter(ii-1,1)) + Spark(ii+5);
    end
    
    if (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) < PEncounter(ii-1,1)) % b_T(7) is max
        XNUD_Encounter(ii-1,1) = {'XU'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) > PEncounter(ii-1,1)) % b_T(7) is max
        XNUD_Encounter(ii-1,1) = {'XD'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) < PEncounter(ii-1,1)) % b_T(7) is min
        XNUD_Encounter(ii-1,1) = {'NU'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) > PEncounter(ii-1,1)) % b_T(7) is min
        XNUD_Encounter(ii-1,1) = {'ND'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) == PEncounter(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Encounter(ii-1,1) = {'XE'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) == PEncounter(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Encounter(ii-1,1) = {'NE'};
    end
    
    if Typ_Pri(Spark(ii+5)) < PEncounter(ii-1,1) % Hunt Up
        Hunt_Encounter(ii-1,1) = {'Hunt_U'};
        if MaxTyp_PriEncounter(ii-1,1) <= Typ_Pri(Spark(ii+5))
            Max_PercentE(ii-1,1) = 0;
        elseif ( Typ_Pri(Spark(ii+5)) < MaxTyp_PriEncounter(ii-1,1) ) ...
                && ( MaxTyp_PriEncounter(ii-1,1) < First_PEncounter(ii-1,1) )
            Max_PercentE(ii-1,1) = abs( MaxTyp_PriEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        elseif MaxTyp_PriEncounter(ii-1,1) >= First_PEncounter(ii-1,1)
            Max_PercentE(ii-1,1) = 100;
        end
        if MinTyp_PriEncounter(ii-1,1) >= Typ_Pri(Spark(ii+5))
            Deviation_E(ii-1,1) = 0;
        elseif MinTyp_PriEncounter(ii-1,1) < Typ_Pri(Spark(ii+5))
            Deviation_E(ii-1,1) = abs( MinTyp_PriEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        end
        %----------------------------------------------
        if (Max_PercentE(ii-1,1) >= 50) && (Deviation_E(ii-1,1) >= 100)
            if x_maxE < x_minE
                TypePosition_E(ii-1,1) = 1; % Max_Percent
            else
                TypePosition_E(ii-1,1) = 2; % Deviation
            end
            TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
            if Typ_Pri(Spark(ii+5)) < Typ_Pri(NDeltaY(ii-1,1))
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        elseif Max_PercentE(ii-1,1) >= 50
            TypePosition_E(ii-1,1) = 1; % Max_Percent
            TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
            if Typ_Pri(Spark(ii+5)) < Typ_Pri(NDeltaY(ii-1,1))
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        elseif Deviation_E(ii-1,1) >= 100
            TypePosition_E(ii-1,1) = 2; % Deviation
            TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
            if Typ_Pri(Spark(ii+5)) < Typ_Pri(NDeltaY(ii-1,1))
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        else
            TypePosition_E(ii-1,1) = 3; % No mode
            if Typ_Pri(Spark(ii+5)) < Typ_Pri(NDeltaY(ii-1,1))
                TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                    abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_MaxPer_Dev_E(ii-1,1) = -abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                    abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark(ii+5)) > PEncounter(ii-1,1) % Hunt Down
        Hunt_Encounter(ii-1,1) = {'Hunt_D'};
        if MinTyp_PriEncounter(ii-1,1) >= Typ_Pri(Spark(ii+5))
            Max_PercentE(ii-1,1) = 0;
        elseif ( Typ_Pri(Spark(ii+5)) > MinTyp_PriEncounter(ii-1,1) ) ...
                && ( MinTyp_PriEncounter(ii-1,1) > First_PEncounter(ii-1,1) )
            Max_PercentE(ii-1,1) = abs( MinTyp_PriEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        elseif MinTyp_PriEncounter(ii-1,1) <= First_PEncounter(ii-1,1)
            Max_PercentE(ii-1,1) = 100;
        end
        if MaxTyp_PriEncounter(ii-1,1) <= Typ_Pri(Spark(ii+5))
            Deviation_E(ii-1,1) = 0;
        elseif MaxTyp_PriEncounter(ii-1,1) > Typ_Pri(Spark(ii+5))
            Deviation_E(ii-1,1) = abs( MaxTyp_PriEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        end
        %----------------------------------------------
        if (Max_PercentE(ii-1,1) >= 50) && (Deviation_E(ii-1,1) >= 100)
            if x_maxE > x_minE
                TypePosition_E(ii-1,1) = 1; % Max_Percent
            else
                TypePosition_E(ii-1,1) = 2; % Deviation
            end
            TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
            if Typ_Pri(Spark(ii+5)) > Typ_Pri(NDeltaY(ii-1,1))
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        elseif Max_PercentE(ii-1,1) >= 50
            TypePosition_E(ii-1,1) = 1; % Max_Percent
            TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
            if Typ_Pri(Spark(ii+5)) > Typ_Pri(NDeltaY(ii-1,1))
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        elseif Deviation_E(ii-1,1) >= 100
            TypePosition_E(ii-1,1) = 2; % Deviation
            TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
            if Typ_Pri(Spark(ii+5)) > Typ_Pri(NDeltaY(ii-1,1))
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        else
            TypePosition_E(ii-1,1) = 3; % No mode
            if Typ_Pri(Spark(ii+5)) > Typ_Pri(NDeltaY(ii-1,1))
                TypePos_MaxPer_Dev_E(ii-1,1) = abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                    abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
                TypePos_E(ii-1,1) = {'Max_Percent'};
            else
                TypePos_MaxPer_Dev_E(ii-1,1) = -abs( Typ_Pri(NDeltaY(ii-1,1)) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                    abs( First_PEncounter(ii-1,1) - Typ_Pri(Spark(ii+5)) );
                TypePos_E(ii-1,1) = {'Deviation'};
            end
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark(ii+5)) == PEncounter(ii-1,1) % Hunt Equal
        Hunt_Encounter(ii-1,1) = {'Hunt_E'};
        Max_PercentE(ii-1,1) = -1;
        Deviation_E(ii-1,1) = -1;
        
        TypePosition_E(ii-1,1) = -1;
        TypePos_MaxPer_Dev_E(ii-1,1) = -1;
        TypePos_E(ii-1,1) = {'NO'};
    end
else
    First_encounter_CT(ii-1,1) = {'Infinity'};
    First_encounter_DT(ii-1,1) = -1;
    Second_encounter_CT(ii-1,1) = {'Infinity'};
    Second_encounter_DT(ii-1,1) = -1;
    Third_encounter_CT(ii-1,1) = {'Infinity'};
    Third_encounter_DT(ii-1,1) = -1;
    
    First_PEncounter(ii-1,1) = -1;
    Second_PEncounter(ii-1,1) = -1;
    Third_PEncounter(ii-1,1) = -1;
    
    XNUD_Encounter(ii-1,1) = {'NO'};
    Hunt_Encounter(ii-1,1) = {'NO'};
    
    MaxTyp_PriEncounter(ii-1,1) = -1;
    MinTyp_PriEncounter(ii-1,1) = -1;
    
    Max_PercentE(ii-1,1) = -1;
    Deviation_E(ii-1,1) = -1;
    
    TypePosition_E(ii-1,1) = -1;
    TypePos_MaxPer_Dev_E(ii-1,1) = -1;
    TypePos_E(ii-1,1) = {'NO'};
end

% Center of gravity
if NumerationE == 3
    Gravity_tE(ii-1,1) = ( NDeltaY(ii-1,1) + NDeltaY(ii-1,2) + NDeltaY(ii-1,3) ) / 3;
    Gravity_PE(ii-1,1) = ( PEncounter(ii-1,1) + PEncounter(ii-1,2) + PEncounter(ii-1,3) ) / 3;
else
    Gravity_tE(ii-1,1) = -1;
    Gravity_PE(ii-1,1) = -1;
end

if (NumerationE==0) || (NumerationE==1) || (NumerationE==2) % The # encounters is not enough.
    MaxTyp_PriceEG(ii-1,1) = -1;
    MinTyp_PriceEG(ii-1,1) = -1;
    
    XNUD_EG(ii-1,1) = {'NO'};
    Hunt_EG(ii-1,1) = {'NO'};
    
    Max_PercentEG(ii-1,1) = -1;
    Deviation_EG(ii-1,1) = -1;
elseif Gravity_tE(ii-1,1) > MM % Because there is no real data.
    MaxTyp_PriceEG(ii-1,1) = -1;
    MinTyp_PriceEG(ii-1,1) = -1;
    
    XNUD_EG(ii-1,1) = {'NO'};
    Hunt_EG(ii-1,1) = {'NO'};
    
    Max_PercentEG(ii-1,1) = -1;
    Deviation_EG(ii-1,1) = -1;
elseif Spark(ii+5) == 0 % The extreme point is not clear.
    MaxTyp_PriceEG(ii-1,1) = -1;
    MinTyp_PriceEG(ii-1,1) = -1;
    
    XNUD_EG(ii-1,1) = {'NO'};
    Hunt_EG(ii-1,1) = {'NO'};
    
    Max_PercentEG(ii-1,1) = -1;
    Deviation_EG(ii-1,1) = -1;
elseif Gravity_tE(ii-1,1) <= Spark(ii+5) % The first encounter is before the first spark.
    MaxTyp_PriceEG(ii-1,1) = -1;
    MinTyp_PriceEG(ii-1,1) = -1;
    
    XNUD_EG(ii-1,1) = {'NO'};
    Hunt_EG(ii-1,1) = {'NO'};
    
    Max_PercentEG(ii-1,1) = -1;
    Deviation_EG(ii-1,1) = -1;
elseif (Gravity_tE(ii-1,1) - Spark(ii+5)) < 1 % There is no data between first encounter and spark timing.
    MaxTyp_PriceEG(ii-1,1) = -1;
    MinTyp_PriceEG(ii-1,1) = -1;
    
    XNUD_EG(ii-1,1) = {'NO'};
    Hunt_EG(ii-1,1) = {'NO'};
    
    Max_PercentEG(ii-1,1) = -1;
    Deviation_EG(ii-1,1) = -1;
elseif (Gravity_tE(ii-1,1) - tt_T(ii+5) + 1) <= (t7_t1(ii-1,1) / 2)
    if Gravity_tE(ii-1,1) <= MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : Gravity_tE(ii-1,1) );
        MaxTyp_PriceEG(ii-1,1) = max(MTyp_Pri); % Find max
        MinTyp_PriceEG(ii-1,1) = min(MTyp_Pri); % Find min
        x_maxEG = find(MTyp_Pri == MaxTyp_PriceEG(ii-1,1)); % Additional
        x_minEG = find(MTyp_Pri == MinTyp_PriceEG(ii-1,1)); % Additional
    elseif Gravity_tE(ii-1,1) > MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : MM );
        MaxTyp_PriceEG(ii-1,1) = max(MTyp_Pri); % Find max
        MinTyp_PriceEG(ii-1,1) = min(MTyp_Pri); % Find min
        x_maxEG = find(MTyp_Pri == MaxTyp_PriceEG(ii-1,1)); % Additional
        x_minEG = find(MTyp_Pri == MinTyp_PriceEG(ii-1,1)); % Additional
    end
    
    if (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) < Gravity_PE(ii-1,1)) % b_T(7) is max
        XNUD_EG(ii-1,1) = {'XU'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) > Gravity_PE(ii-1,1)) % b_T(7) is max
        XNUD_EG(ii-1,1) = {'XD'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) < Gravity_PE(ii-1,1)) % b_T(7) is min
        XNUD_EG(ii-1,1) = {'NU'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) > Gravity_PE(ii-1,1)) % b_T(7) is min
        XNUD_EG(ii-1,1) = {'ND'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) == Gravity_PE(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_EG(ii-1,1) = {'XE'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) == Gravity_PE(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_EG(ii-1,1) = {'NE'};
    end
    
    if Typ_Pri(Spark(ii+5)) < Gravity_PE(ii-1,1) % Hunt Up
        Hunt_EG(ii-1,1) = {'Hunt_U'};
        if MaxTyp_PriceEG(ii-1,1) <= Typ_Pri(Spark(ii+5))
            Max_PercentEG(ii-1,1) = 0;
        elseif ( Typ_Pri(Spark(ii+5)) < MaxTyp_PriceEG(ii-1,1) ) ...
                && ( MaxTyp_PriceEG(ii-1,1) < Gravity_PE(ii-1,1) )
            Max_PercentEG(ii-1,1) = abs( MaxTyp_PriceEG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PE(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        elseif MaxTyp_PriceEG(ii-1,1) >= Gravity_PE(ii-1,1)
            Max_PercentEG(ii-1,1) = 100;
        end
        if MinTyp_PriceEG(ii-1,1) >= Typ_Pri(Spark(ii+5))
            Deviation_EG(ii-1,1) = 0;
        elseif MinTyp_PriceEG(ii-1,1) < Typ_Pri(Spark(ii+5))
            Deviation_EG(ii-1,1) = abs( MinTyp_PriceEG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PE(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        end
    elseif Typ_Pri(Spark(ii+5)) > Gravity_PE(ii-1,1) % Hunt Down
        Hunt_EG(ii-1,1) = {'Hunt_D'};
        if MinTyp_PriceEG(ii-1,1) >= Typ_Pri(Spark(ii+5))
            Max_PercentEG(ii-1,1) = 0;
        elseif ( Typ_Pri(Spark(ii+5)) > MinTyp_PriceEG(ii-1,1) ) ...
                && ( MinTyp_PriceEG(ii-1,1) > Gravity_PE(ii-1,1) )
            Max_PercentEG(ii-1,1) = abs( MinTyp_PriceEG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PE(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        elseif MinTyp_PriceEG(ii-1,1) <= Gravity_PE(ii-1,1)
            Max_PercentEG(ii-1,1) = 100;
        end
        if MaxTyp_PriceEG(ii-1,1) <= Typ_Pri(Spark(ii+5))
            Deviation_EG(ii-1,1) = 0;
        elseif MaxTyp_PriceEG(ii-1,1) > Typ_Pri(Spark(ii+5))
            Deviation_EG(ii-1,1) = abs( MaxTyp_PriceEG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PE(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        end
    elseif Typ_Pri(Spark(ii+5)) == Gravity_PE(ii-1,1) % Hunt Equal
        Hunt_EG(ii-1,1) = {'Hunt_E'};
        Max_PercentEG(ii-1,1) = -1;
        Deviation_EG(ii-1,1) = -1;
    end
else
    MaxTyp_PriceEG(ii-1,1) = -1;
    MinTyp_PriceEG(ii-1,1) = -1;
    
    XNUD_EG(ii-1,1) = {'NO'};
    Hunt_EG(ii-1,1) = {'NO'};
    
    Max_PercentEG(ii-1,1) = -1;
    Deviation_EG(ii-1,1) = -1;
end
%--------------------------------  With noise  ----------------------------
First_t = tt_T(ii+5);
Numeration = 0;
DeltaY3 = a1_Previous * First_t + a0_Previous ...
    + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * First_t + Sol_Total(ii-1,17)) ...
    + Sol_Total(ii-1,18) * sin(Sol_Total(ii-1,19) * First_t + Sol_Total(ii-1,20)) ...
    - a1 * First_t - a0 ...
    - Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26)) ...
    - Sol_Total(ii-1,27) * sin(Sol_Total(ii-1,28) * First_t + Sol_Total(ii-1,29));

while (Numeration < 3) && ( First_t <= (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) )
    First_t = First_t + 1;
    DeltaY4 = a1_Previous * First_t + a0_Previous ...
        + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * First_t + Sol_Total(ii-1,17)) ...
        + Sol_Total(ii-1,18) * sin(Sol_Total(ii-1,19) * First_t + Sol_Total(ii-1,20)) ...
        - a1 * First_t - a0 ...
        - Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26)) ...
        - Sol_Total(ii-1,27) * sin(Sol_Total(ii-1,28) * First_t + Sol_Total(ii-1,29));
    
    if ( sign(DeltaY3) ~= sign(DeltaY4)) && (First_t <= MM)
        Numeration = Numeration + 1;
        NY2(ii-1,Numeration) = Row(First_t,1);
        YForecast(ii-1,Numeration) = First_t;
        PForecast(ii-1,Numeration) = a1 * First_t + a0 + Sol_Total(ii-1,24) ...
            * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26)) + Sol_Total(ii-1,27) ...
            * sin(Sol_Total(ii-1,28) * First_t + Sol_Total(ii-1,29));
        
        if (Spark(ii+5) ~= 0) && (Typ_Pri(Spark(ii+5)) <= PForecast(ii-1,Numeration))
            PForecast_Rev(ii-1,Numeration) = Typ_Pri(Spark(ii+5)) - abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,Numeration));
        elseif (Spark(ii+5) ~= 0) && (Typ_Pri(Spark(ii+5)) > PForecast(ii-1,Numeration))
            PForecast_Rev(ii-1,Numeration) = Typ_Pri(Spark(ii+5)) + abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,Numeration));
        elseif Spark(ii+5) == 0
            PForecast_Rev(ii-1,Numeration) = 0;
        end
    elseif ( sign(DeltaY3) ~= sign(DeltaY4)) && (First_t > MM)
        Numeration = Numeration + 1;
        NY2(ii-1,Numeration) = {'Infinity'};
        YForecast(ii-1,Numeration) = First_t;
        PForecast(ii-1,Numeration) = a1 * First_t + a0 + Sol_Total(ii-1,24) ...
            * sin(Sol_Total(ii-1,25) * First_t + Sol_Total(ii-1,26)) + Sol_Total(ii-1,27) ...
            * sin(Sol_Total(ii-1,28) * First_t + Sol_Total(ii-1,29));
        
        if (Spark(ii+5) ~= 0) && (Typ_Pri(Spark(ii+5)) <= PForecast(ii-1,Numeration))
            PForecast_Rev(ii-1,Numeration) = Typ_Pri(Spark(ii+5)) - abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,Numeration));
        elseif (Spark(ii+5) ~= 0) && (Typ_Pri(Spark(ii+5)) > PForecast(ii-1,Numeration))
            PForecast_Rev(ii-1,Numeration) = Typ_Pri(Spark(ii+5)) + abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,Numeration));
        elseif Spark(ii+5) == 0
            PForecast_Rev(ii-1,Numeration) = 0;
        end
    end
    DeltaY3 = DeltaY4;
end

if Spark(ii+5) ~= 0
    Hunt_Price(ii-1,1) = Typ_Pri(Spark(ii+5));
else
    Hunt_Price(ii-1,1) = 0;
end

if (Numeration == 0) || (Numeration == 1) || (Numeration == 2) % The # encounters is not enough.
    First_Forecast_CT(ii-1,1) = {'Infinity'};
    First_Forecast_DT(ii-1,1) = -1;
    Second_Forecast_CT(ii-1,1) = {'Infinity'};
    Second_Forecast_DT(ii-1,1) = -1;
    Third_Forecast_CT(ii-1,1) = {'Infinity'};
    Third_Forecast_DT(ii-1,1) = -1;
    
    First_PForecast(ii-1,1) = -1;
    Second_PForecast(ii-1,1) = -1;
    Third_PForecast(ii-1,1) = -1;
    
    XNUD_Forecast(ii-1,1) = {'NO'};
    Hunt_Forecast(ii-1,1) = {'NO'};
    SellBuy(ii-1,1) = {'No'};
    
    MaxTyp_PriForecast(ii-1,1) = -1;
    MinTyp_PriForecast(ii-1,1) = -1;
    
%     Max_PercentF(ii-1,1) = -1;
%     Deviation_F(ii-1,1) = -1;
    
    TypePosition_F(ii-1,1) = -1;
    TypePos_MaxPer_Dev_F(ii-1,1) = -1;
    TypePos_F(ii-1,1) = {'NO'};
    
    FinalPrice_F(ii-1,1) = -1;
    PiP(ii-1,1) = -1;
    NetProfit(ii-1,1) = -1;
    PerPosition(ii-1,1) = -1;
    
    TimeOpen(ii-1,1) = -1;
    PriceOpen(ii-1,1) = -1;
    TimeClose(ii-1,1) = -1;
    PriceClose(ii-1,1) = -1;
    
    New_lineP(ii-1,1) = -1;
    New_lineP_D(ii-1,1) = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1; Cut_Number(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
    
elseif YForecast(ii-1,1) > MM % Because there is no real data.
    First_Forecast_CT(ii-1,1) = {'Infinity'};
    First_Forecast_DT(ii-1,1) = -1;
    Second_Forecast_CT(ii-1,1) = {'Infinity'};
    Second_Forecast_DT(ii-1,1) = -1;
    Third_Forecast_CT(ii-1,1) = {'Infinity'};
    Third_Forecast_DT(ii-1,1) = -1;
    
    First_PForecast(ii-1,1) = -1;
    Second_PForecast(ii-1,1) = -1;
    Third_PForecast(ii-1,1) = -1;
    
    XNUD_Forecast(ii-1,1) = {'NO'};
    Hunt_Forecast(ii-1,1) = {'NO'};
    SellBuy(ii-1,1) = {'No'};
    
    MaxTyp_PriForecast(ii-1,1) = -1;
    MinTyp_PriForecast(ii-1,1) = -1;
    
%     Max_PercentF(ii-1,1) = -1;
%     Deviation_F(ii-1,1) = -1;
    
    TypePosition_F(ii-1,1) = -1;
    TypePos_MaxPer_Dev_F(ii-1,1) = -1;
    TypePos_F(ii-1,1) = {'NO'};
    
    FinalPrice_F(ii-1,1) = -1;
    PiP(ii-1,1) = -1;
    NetProfit(ii-1,1) = -1;
    PerPosition(ii-1,1) = -1;
    
    TimeOpen(ii-1,1) = -1;
    PriceOpen(ii-1,1) = -1;
    TimeClose(ii-1,1) = -1;
    PriceClose(ii-1,1) = -1;
    
    New_lineP(ii-1,1) = -1;
    New_lineP_D(ii-1,1) = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1; Cut_Number(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
    
elseif Spark(ii+5) == 0 % The extreme point is not clear.
    First_Forecast_CT(ii-1,1) = {'Infinity'};
    First_Forecast_DT(ii-1,1) = -1;
    Second_Forecast_CT(ii-1,1) = {'Infinity'};
    Second_Forecast_DT(ii-1,1) = -1;
    Third_Forecast_CT(ii-1,1) = {'Infinity'};
    Third_Forecast_DT(ii-1,1) = -1;
    
    First_PForecast(ii-1,1) = -1;
    Second_PForecast(ii-1,1) = -1;
    Third_PForecast(ii-1,1) = -1;
    
    XNUD_Forecast(ii-1,1) = {'NO'};
    Hunt_Forecast(ii-1,1) = {'NO'};
    SellBuy(ii-1,1) = {'No'};
    
    MaxTyp_PriForecast(ii-1,1) = -1;
    MinTyp_PriForecast(ii-1,1) = -1;
    
%     Max_PercentF(ii-1,1) = -1;
%     Deviation_F(ii-1,1) = -1;
    
    TypePosition_F(ii-1,1) = -1;
    TypePos_MaxPer_Dev_F(ii-1,1) = -1;
    TypePos_F(ii-1,1) = {'NO'};
    
    FinalPrice_F(ii-1,1) = -1;
    PiP(ii-1,1) = -1;
    NetProfit(ii-1,1) = -1;
    PerPosition(ii-1,1) = -1;
    
    TimeOpen(ii-1,1) = -1;
    PriceOpen(ii-1,1) = -1;
    TimeClose(ii-1,1) = -1;
    PriceClose(ii-1,1) = -1;
    
    New_lineP(ii-1,1) = -1;
    New_lineP_D(ii-1,1) = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1; Cut_Number(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
    
elseif YForecast(ii-1,1) <= Spark(ii+5) % The first encounter is before the first spark.
    First_Forecast_CT(ii-1,1) = {'Infinity'};
    First_Forecast_DT(ii-1,1) = -1;
    Second_Forecast_CT(ii-1,1) = {'Infinity'};
    Second_Forecast_DT(ii-1,1) = -1;
    Third_Forecast_CT(ii-1,1) = {'Infinity'};
    Third_Forecast_DT(ii-1,1) = -1;
    
    First_PForecast(ii-1,1) = -1;
    Second_PForecast(ii-1,1) = -1;
    Third_PForecast(ii-1,1) = -1;
    
    XNUD_Forecast(ii-1,1) = {'NO'};
    Hunt_Forecast(ii-1,1) = {'NO'};
    SellBuy(ii-1,1) = {'No'};
    
    MaxTyp_PriForecast(ii-1,1) = -1;
    MinTyp_PriForecast(ii-1,1) = -1;
    
%     Max_PercentF(ii-1,1) = -1;
%     Deviation_F(ii-1,1) = -1;
    
    TypePosition_F(ii-1,1) = -1;
    TypePos_MaxPer_Dev_F(ii-1,1) = -1;
    TypePos_F(ii-1,1) = {'NO'};
    
    FinalPrice_F(ii-1,1) = -1;
    PiP(ii-1,1) = -1;
    NetProfit(ii-1,1) = -1;
    PerPosition(ii-1,1) = -1;
    
    TimeOpen(ii-1,1) = -1;
    PriceOpen(ii-1,1) = -1;
    TimeClose(ii-1,1) = -1;
    PriceClose(ii-1,1) = -1;
    
    New_lineP(ii-1,1) = -1;
    New_lineP_D(ii-1,1) = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1; Cut_Number(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
    
elseif (YForecast(ii-1,1) - Spark(ii+5)) < 1 % There is no data between first encounter and spark timing.
    First_Forecast_CT(ii-1,1) = {'Infinity'};
    First_Forecast_DT(ii-1,1) = -1;
    Second_Forecast_CT(ii-1,1) = {'Infinity'};
    Second_Forecast_DT(ii-1,1) = -1;
    Third_Forecast_CT(ii-1,1) = {'Infinity'};
    Third_Forecast_DT(ii-1,1) = -1;
    
    First_PForecast(ii-1,1) = -1;
    Second_PForecast(ii-1,1) = -1;
    Third_PForecast(ii-1,1) = -1;
    
    XNUD_Forecast(ii-1,1) = {'NO'};
    Hunt_Forecast(ii-1,1) = {'NO'};
    SellBuy(ii-1,1) = {'No'};
    
    MaxTyp_PriForecast(ii-1,1) = -1;
    MinTyp_PriForecast(ii-1,1) = -1;
    
%     Max_PercentF(ii-1,1) = -1;
%     Deviation_F(ii-1,1) = -1;
    
    TypePosition_F(ii-1,1) = -1;
    TypePos_MaxPer_Dev_F(ii-1,1) = -1;
    TypePos_F(ii-1,1) = {'NO'};
    
    FinalPrice_F(ii-1,1) = -1;
    PiP(ii-1,1) = -1;
    NetProfit(ii-1,1) = -1;
    PerPosition(ii-1,1) = -1;
    
    TimeOpen(ii-1,1) = -1;
    PriceOpen(ii-1,1) = -1;
    TimeClose(ii-1,1) = -1;
    PriceClose(ii-1,1) = -1;
    
    New_lineP(ii-1,1) = -1;
    New_lineP_D(ii-1,1) = -1;
    
    % RR=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1; Cut_Number(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
    
    % RR=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
    
    % RR=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
    
elseif (YForecast(ii-1,1) - tt_T(ii+5) + 1) <= (t7_t1(ii-1,1) / 2)
    First_Forecast_CT(ii-1,1) = NY2(ii-1,1);
    First_Forecast_DT(ii-1,1) = YForecast(ii-1,1);
    Second_Forecast_CT(ii-1,1) = NY2(ii-1,2);
    Second_Forecast_DT(ii-1,1) = YForecast(ii-1,2);
    Third_Forecast_CT(ii-1,1) = NY2(ii-1,3);
    Third_Forecast_DT(ii-1,1) = YForecast(ii-1,3);
    
    First_PForecast(ii-1,1) = PForecast(ii-1,1);
    Second_PForecast(ii-1,1) = PForecast(ii-1,2);
    Third_PForecast(ii-1,1) = PForecast(ii-1,3);
    
    if YForecast(ii-1,1) <= MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : YForecast(ii-1,1) );
        MaxTyp_PriForecast(ii-1,1) = max(MTyp_Pri); % Find max
        MinTyp_PriForecast(ii-1,1) = min(MTyp_Pri); % Find min
        x_maxF = find(MTyp_Pri == MaxTyp_PriForecast(ii-1,1)) + Spark(ii+5);
        x_minF = find(MTyp_Pri == MinTyp_PriForecast(ii-1,1)) + Spark(ii+5);
    elseif YForecast(ii-1,1) > MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : MM );
        MaxTyp_PriForecast(ii-1,1) = max(MTyp_Pri);
        MinTyp_PriForecast(ii-1,1) = min(MTyp_Pri);
        x_maxF = find(MTyp_Pri == MaxTyp_PriForecast(ii-1,1)) + Spark(ii+5);
        x_minF = find(MTyp_Pri == MinTyp_PriForecast(ii-1,1)) + Spark(ii+5);
    end
    
    if (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) < PForecast(ii-1,1)) % b_T(7) is max
        XNUD_Forecast(ii-1,1) = {'XU'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) > PForecast(ii-1,1)) % b_T(7) is max
        XNUD_Forecast(ii-1,1) = {'XD'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) < PForecast(ii-1,1)) % b_T(7) is min
        XNUD_Forecast(ii-1,1) = {'NU'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) > PForecast(ii-1,1)) % b_T(7) is min
        XNUD_Forecast(ii-1,1) = {'ND'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) == PForecast(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Forecast(ii-1,1) = {'XE'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) == PForecast(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_Forecast(ii-1,1) = {'NE'};
    end
    
    TimeOpen(ii-1,1) = Spark(ii+5);
    PriceOpen(ii-1,1) = Typ_Pri(Spark(ii+5));
    
    NumberTime = 0; % Position at the same time
    for w = 1 : ii - 2
        if TimeOpen(ii-1,1) < TimeClose(ii-1-w,1)
            NumberTime = NumberTime + 1;
        end
    end
    
    if NumberTime == 0
        if Typ_Pri(Spark(ii+5)) < First_PForecast(ii-1,1) % Hunt Up
            Hunt_Forecast(ii-1,1) = {'Hunt_U'};
            SellBuy(ii-1,1) = {'Buy'};
            TypePosition_F(ii-1,1) = -1; % No mode
            New_lineP(ii-1,1) = First_PForecast(ii-1,1) + 4 * 10^(-(Digit-1));
            New_lineP_D(ii-1,1) = Typ_Pri(Spark(ii+5)) - abs(First_PForecast(ii-1,1) - Typ_Pri(Spark(ii+5))) + 4 * 10^(-(Digit-1));
            
            for i = Spark(ii+5) : YForecast(ii-1,1) - 1
                if (((Typ_Pri(i) >= New_lineP(ii-1,1)) && (Typ_Pri(i+1) < New_lineP(ii-1,1))) ...
                        || ((Typ_Pri(i) <= New_lineP(ii-1,1)) && (Typ_Pri(i+1) > New_lineP(ii-1,1)))) ...
                        && (((Typ_Pri(i) <= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) > New_lineP_D(ii-1,1))) ...
                        || ((Typ_Pri(i) >= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) < New_lineP_D(ii-1,1))))
                    TypePosition_F(ii-1,1) = 2; % Deviation
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose(ii-1,1) = New_lineP_D(ii-1,1);
                    break
                elseif ((Typ_Pri(i) >= New_lineP(ii-1,1)) && (Typ_Pri(i+1) < New_lineP(ii-1,1))) ...
                        || ((Typ_Pri(i) <= New_lineP(ii-1,1)) && (Typ_Pri(i+1) > New_lineP(ii-1,1)))
                    TypePosition_F(ii-1,1) = 1; % Max_Percent
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 1;
                    TimeClose(ii-1,1) = i;
                    PriceClose(ii-1,1) = New_lineP(ii-1,1);
                    break
                elseif ((Typ_Pri(i) <= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) > New_lineP_D(ii-1,1))) ...
                        || ((Typ_Pri(i) >= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) < New_lineP_D(ii-1,1)))
                    TypePosition_F(ii-1,1) = 2; % Deviation
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose(ii-1,1) = New_lineP_D(ii-1,1);
                    break
                end
            end
            i = i + 1;
            if TypePosition_F(ii-1,1) == -1 % No mode
                if Typ_Pri(YForecast(ii-1,1)) >= Typ_Pri(Spark(ii+5)) + 4 * 10^(-(Digit-1))
                    TypePosition_F(ii-1,1) = 3; % No mode & in profit
                    TimeClose(ii-1,1) = YForecast(ii-1,1);
                    PriceClose(ii-1,1) = Typ_Pri(YForecast(ii-1,1));
                else
                    while (Typ_Pri(i) < Typ_Pri(Spark(ii+5)) + 4 * 10^(-(Digit-1))) ...
                            && (Typ_Pri(i) > New_lineP_D(ii-1,1)) && (i < MM)
                        i = i + 1;
                    end
                    
                    if i > MM
                        TypePosition_F(ii-1,1) = 100; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose(ii-1,1) = -1;
                    elseif Typ_Pri(i) >= Typ_Pri(Spark(ii+5)) + 4 * 10^(-(Digit-1))
                        TypePosition_F(ii-1,1) = 32; % No mode & Risk Free
                        TimeClose(ii-1,1) = i;
                        PriceClose(ii-1,1) = Typ_Pri(Spark(ii+5)) + 4 * 10^(-(Digit-1));
                    elseif Typ_Pri(i) <= New_lineP_D(ii-1,1)
                        TypePosition_F(ii-1,1) = 33; % No mode & Deviation
                        NSZO = NSZO + 1;
                        SequenceZeroOne(NSZO) = 0;
                        TimeClose(ii-1,1) = i;
                        PriceClose(ii-1,1) = New_lineP_D(ii-1,1);
                    else
                        TypePosition_F(ii-1,1) = 50; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose(ii-1,1) = -1;
                    end
                end
            end
            %----------------------------------------------
        elseif Typ_Pri(Spark(ii+5)) > First_PForecast(ii-1,1) % Hunt Down
            Hunt_Forecast(ii-1,1) = {'Hunt_D'};
            SellBuy(ii-1,1) = {'Sell'};
            TypePosition_F(ii-1,1) = -1; % No mode
            New_lineP(ii-1,1) = First_PForecast(ii-1,1) - 4 * 10^(-(Digit-1));
            New_lineP_D(ii-1,1) = Typ_Pri(Spark(ii+5)) + abs(First_PForecast(ii-1,1) - Typ_Pri(Spark(ii+5))) - 4 * 10^(-(Digit-1));
            
            for i = Spark(ii+5) : YForecast(ii-1,1) - 1
                if (((Typ_Pri(i) >= New_lineP(ii-1,1)) && (Typ_Pri(i+1) < New_lineP(ii-1,1))) ...
                        || ((Typ_Pri(i) <= New_lineP(ii-1,1)) && (Typ_Pri(i+1) > New_lineP(ii-1,1)))) ...
                        && (((Typ_Pri(i) <= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) > New_lineP_D(ii-1,1))) ...
                        || ((Typ_Pri(i) >= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) < New_lineP_D(ii-1,1))))
                    TypePosition_F(ii-1,1) = 2; % Deviation
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose(ii-1,1) = New_lineP_D(ii-1,1);
                    break
                elseif ((Typ_Pri(i) >= New_lineP(ii-1,1)) && (Typ_Pri(i+1) < New_lineP(ii-1,1))) ...
                        || ((Typ_Pri(i) <= New_lineP(ii-1,1)) && (Typ_Pri(i+1) > New_lineP(ii-1,1)))
                    TypePosition_F(ii-1,1) = 1; % Max_Percent
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 1;
                    TimeClose(ii-1,1) = i;
                    PriceClose(ii-1,1) = New_lineP(ii-1,1);
                    break
                elseif ((Typ_Pri(i) <= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) > New_lineP_D(ii-1,1))) ...
                        || ((Typ_Pri(i) >= New_lineP_D(ii-1,1)) && (Typ_Pri(i+1) < New_lineP_D(ii-1,1)))
                    TypePosition_F(ii-1,1) = 2; % Deviation
                    NSZO = NSZO + 1;
                    SequenceZeroOne(NSZO) = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose(ii-1,1) = New_lineP_D(ii-1,1);
                    break
                end
            end
            i = i + 1;
            if TypePosition_F(ii-1,1) == -1 % No mode
                if Typ_Pri(YForecast(ii-1,1)) <= Typ_Pri(Spark(ii+5)) - 4 * 10^(-(Digit-1))
                    TypePosition_F(ii-1,1) = 3; % No mode & in profit
                    TimeClose(ii-1,1) = YForecast(ii-1,1);
                    PriceClose(ii-1,1) = Typ_Pri(YForecast(ii-1,1));
                else
                    while (Typ_Pri(i) > Typ_Pri(Spark(ii+5)) - 4 * 10^(-(Digit-1))) ...
                            && (Typ_Pri(i) < New_lineP_D(ii-1,1)) && (i < MM)
                        i = i + 1;
                    end
                    
                    if i > MM
                        TypePosition_F(ii-1,1) = 100; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose(ii-1,1) = -1;
                    elseif Typ_Pri(i) <= Typ_Pri(Spark(ii+5)) - 4 * 10^(-(Digit-1))
                        TypePosition_F(ii-1,1) = 32; % No mode & Risk Free
                        TimeClose(ii-1,1) = i;
                        PriceClose(ii-1,1) = Typ_Pri(Spark(ii+5)) - 4 * 10^(-(Digit-1));
                    elseif Typ_Pri(i) >= New_lineP_D(ii-1,1)
                        TypePosition_F(ii-1,1) = 33; % No mode & Deviation
                        NSZO = NSZO + 1;
                        SequenceZeroOne(NSZO) = 0;
                        TimeClose(ii-1,1) = i;
                        PriceClose(ii-1,1) = New_lineP_D(ii-1,1);
                    else
                        TypePosition_F(ii-1,1) = 100; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose(ii-1,1) = -1;
                    end
                end
            end
            %----------------------------------------------
        elseif Typ_Pri(Spark(ii+5)) == First_PForecast(ii-1,1) % Hunt Equal
            Hunt_Forecast(ii-1,1) = {'Hunt_E'};
            SellBuy(ii-1,1) = {'No'};
            %         Max_PercentF(ii-1,1) = -1;
            %         Deviation_F(ii-1,1) = -1;
            
            TypePosition_F(ii-1,1) = -1;
            TypePos_MaxPer_Dev_F(ii-1,1) = -1;
            TypePos_F(ii-1,1) = {'NO'};
            
            FinalPrice_F(ii-1,1) = Typ_Pri(YForecast(ii-1,1));
            PiP(ii-1,1) = -1;
            NetProfit(ii-1,1) = -1;
            PerPosition(ii-1,1) = -1;
            
            TimeClose(ii-1,1) = -1;
            PriceClose(ii-1,1) = -1;
            
            New_lineP(ii-1,1) = -1;
            New_lineP_D(ii-1,1) = -1;
        end
    else
        Hunt_Forecast(ii-1,1) = {'Position at the same time'};
        SellBuy(ii-1,1) = {'Position at the same time'};
        % Max_PercentF(ii-1,1) = -1;
        % Deviation_F(ii-1,1) = -1;
        
        TypePosition_F(ii-1,1) = -1;
        TypePos_MaxPer_Dev_F(ii-1,1) = -1;
        TypePos_F(ii-1,1) = {'Position at the same time'};
        
        FinalPrice_F(ii-1,1) = Typ_Pri(YForecast(ii-1,1));
        PiP(ii-1,1) = -1;
        NetProfit(ii-1,1) = -1;
        PerPosition(ii-1,1) = -1;
        
        TimeClose(ii-1,1) = -1;
        PriceClose(ii-1,1) = -1;
        
        New_lineP(ii-1,1) = -1;
        New_lineP_D(ii-1,1) = -1;
    end
    %--------------------------  Select cut 1,2,3  ------------------------
    x1 = Sol_Total(ii-1,1) : YForecast(ii-1,3); % First Curves
    x_curve2 = Sol_Total(ii-1,3) : YForecast(ii-1,3); % Second Curves
    
    r1 = a1_Previous * x1 + a0_Previous; % First Curves
    y10 = r1 + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * x1 + Sol_Total(ii-1,17)) ...
        + Sol_Total(ii-1,18) * sin(Sol_Total(ii-1,19) * x1 + Sol_Total(ii-1,20)); % First Curves
    
    r2 = a1 * x_curve2 + a0; % Second Curves
    y11 = r2 + Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * x_curve2 + Sol_Total(ii-1,26)) ...
        + Sol_Total(ii-1,27) * sin(Sol_Total(ii-1,28) * x_curve2 + Sol_Total(ii-1,29)); % Second Curves
    
%     figure()
%     plot(x1,y10,'b.')
    S_E_1 = Select_Extreme(y10, Sol_Total(ii-1,1), YForecast(ii-1,3)); % First Curves
%     hold on
%     plot(x_curve2,y11,'g.')
    S_E_2 = Select_Extreme(y11, Sol_Total(ii-1,3), YForecast(ii-1,3)); % Second Curves
    
    if (S_E_1(1,1) ~= 0) && (S_E_2(1,1) ~= 0)
        % First Curves
        Ext1_Cut_3_time = S_E_1(size(S_E_1,1),1); % Cut 3
        %Ext1_Cut_3_Price = S_E_1(size(S_E_1,1),2);
        Ext1_Cut_3_MinMax = S_E_1(size(S_E_1,1),3);
        for i = size(S_E_1,1) : -1 : 1 % Cut 2
            if YForecast(ii-1,2) > S_E_1(i,1)
                Ext1_Cut_2_time = S_E_1(i,1);
                %Ext1_Cut_2_Price = S_E_1(i,2);
                Ext1_Cut_2_MinMax = S_E_1(i,3);
                break
            end
        end
        for i = size(S_E_1,1) : -1 : 1 % Cut 1
            if YForecast(ii-1,1) > S_E_1(i,1)
                Ext1_Cut_1_time = S_E_1(i,1);
                %Ext1_Cut_1_Price = S_E_1(i,2);
                Ext1_Cut_1_MinMax = S_E_1(i,3);
                break
            end
        end
        % Second Curves
        Ext2_Cut_3_time = S_E_2(size(S_E_2,1),1); % Cut 3
        %Ext2_Cut_3_Price = S_E_2(size(S_E_2,1),2);
        Ext2_Cut_3_MinMax = S_E_2(size(S_E_2,1),3);
        for i = size(S_E_2,1) : -1 : 1 % Cut 2
            if YForecast(ii-1,2) > S_E_2(i,1)
                Ext2_Cut_2_time = S_E_2(i,1);
                %Ext2_Cut_2_Price = S_E_2(i,2);
                Ext2_Cut_2_MinMax = S_E_2(i,3);
                break
            end
        end
        for i = size(S_E_2,1) : -1 : 1 % Cut 1
            if YForecast(ii-1,1) > S_E_2(i,1)
                Ext2_Cut_1_time = S_E_2(i,1);
                %Ext2_Cut_1_Price = S_E_2(i,2);
                Ext2_Cut_1_MinMax = S_E_2(i,3);
                break
            end
        end
        ii3 = ii3 + 1;
        Trend_Curves_Cut_1(ii3,1) = Result_Trend_Curves(Ext1_Cut_1_time, Ext2_Cut_1_time, ...
            Ext1_Cut_1_MinMax, Ext2_Cut_1_MinMax, YForecast(ii-1,1), S_E_1, S_E_2); % Cut 1
        
        Trend_Curves_Cut_2(ii3,1) = Result_Trend_Curves(Ext1_Cut_2_time, Ext2_Cut_2_time, ...
            Ext1_Cut_2_MinMax, Ext2_Cut_2_MinMax, YForecast(ii-1,2), S_E_1, S_E_2); % Cut 2
        
        Trend_Curves_Cut_3(ii3,1) = Result_Trend_Curves(Ext1_Cut_3_time, Ext2_Cut_3_time, ...
            Ext1_Cut_3_MinMax, Ext2_Cut_3_MinMax, YForecast(ii-1,3), S_E_1, S_E_2); % Cut 3
        
        %-----------------------  Three reward to risk  -------------------
        if (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==0)
            Cut_Number(ii-1,1) = -1;
            % R/R=1, TP=100, SL=100
            XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
            TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
            PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
            %  Reverse
            XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
            TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
            PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
            
            % R/R=1, TP=50, SL=50
            XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
            TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
            PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
            %  Reverse
            XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
            TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
            PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
            
            % R/R=3, TP=150, SL=50
            XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
            TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
            PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
            %  Reverse
            XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
            TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
            PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
            
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==0) % Cut 1
            if abs(PForecast(ii-1,1) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                Cut_Number(ii-1,1) = 1;
                % R/R=1, TP=100, SL=100
                [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                    TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                    PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                    NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                
                if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                    TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                    PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                    NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                
                if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                end
                
                % R/R=1, TP=50, SL=50
                [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                    TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                    PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                    NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                
                if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                    TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                    PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                    NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                
                if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                end
                
                % R/R=3, TP=150, SL=50
                [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                    TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                    PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                    NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                
                if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                    TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                    PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                    NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                
                if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                end
            else
                Cut_Number(ii-1,1) = -1;
                % R/R=1, TP=100, SL=100
                XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                
                % R/R=1, TP=50, SL=50
                XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                
                % R/R=3, TP=150, SL=50
                XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==0) % Cut 2
            if abs(PForecast(ii-1,2) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                Cut_Number(ii-1,1) = 2;
                % R/R=1, TP=100, SL=100
                [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                    TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                    PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                    NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                
                if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                    TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                    PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                    NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                
                if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                end
                
                % R/R=1, TP=50, SL=50
                [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                    TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                    PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                    NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                
                if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                    TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                    PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                    NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                
                if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                end
                
                % R/R=3, TP=150, SL=50
                [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                    TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                    PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                    NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                
                if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                    TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                    PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                    NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                
                if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                end
            else
                Cut_Number(ii-1,1) = -1;
                % R/R=1, TP=100, SL=100
                XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
                TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                
                % R/R=1, TP=50, SL=50
                XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                
                % R/R=3, TP=150, SL=50
                XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==1) % Cut 3
            if abs(PForecast(ii-1,3) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                Cut_Number(ii-1,1) = 3;
                % R/R=1, TP=100, SL=100
                [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                    TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                    PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                    NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                
                if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                    TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                    PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                    NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                
                if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                end
                
                % R/R=1, TP=50, SL=50
                [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                    TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                    PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                    NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                
                if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                    TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                    PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                    NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                
                if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                end
                
                % R/R=3, TP=150, SL=50
                [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                    TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                    PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                    NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                
                if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                    SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                end
                %  Reverse
                [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                    TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                    PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                    = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                    NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                
                if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                    SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                end
            else
                Cut_Number(ii-1,1) = -1;
                % R/R=1, TP=100, SL=100
                XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                
                % R/R=1, TP=50, SL=50
                XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                
                % R/R=3, TP=150, SL=50
                XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                %  Reverse
                XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==0)
            if abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) % Cut 1
                if abs(PForecast(ii-1,1) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 1;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            elseif abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) > abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) % Cut 2
                if abs(PForecast(ii-1,2) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 2;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==0) && (Trend_Curves_Cut_3(ii3,1)==1)
            if abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) % Cut 1
                if abs(PForecast(ii-1,1) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 1;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            elseif abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) > abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) % Cut 3
                if abs(PForecast(ii-1,3) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 3;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==0) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==1)
            if abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) % Cut 2
                if abs(PForecast(ii-1,2) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 2;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            elseif abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) > abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) % Cut 3
                if abs(PForecast(ii-1,3) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 3;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                   %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                   %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                   %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            end
        %------------------------------------------------------------------
        elseif (Trend_Curves_Cut_1(ii3,1)==1) && (Trend_Curves_Cut_2(ii3,1)==1) && (Trend_Curves_Cut_3(ii3,1)==1)
            if ( abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) ) && ...
                    ( abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) ) % Cut 1
                if abs(PForecast(ii-1,1) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 1;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast(ii-1,1), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,1), PForecast_Rev(ii-1,1), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            elseif ( abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) ) && ...
                    ( abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) ) % Cut 2
                if abs(PForecast(ii-1,2) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 2;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast(ii-1,2), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,2), PForecast_Rev(ii-1,2), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            elseif ( abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,1)) ) && ...
                    ( abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,3)) <= abs(Typ_Pri(Spark(ii+5)) - PForecast(ii-1,2)) ) % Cut 3
                if abs(PForecast(ii-1,3) - Typ_Pri(Spark(ii+5))) >= Dist_HC
                    Cut_Number(ii-1,1) = 3;
                    % R/R=1, TP=100, SL=100
                    [XNUD_Forecast_Sel_R1(ii-1,1), Hunt_Forecast_Sel_R1(ii-1,1), SellBuy_Sel_R1(ii-1,1), ...
                        TimeOpen_Sel_R1(ii-1,1), PriceOpen_Sel_R1(ii-1,1), TimeClose_Sel_R1, ...
                        PriceClose_Sel_R1(ii-1,1), TypePosition_Sel_R1(ii-1,1), NSZO_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R1, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 1;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    elseif (NSZO_R1 > 0) && (TypePosition_Sel_R1(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R1(NSZO_R1) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R1(ii-1,1), Hunt_Forecast_Rev_R1(ii-1,1), SellBuy_Rev_R1(ii-1,1), ...
                        TimeOpen_Rev_R1(ii-1,1), PriceOpen_Rev_R1(ii-1,1), TimeClose_Rev_R1, ...
                        PriceClose_Rev_R1(ii-1,1), TypePosition_Rev_R1(ii-1,1), NSZO_Rev_R1] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R1, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R1, ii, Digit, MM, 1, Stop_HD);
                    
                    if (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 1;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    elseif (NSZO_Rev_R1 > 0) && (TypePosition_Rev_R1(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R1(NSZO_Rev_R1) = 0;
                    end
                    
                    % R/R=1, TP=50, SL=50
                    [XNUD_Forecast_Sel_R2(ii-1,1), Hunt_Forecast_Sel_R2(ii-1,1), SellBuy_Sel_R2(ii-1,1), ...
                        TimeOpen_Sel_R2(ii-1,1), PriceOpen_Sel_R2(ii-1,1), TimeClose_Sel_R2, ...
                        PriceClose_Sel_R2(ii-1,1), TypePosition_Sel_R2(ii-1,1), NSZO_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R2, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 1;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    elseif (NSZO_R2 > 0) && (TypePosition_Sel_R2(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R2(NSZO_R2) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R2(ii-1,1), Hunt_Forecast_Rev_R2(ii-1,1), SellBuy_Rev_R2(ii-1,1), ...
                        TimeOpen_Rev_R2(ii-1,1), PriceOpen_Rev_R2(ii-1,1), TimeClose_Rev_R2, ...
                        PriceClose_Rev_R2(ii-1,1), TypePosition_Rev_R2(ii-1,1), NSZO_Rev_R2] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R2, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R2, ii, Digit, MM, 2, Stop_HD);
                    
                    if (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 1;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    elseif (NSZO_Rev_R2 > 0) && (TypePosition_Rev_R2(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R2(NSZO_Rev_R2) = 0;
                    end
                    
                    % R/R=3, TP=150, SL=50
                    [XNUD_Forecast_Sel_R3(ii-1,1), Hunt_Forecast_Sel_R3(ii-1,1), SellBuy_Sel_R3(ii-1,1), ...
                        TimeOpen_Sel_R3(ii-1,1), PriceOpen_Sel_R3(ii-1,1), TimeClose_Sel_R3, ...
                        PriceClose_Sel_R3(ii-1,1), TypePosition_Sel_R3(ii-1,1), NSZO_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast(ii-1,3), Spark(ii+5), ...
                        NSZO_R3, b_T(ii+4), b_T(ii+5), TimeClose_Sel_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 1)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 1;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 2)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    elseif (NSZO_R3 > 0) && (TypePosition_Sel_R3(ii-1,1) == 33)
                        SequenceZeroOne_Sel_R3(NSZO_R3) = 0;
                    end
                    %  Reverse
                    [XNUD_Forecast_Rev_R3(ii-1,1), Hunt_Forecast_Rev_R3(ii-1,1), SellBuy_Rev_R3(ii-1,1), ...
                        TimeOpen_Rev_R3(ii-1,1), PriceOpen_Rev_R3(ii-1,1), TimeClose_Rev_R3, ...
                        PriceClose_Rev_R3(ii-1,1), TypePosition_Rev_R3(ii-1,1), NSZO_Rev_R3] ...
                        = Select_Profit_S(Typ_Pri, YForecast(ii-1,3), PForecast_Rev(ii-1,3), Spark(ii+5), ...
                        NSZO_Rev_R3, b_T(ii+4), b_T(ii+5), TimeClose_Rev_R3, ii, Digit, MM, 3, Stop_HD);
                    
                    if (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 1)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 1;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 2)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    elseif (NSZO_Rev_R3 > 0) && (TypePosition_Rev_R3(ii-1,1) == 33)
                        SequenceZeroOne_Rev_R3(NSZO_Rev_R3) = 0;
                    end
                else
                    Cut_Number(ii-1,1) = -1;
                    % R/R=1, TP=100, SL=100
                    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'};
                    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
                    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
                    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
                    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
                    
                    % R/R=1, TP=50, SL=50
                    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
                    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
                    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
                    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
                    
                    % R/R=3, TP=150, SL=50
                    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
                    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
                    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
                    %  Reverse
                    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
                    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
                    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
                end
            end
        end
    else
        Cut_Number(ii-1,1) = -1;
        % R/R=1, TP=100, SL=100
        XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
        TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
        PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
        %  Reverse
        XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
        TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
        PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
        
        % R/R=1, TP=50, SL=50
        XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
        TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
        PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
        %  Reverse
        XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
        TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
        PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
        
        % R/R=3, TP=150, SL=50
        XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
        TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
        PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
        %  Reverse
        XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
        TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
        PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
    end
else
    First_Forecast_CT(ii-1,1) = {'Infinity'};
    First_Forecast_DT(ii-1,1) = -1;
    Second_Forecast_CT(ii-1,1) = {'Infinity'};
    Second_Forecast_DT(ii-1,1) = -1;
    Third_Forecast_CT(ii-1,1) = {'Infinity'};
    Third_Forecast_DT(ii-1,1) = -1;
    
    First_PForecast(ii-1,1) = -1;
    Second_PForecast(ii-1,1) = -1;
    Third_PForecast(ii-1,1) = -1;
    
    XNUD_Forecast(ii-1,1) = {'NO'};
    Hunt_Forecast(ii-1,1) = {'NO'};
    SellBuy(ii-1,1) = {'No'};
    
    MaxTyp_PriForecast(ii-1,1) = -1;
    MinTyp_PriForecast(ii-1,1) = -1;
    
%     Max_PercentF(ii-1,1) = -1;
%     Deviation_F(ii-1,1) = -1;
    
    TypePosition_F(ii-1,1) = -1;
    TypePos_MaxPer_Dev_F(ii-1,1) = -1;
    TypePos_F(ii-1,1) = {'NO'};
    
    FinalPrice_F(ii-1,1) = -1;
    PiP(ii-1,1) = -1;
    NetProfit(ii-1,1) = -1;
    PerPosition(ii-1,1) = -1;
    
    TimeOpen(ii-1,1) = -1;
    PriceOpen(ii-1,1) = -1;
    TimeClose(ii-1,1) = -1;
    PriceClose(ii-1,1) = -1;
    
    New_lineP(ii-1,1) = -1;
    New_lineP_D(ii-1,1) = -1;
    
    Cut_Number(ii-1,1)=-1;
    % R/R=1, TP=100, SL=100
    XNUD_Forecast_Sel_R1(ii-1,1)={'No'}; Hunt_Forecast_Sel_R1(ii-1,1)={'No'}; SellBuy_Sel_R1(ii-1,1)={'No'}; 
    TimeOpen_Sel_R1(ii-1,1)=-1; PriceOpen_Sel_R1(ii-1,1)=-1; TimeClose_Sel_R1(ii-1,1)=-1; 
    PriceClose_Sel_R1(ii-1,1)=-1; TypePosition_Sel_R1(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R1(ii-1,1)={'No'}; Hunt_Forecast_Rev_R1(ii-1,1)={'No'}; SellBuy_Rev_R1(ii-1,1)={'No'};
    TimeOpen_Rev_R1(ii-1,1)=-1; PriceOpen_Rev_R1(ii-1,1)=-1; TimeClose_Rev_R1(ii-1,1)=-1; 
    PriceClose_Rev_R1(ii-1,1)=-1; TypePosition_Rev_R1(ii-1,1)=-1;
    
    % R/R=1, TP=50, SL=50
    XNUD_Forecast_Sel_R2(ii-1,1)={'No'}; Hunt_Forecast_Sel_R2(ii-1,1)={'No'}; SellBuy_Sel_R2(ii-1,1)={'No'}; 
    TimeOpen_Sel_R2(ii-1,1)=-1; PriceOpen_Sel_R2(ii-1,1)=-1; TimeClose_Sel_R2(ii-1,1)=-1; 
    PriceClose_Sel_R2(ii-1,1)=-1; TypePosition_Sel_R2(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R2(ii-1,1)={'No'}; Hunt_Forecast_Rev_R2(ii-1,1)={'No'}; SellBuy_Rev_R2(ii-1,1)={'No'};
    TimeOpen_Rev_R2(ii-1,1)=-1; PriceOpen_Rev_R2(ii-1,1)=-1; TimeClose_Rev_R2(ii-1,1)=-1; 
    PriceClose_Rev_R2(ii-1,1)=-1; TypePosition_Rev_R2(ii-1,1)=-1;
    
    % R/R=3, TP=150, SL=50
    XNUD_Forecast_Sel_R3(ii-1,1)={'No'}; Hunt_Forecast_Sel_R3(ii-1,1)={'No'}; SellBuy_Sel_R3(ii-1,1)={'No'}; 
    TimeOpen_Sel_R3(ii-1,1)=-1; PriceOpen_Sel_R3(ii-1,1)=-1; TimeClose_Sel_R3(ii-1,1)=-1; 
    PriceClose_Sel_R3(ii-1,1)=-1; TypePosition_Sel_R3(ii-1,1)=-1;
    %  Reverse
    XNUD_Forecast_Rev_R3(ii-1,1)={'No'}; Hunt_Forecast_Rev_R3(ii-1,1)={'No'}; SellBuy_Rev_R3(ii-1,1)={'No'};
    TimeOpen_Rev_R3(ii-1,1)=-1; PriceOpen_Rev_R3(ii-1,1)=-1; TimeClose_Rev_R3(ii-1,1)=-1; 
    PriceClose_Rev_R3(ii-1,1)=-1; TypePosition_Rev_R3(ii-1,1)=-1;
end

% Center of gravity
if Numeration == 3
    Gravity_tF(ii-1,1) = ( YForecast(ii-1,1) + YForecast(ii-1,2) + YForecast(ii-1,3) ) / 3;
    Gravity_PF(ii-1,1) = ( PForecast(ii-1,1) + PForecast(ii-1,2) + PForecast(ii-1,3) ) / 3;
else
    Gravity_tF(ii-1,1) = -1;
    Gravity_PF(ii-1,1) = -1;
end

if (Numeration == 0) || (Numeration == 1) || (Numeration == 2) % The # encounters is not enough.
    MaxTyp_PriceFG(ii-1,1) = -1;
    MinTyp_PriceFG(ii-1,1) = -1;
    
    XNUD_FG(ii-1,1) = {'NO'};
    Hunt_FG(ii-1,1) = {'NO'};
    
    Max_PercentFG(ii-1,1) = -1;
    Deviation_FG(ii-1,1) = -1;
elseif Gravity_tF(ii-1,1) > MM % Because there is no real data.
    MaxTyp_PriceFG(ii-1,1) = -1;
    MinTyp_PriceFG(ii-1,1) = -1;
    
    XNUD_FG(ii-1,1) = {'NO'};
    Hunt_FG(ii-1,1) = {'NO'};
    
    Max_PercentFG(ii-1,1) = -1;
    Deviation_FG(ii-1,1) = -1;
elseif Spark(ii+5) == 0 % The extreme point is not clear.
    MaxTyp_PriceFG(ii-1,1) = -1;
    MinTyp_PriceFG(ii-1,1) = -1;
    
    XNUD_FG(ii-1,1) = {'NO'};
    Hunt_FG(ii-1,1) = {'NO'};
    
    Max_PercentFG(ii-1,1) = -1;
    Deviation_FG(ii-1,1) = -1;
elseif Gravity_tF(ii-1,1) <= Spark(ii+5) % The first encounter is before the first spark.
    MaxTyp_PriceFG(ii-1,1) = -1;
    MinTyp_PriceFG(ii-1,1) = -1;
    
    XNUD_FG(ii-1,1) = {'NO'};
    Hunt_FG(ii-1,1) = {'NO'};
    
    Max_PercentFG(ii-1,1) = -1;
    Deviation_FG(ii-1,1) = -1;
elseif (Gravity_tF(ii-1,1) - Spark(ii+5)) < 1 % There is no data between first encounter and spark timing.
    MaxTyp_PriceFG(ii-1,1) = -1;
    MinTyp_PriceFG(ii-1,1) = -1;
    
    XNUD_FG(ii-1,1) = {'NO'};
    Hunt_FG(ii-1,1) = {'NO'};
    
    Max_PercentFG(ii-1,1) = -1;
    Deviation_FG(ii-1,1) = -1;
elseif (Gravity_tF(ii-1,1) - tt_T(ii+5) + 1) <= (t7_t1(ii-1,1) / 2)
    if Gravity_tF(ii-1,1) <= MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : Gravity_tF(ii-1,1) );
        MaxTyp_PriceFG(ii-1,1) = max(MTyp_Pri); % Find max
        MinTyp_PriceFG(ii-1,1) = min(MTyp_Pri); % Find min
        x_maxFG = find(MTyp_Pri == MaxTyp_PriceFG(ii-1,1)); % Additional
        x_minFG = find(MTyp_Pri == MinTyp_PriceFG(ii-1,1)); % Additional
    elseif Gravity_tF(ii-1,1) > MM
        MTyp_Pri = Typ_Pri( Spark(ii+5)+1 : MM );
        MaxTyp_PriceFG(ii-1,1) = max(MTyp_Pri);
        MinTyp_PriceFG(ii-1,1) = min(MTyp_Pri);
        x_maxFG = find(MTyp_Pri == MaxTyp_PriceFG(ii-1,1)); % Additional
        x_minFG = find(MTyp_Pri == MinTyp_PriceFG(ii-1,1)); % Additional
    end
    
    if (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) < Gravity_PF(ii-1,1)) % b_T(7) is max
        XNUD_FG(ii-1,1) = {'XU'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) > Gravity_PF(ii-1,1)) % b_T(7) is max
        XNUD_FG(ii-1,1) = {'XD'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) < Gravity_PF(ii-1,1)) % b_T(7) is min
        XNUD_FG(ii-1,1) = {'NU'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) > Gravity_PF(ii-1,1)) % b_T(7) is min
        XNUD_FG(ii-1,1) = {'ND'};
    elseif (b_T(ii+4) < b_T(ii+5)) && (b_T(ii+5) == Gravity_PF(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_FG(ii-1,1) = {'XE'};
    elseif (b_T(ii+4) > b_T(ii+5)) && (b_T(ii+5) == Gravity_PF(ii-1,1)) % b_T(7) is equal to PForecast(1,1).
        XNUD_FG(ii-1,1) = {'NE'};
    end
    
    if Typ_Pri(Spark(ii+5)) < Gravity_PF(ii-1,1) % Hunt Up
        Hunt_FG(ii-1,1) = {'Hunt_U'};
        if MaxTyp_PriceFG(ii-1,1) <= Typ_Pri(Spark(ii+5))
            Max_PercentFG(ii-1,1) = 0;
        elseif ( Typ_Pri(Spark(ii+5)) < MaxTyp_PriceFG(ii-1,1) ) ...
                && ( MaxTyp_PriceFG(ii-1,1) < Gravity_PF(ii-1,1) )
            Max_PercentFG(ii-1,1) = abs( MaxTyp_PriceFG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PF(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        elseif MaxTyp_PriceFG(ii-1,1) >= Gravity_PF(ii-1,1)
            Max_PercentFG(ii-1,1) = 100;
        end
        if MinTyp_PriceFG(ii-1,1) >= Typ_Pri(Spark(ii+5))
            Deviation_FG(ii-1,1) = 0;
        elseif MinTyp_PriceFG(ii-1,1) < Typ_Pri(Spark(ii+5))
            Deviation_FG(ii-1,1) = abs( MinTyp_PriceFG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PF(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        end
    elseif Typ_Pri(Spark(ii+5)) > Gravity_PF(ii-1,1) % Hunt Down
        Hunt_FG(ii-1,1) = {'Hunt_D'};
        if MinTyp_PriceFG(ii-1,1) >= Typ_Pri(Spark(ii+5))
            Max_PercentFG(ii-1,1) = 0;
        elseif ( Typ_Pri(Spark(ii+5)) > MinTyp_PriceFG(ii-1,1) ) ...
                && ( MinTyp_PriceFG(ii-1,1) > Gravity_PF(ii-1,1) )
            Max_PercentFG(ii-1,1) = abs( MinTyp_PriceFG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PF(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        elseif MinTyp_PriceFG(ii-1,1) <= Gravity_PF(ii-1,1)
            Max_PercentFG(ii-1,1) = 100;
        end
        if MaxTyp_PriceFG(ii-1,1) <= Typ_Pri(Spark(ii+5))
            Deviation_FG(ii-1,1) = 0;
        elseif MaxTyp_PriceFG(ii-1,1) > Typ_Pri(Spark(ii+5))
            Deviation_FG(ii-1,1) = abs( MaxTyp_PriceFG(ii-1,1) - Typ_Pri(Spark(ii+5)) ) * 100 / ...
                abs( Gravity_PF(ii-1,1) - Typ_Pri(Spark(ii+5)) );
        end
    elseif Typ_Pri(Spark(ii+5)) == Gravity_PF(ii-1,1) % Hunt Equal
        Hunt_FG(ii-1,1) = {'Hunt_E'};
        Max_PercentFG(ii-1,1) = -1;
        Deviation_FG(ii-1,1) = -1;
    end
else
    MaxTyp_PriceFG(ii-1,1) = -1;
    MinTyp_PriceFG(ii-1,1) = -1;
    
    XNUD_FG(ii-1,1) = {'NO'};
    Hunt_FG(ii-1,1) = {'NO'};
    
    Max_PercentFG(ii-1,1) = -1;
    Deviation_FG(ii-1,1) = -1;
end
%---------------------------------  No noise  -----------------------------
% if NumerationE == 3
% x1 = Sol_Total(ii-1,1) : NDeltaY(ii-1,3);
% x_curve2 = Sol_Total(ii-1,3) : NDeltaY(ii-1,3);
% 
% r1 = a1_Previous * x1 + a0_Previous;
% y6 = r1 + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * x1 + Sol_Total(ii-1,17));
% 
% r2 = a1 * x_curve2 + a0;
% y7 = r2 + Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * x_curve2 + Sol_Total(ii-1,26));
% if NDeltaY(ii-1,3) > MM
%     figure()
%     p1 = plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
%     hold on
%     p2 = plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
%     hold on
% elseif Second_encounter_DT(ii-1,1) > MM
%     figure()
%     p1 = plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
%     hold on
%     p2 = plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
%     hold on
% elseif First_encounter_DT(ii-1,1) > MM
%     figure()
%     p1 = plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
%     hold on
%     p2 = plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
%     hold on
% else
%     figure()
%     p1 = plot(x1,Typ_Pri(Sol_Total(ii-1,1):NDeltaY(ii-1,3)),'k.');
%     hold on
%     bi = ii + 5;
%     while ( bi < length(tt_T) ) && ( tt_T(bi) < NDeltaY(ii-1,3) )
%         bi = bi + 1;
%     end
%     p2 = plot(tt_T(ii-1:bi),b_T(ii-1:bi),'mo','linewidth',2);
%     hold on
% end
% p3 = plot(x1,y6,'b');
% hold on
% p4 = plot(x_curve2,y7,'g');
% title('Systems 1 and 2: r+a_1*sin(w_1t+\theta_1)')
% %legend('Real Data','Six Points','First System','Second System')
% %axis tight
% Y_1 = ylim;
% 
% for i = 1 : 3 % Additional
%     hold on
%     X_1 = [NDeltaY(ii-1,i), NDeltaY(ii-1,i)];
%     p5 = line(X_1,Y_1,'Color',[.5 .3 1.0]);%[.6 .1 .9]
%     hold on
%     p6 = plot(NDeltaY(ii-1,i),PEncounter(ii-1,i),'ch','linewidth',2); % Additional
% end
% hold on
% p7 = plot(x_maxE,MaxTyp_PriEncounter(ii-1,1),'gp','linewidth',1.5); % Additional
% hold on
% plot(x_minE,MinTyp_PriEncounter(ii-1,1),'gp','linewidth',1.5) % Additional
% hold on
% X_1 = [Spark(ii+5), Spark(ii+5)];
% p8 = line(X_1,Y_1,'Color',[.3 .9 .9]); % Additional
% hold on
% p9 = plot(Gravity_tE(ii-1,1),Gravity_PE(ii-1,1),'r>','linewidth',1.5);%p9.Marker = '>';
% hold on
% X_1 = [Gravity_tE(ii-1,1), Gravity_tE(ii-1,1)];
% p10 = line(X_1,Y_1,'Color',[1.0 .2 1.0]);%[.9 .1 .9]
% p10.LineStyle = '--';
% hold on
% p11 = plot(x_maxEG+Spark(ii+5),MaxTyp_PriceEG(ii-1,1),'g^','linewidth',1.5); % Additional
% hold on
% plot(x_minEG+Spark(ii+5),MinTyp_PriceEG(ii-1,1),'g^','linewidth',1.5); % Additional
% legend([p1,p2,p3,p4,p5,p6,p7(1,1),p8,p9,p10,p11(1,1)],'Real Data','Six Points','First System',...
%     'Second System','Encounter of Two System','Encounter of Two System',...
%     'Min & Max Real Data','Hunt','Center of gravity','Center of gravity',...
%     'Min & Max Real Data (gravity)')
% hold on
% if Spark(ii+5) ~= 0
%     plot(Spark(ii+5),Typ_Pri(Spark(ii+5)),'rh','linewidth',2) % Additional
% end
% 
% end

%---------------------------  Life time (no noise)  -----------------------
%x_curve1 = Sol_Total(1,11) : NDeltaY(1,3);
x_curve2 = Spark(ii+5) + 1 : (tt_T(ii+5) + 2 * t7_t1(ii-1,1));

y16 =@(t_1) a1_Previous * t_1 + a0_Previous ...
    + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * t_1 + Sol_Total(ii-1,17));

y17 =@(t_1) a1 * t_1 + a0 ...
    + Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * t_1 + Sol_Total(ii-1,26));

dis1(1) = abs( y16(Sol_Total(ii-1,1)) - Sol_Total(ii-1,2) );
dis1(2) = abs( y16(Sol_Total(ii-1,3)) - Sol_Total(ii-1,4) );
dis1(3) = abs( y16(Sol_Total(ii-1,5)) - Sol_Total(ii-1,6) );
dis1(4) = abs( y16(Sol_Total(ii-1,7)) - Sol_Total(ii-1,8) );
dis1(5) = abs( y16(Sol_Total(ii-1,9)) - Sol_Total(ii-1,10) );
dis1(6) = abs( y16(Sol_Total(ii-1,11)) - Sol_Total(ii-1,12) );
max_dis1 = max(dis1);
average1 =  ( dis1(1) + dis1(2) + dis1(3) + dis1(4) + dis1(5) + dis1(6) ) / 6;

dis2(1) = abs( y17(Sol_Total(ii-1,3)) - Sol_Total(ii-1,4) );
dis2(2) = abs( y17(Sol_Total(ii-1,5)) - Sol_Total(ii-1,6) );
dis2(3) = abs( y17(Sol_Total(ii-1,7)) - Sol_Total(ii-1,8) );
dis2(4) = abs( y17(Sol_Total(ii-1,9)) - Sol_Total(ii-1,10) );
dis2(5) = abs( y17(Sol_Total(ii-1,11)) - Sol_Total(ii-1,12) );
dis2(6) = abs( y17(Sol_Total(ii-1,13)) - Sol_Total(ii-1,14) );
max_dis2 = max(dis2);
average2 =  ( dis2(1) + dis2(2) + dis2(3) + dis2(4) + dis2(5) + dis2(6) ) / 6;

y20 = y16(x_curve2) + (max_dis1 + average1) / 2;
y21 = y16(x_curve2) - (max_dis1 + average1) / 2;

y22 = y17(x_curve2) + (max_dis2 + average2) / 2;
y23 = y17(x_curve2) - (max_dis2 + average2) / 2;

LifeTime1 = -1; LifeTime2 = -1;
if (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) <= MM
    for i = Spark(ii+5) + 1 : (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) - 1
        if ( Typ_Pri(i+1) >= y20(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y20(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y20(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y20(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y21(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y21(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y21(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y21(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        end
    end
else
    for i = Spark(ii+5) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y20(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y20(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y20(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y20(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y21(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y21(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y21(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y21(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        end
    end
end
if (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) <= MM
    for i = Spark(ii+5) + 1 : (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) - 1
        if ( Typ_Pri(i+1) >= y22(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y22(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y22(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y22(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y23(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y23(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y23(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y23(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        end
    end
else
    for i = Spark(ii+5) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y22(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y22(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y22(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y22(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y23(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y23(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y23(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y23(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        end
    end
end
if (LifeTime1 == -1) && (LifeTime2 == -1)
    LifeTimeA(ii-1,1) = -1;
elseif (LifeTime1 ~= -1) && (LifeTime2 == -1)
    LifeTimeA(ii-1,1) = LifeTime1; % Life time for no nise
elseif (LifeTime1 == -1) && (LifeTime2 ~= -1)
    LifeTimeA(ii-1,1) = LifeTime2; % Life time for no nise
elseif (LifeTime1 ~= -1) && (LifeTime2 ~= -1)
    LifeTimeA(ii-1,1) = max(LifeTime1,LifeTime2); % Life time for no nise
else
    LifeTimeA(ii-1,1) = -1;
end

% hold on
% plot(x_curve2,y20,'b-.',x_curve2,y21,'b-.');
% hold on
% plot(x_curve2,y22,'g-.',x_curve2,y23,'g-.');
% hold on
% if LifeTime1 ~= -1
%     p12 = line(LifeTime1,Typ_Pri(LifeTime1),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p12.Marker = 's';
% end
% hold on
% if LifeTime2 ~= -1
%     p13 = line(LifeTime2,Typ_Pri(LifeTime2),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p13.Marker = 's';
% end

%--------------------------------  With noise  ----------------------------

if (Numeration == 3) && ( (TypePosition_Sel_R2(ii-1,1) == 1) || ...
        (TypePosition_Sel_R2(ii-1,1) == 2) || (TypePosition_Sel_R2(ii-1,1) == 33) )
x1 = Sol_Total(ii-1,1) : YForecast(ii-1,3);
x_curve2 = Sol_Total(ii-1,3) : YForecast(ii-1,3);

r1 = a1_Previous * x1 + a0_Previous;
y10 = r1 + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * x1 + Sol_Total(ii-1,17)) ...
    + Sol_Total(ii-1,18) * sin(Sol_Total(ii-1,19) * x1 + Sol_Total(ii-1,20));

r2 = a1 * x_curve2 + a0;
y11 = r2 + Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * x_curve2 + Sol_Total(ii-1,26)) ...
    + Sol_Total(ii-1,27) * sin(Sol_Total(ii-1,28) * x_curve2 + Sol_Total(ii-1,29));

if YForecast(ii-1,3) > MM
    figure()
    p1 = plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
    hold on
    p2 = plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
    hold on
elseif Second_Forecast_DT(ii-1,1) > MM
    figure()
    p1 = plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
    hold on
    p2 = plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
    hold on
elseif First_Forecast_DT(ii-1,1) > MM
    figure()
    p1 = plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
    hold on
    p2 = plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
    hold on
else
    figure()
    p1 = plot(x1,Typ_Pri(Sol_Total(ii-1,1):YForecast(ii-1,3)),'k.');
    hold on
    bi = ii + 5;
    while ( bi < length(tt_T) ) && ( tt_T(bi) < YForecast(ii-1,3) )
        bi = bi + 1;
    end
    p2 = plot(tt_T(ii-1:bi),b_T(ii-1:bi),'mo','linewidth',2);
    hold on
end
p3 = plot(x1,y10,'b.');
hold on
p4 = plot(x_curve2,y11,'g.');
title('Systems 1 and 2: r+a_1*sin(w_1t+\theta_1)+a_2*sin(w_2t+\theta_2)')
%legend('Real Data','Six Points','First System','Second System')
%axis tight
Y_1 = ylim;

for i = 1 : 3 % Additional
    hold on
    X_1 = [YForecast(ii-1,i), YForecast(ii-1,i)];
    p5 = line(X_1,Y_1,'Color',[.5 .3 1.0]);%[.6 .1 .9]
    hold on
    p6 = plot(YForecast(ii-1,i),PForecast(ii-1,i),'ch','linewidth',2); % Additional
end
hold on
p7 = plot(x_maxF,MaxTyp_PriForecast(ii-1,1),'gp','linewidth',1.5); % Additional
plot(x_minF,MinTyp_PriForecast(ii-1,1),'gp','linewidth',1.5) % Additional
hold on
X_1 = [Spark(ii+5), Spark(ii+5)];
p8 = line(X_1,Y_1,'Color',[.3 .9 .9]); % Additional
hold on
% p9 = plot(Gravity_tF(ii-1,1),Gravity_PF(ii-1,1),'r>','linewidth',1.5);%p9.Marker = '>';
% hold on
% X_1 = [Gravity_tF(ii-1,1), Gravity_tF(ii-1,1)];
% p10 = line(X_1,Y_1,'Color',[1.0 .2 1.0]);%[.9 .1 .9]
% p10.LineStyle = '--';
% hold on
% p11 = plot(x_maxFG+Spark(ii+5),MaxTyp_PriceFG(ii-1,1),'g^','linewidth',1.5); % Additional
% hold on
%plot(x_minFG+Spark(ii+5),MinTyp_PriceFG(ii-1,1),'g^','linewidth',1.5); % Additional
legend([p1,p2,p3,p4,p5,p6,p7(1,1),p8],'Real Data','Six Points','First System',...
    'Second System','Encounter of Two System','Encounter of Two System',...
    'Min & Max Real Data','Hunt')
% legend([p1,p2,p3,p4,p5,p6,p7(1,1),p8,p9,p10,p11(1,1)],'Real Data','Six Points','First System',...
%     'Second System','Encounter of Two System','Encounter of Two System',...
%     'Min & Max Real Data','Hunt','Center of gravity','Center of gravity',...
%     'Min & Max Real Data (gravity)')
hold on
if Spark(ii+5) ~= 0
    plot(Spark(ii+5),Typ_Pri(Spark(ii+5)),'rh','linewidth',2) % Additional
end

for i = 1 : 3 % Additional
    hold on
    plot(YForecast(ii-1,i),PForecast_Rev(ii-1,i),'cs','linewidth',2); % Additional
end

if (Spark(ii+5) ~= 0) && (TypePosition_Sel_R2(ii-1,1) ~= -1) && (Typ_Pri(Spark(ii+5)) < PForecast(ii-1,Cut_Number(ii-1,1)))
    X_2 = xlim;
    New_lineP(ii-1,1) = (PForecast(ii-1,Cut_Number(ii-1,1)) + Typ_Pri(Spark(ii+5))) / 2  + 4 * 10^(-(Digit-1)); % TP = 50
    New_lineP_D(ii-1,1) = Typ_Pri(Spark(ii+5)) - abs(PForecast(ii-1,Cut_Number(ii-1,1)) - Typ_Pri(Spark(ii+5))) / 2 + 4 * 10^(-(Digit-1)); % SL = 50
    
    Y_2 = [New_lineP(ii-1,1), New_lineP(ii-1,1)];
    Y_3 = [New_lineP_D(ii-1,1), New_lineP_D(ii-1,1)];
    Y_4 = [Hunt_Price(ii-1,1) + 4 * 10^(-(Digit-1)), Hunt_Price(ii-1,1) + 4 * 10^(-(Digit-1))];
    hold on
    line(X_2,Y_2,'Color',[.5 .3 1.0]); % Hunt Up
    hold on
    line(X_2,Y_3,'Color',[1.0 .2 1.0]);
    hold on
    line(X_2,Y_4,'Color',[.3 .9 .9]);
    hold on
    plot(TimeClose_Sel_R2(ii-1,1),PriceClose_Sel_R2(ii-1,1),'b^','linewidth',1.5); % Additional
elseif (Spark(ii+5) ~= 0) && (TypePosition_Sel_R2(ii-1,1) ~= -1) && (Typ_Pri(Spark(ii+5)) > PForecast(ii-1,Cut_Number(ii-1,1)))
    X_2 = xlim;
    New_lineP(ii-1,1) = (PForecast(ii-1,Cut_Number(ii-1,1)) + Typ_Pri(Spark(ii+5))) / 2  - 4 * 10^(-(Digit-1)); % TP = 50
    New_lineP_D(ii-1,1) = Typ_Pri(Spark(ii+5)) + abs(PForecast(ii-1,Cut_Number(ii-1,1)) - Typ_Pri(Spark(ii+5))) / 2 - 4 * 10^(-(Digit-1)); % SL = 50
    
    Y_2 = [New_lineP(ii-1,1), New_lineP(ii-1,1)];
    Y_3 = [New_lineP_D(ii-1,1), New_lineP_D(ii-1,1)];
    Y_4 = [Hunt_Price(ii-1,1) - 4 * 10^(-(Digit-1)), Hunt_Price(ii-1,1) - 4 * 10^(-(Digit-1))];
    hold on
    line(X_2,Y_2,'Color',[.5 .3 1.0]); % Hunt Down
    hold on
    line(X_2,Y_3,'Color',[1.0 .2 1.0]);
    hold on
    line(X_2,Y_4,'Color',[.3 .9 .9]);
    hold on
    plot(TimeClose_Sel_R2(ii-1,1),PriceClose_Sel_R2(ii-1,1),'b^','linewidth',1.5); % Additional
elseif (TypePosition_Sel_R2(ii-1,1) == 100)
    X_2 = xlim;
    New_lineP(ii-1,1) = (PForecast(ii-1,Cut_Number(ii-1,1)) + Typ_Pri(Spark(ii+5))) / 2  - 4 * 10^(-(Digit-1)); % TP = 50
    New_lineP_D(ii-1,1) = Typ_Pri(Spark(ii+5)) + abs(PForecast(ii-1,Cut_Number(ii-1,1)) - Typ_Pri(Spark(ii+5))) / 2 - 4 * 10^(-(Digit-1)); % SL = 50
    
    Y_2 = [New_lineP(ii-1,1), New_lineP(ii-1,1)];
    Y_3 = [New_lineP_D(ii-1,1), New_lineP_D(ii-1,1)];
    Y_4 = [Hunt_Price(ii-1,1) - 4 * 10^(-(Digit-1)), Hunt_Price(ii-1,1) - 4 * 10^(-(Digit-1))];
    hold on
    line(X_2,Y_2,'Color',[.5 .3 1.0]); % Hunt Down
    hold on
    line(X_2,Y_3,'Color',[1.0 .2 1.0]);
    hold on
    line(X_2,Y_4,'Color',[.3 .9 .9]);
    hold on
    plot(TimeClose_Sel_R2(ii-1,1),PriceClose_Sel_R2(ii-1,1),'b^','linewidth',1.5); % Additional
end

end

%---------------------------  Life time (with noise)  -----------------------
%x_curve1 = Sol_Total(1,11) : YForecast(1,3);
x_curve2 = Spark(ii+5) + 1 : (tt_T(ii+5) + 2 * t7_t1(ii-1,1));

y18 =@(t_1) a1_Previous * t_1 + a0_Previous ...
    + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * t_1 + Sol_Total(ii-1,17)) ...
    + Sol_Total(ii-1,18) * sin(Sol_Total(ii-1,19) * t_1 + Sol_Total(ii-1,20));

y19 =@(t_1) a1 * t_1 + a0 ...
    + Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * t_1 + Sol_Total(ii-1,26))...
    + Sol_Total(ii-1,27) * sin(Sol_Total(ii-1,28) * t_1 + Sol_Total(ii-1,29));

dis1(1) = abs( y18(Sol_Total(ii-1,1)) - Sol_Total(ii-1,2) );
dis1(2) = abs( y18(Sol_Total(ii-1,3)) - Sol_Total(ii-1,4) );
dis1(3) = abs( y18(Sol_Total(ii-1,5)) - Sol_Total(ii-1,6) );
dis1(4) = abs( y18(Sol_Total(ii-1,7)) - Sol_Total(ii-1,8) );
dis1(5) = abs( y18(Sol_Total(ii-1,9)) - Sol_Total(ii-1,10) );
dis1(6) = abs( y18(Sol_Total(ii-1,11)) - Sol_Total(ii-1,12) );
max_dis1 = max(dis1);
average1 =  ( dis1(1) + dis1(2) + dis1(3) + dis1(4) + dis1(5) + dis1(6) ) / 6;

dis2(1) = abs( y19(Sol_Total(ii-1,3)) - Sol_Total(ii-1,4) );
dis2(2) = abs( y19(Sol_Total(ii-1,5)) - Sol_Total(ii-1,6) );
dis2(3) = abs( y19(Sol_Total(ii-1,7)) - Sol_Total(ii-1,8) );
dis2(4) = abs( y19(Sol_Total(ii-1,9)) - Sol_Total(ii-1,10) );
dis2(5) = abs( y19(Sol_Total(ii-1,11)) - Sol_Total(ii-1,12) );
dis2(6) = abs( y19(Sol_Total(ii-1,13)) - Sol_Total(ii-1,14) );
max_dis2 = max(dis2);
average2 =  ( dis2(1) + dis2(2) + dis2(3) + dis2(4) + dis2(5) + dis2(6) ) / 6;

y24 = y18(x_curve2) + (max_dis1 + average1) / 2;
y25 = y18(x_curve2) - (max_dis1 + average1) / 2;

y26 = y19(x_curve2) + (max_dis2 + average2) / 2;
y27 = y19(x_curve2) - (max_dis2 + average2) / 2;

LifeTime1 = -1; LifeTime2 = -1;
if (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) <= MM
    for i = Spark(ii+5) + 1 : (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) - 1
        if ( Typ_Pri(i+1) >= y24(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y24(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y24(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y24(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y25(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y25(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y25(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y25(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        end
    end
else
    for i = Spark(ii+5) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y24(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y24(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y24(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y24(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y25(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y25(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y25(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y25(i-Spark(ii+5)) )
            LifeTime1 = i+1;
            break
        end
    end
end
if (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) <= MM
    for i = Spark(ii+5) + 1 : (tt_T(ii+5) + 2 * t7_t1(ii-1,1)) - 1
        if ( Typ_Pri(i+1) >= y26(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y26(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y26(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y26(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y27(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y27(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y27(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y27(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        end
    end
else
    for i = Spark(ii+5) + 1 : MM - 1
        if ( Typ_Pri(i+1) >= y26(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y26(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y26(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y26(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) <= y27(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) > y27(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        elseif ( Typ_Pri(i+1) >= y27(i-Spark(ii+5)+1) ) && ( Typ_Pri(i) < y27(i-Spark(ii+5)) )
            LifeTime2 = i+1;
            break
        end
    end
end
if (LifeTime1 == -1) && (LifeTime2 == -1)
    LifeTimeB(ii-1,1) = -1;
elseif (LifeTime1 ~= -1) && (LifeTime2 == -1)
    LifeTimeB(ii-1,1) = LifeTime1; % Life time for with nise
elseif (LifeTime1 == -1) && (LifeTime2 ~= -1)
    LifeTimeB(ii-1,1) = LifeTime2; % Life time for with nise
elseif (LifeTime1 ~= -1) && (LifeTime2 ~= -1)
    LifeTimeB(ii-1,1) = max(LifeTime1,LifeTime2); % Life time for with nise
else
    LifeTimeB(ii-1,1) = -1;
end

% hold on
% plot(x_curve2,y24,'b-.',x_curve2,y25,'b-.');
% hold on
% plot(x_curve2,y26,'g-.',x_curve2,y27,'g-.');
% hold on
% if LifeTime1 ~= -1
%     p12 = line(LifeTime1,Typ_Pri(LifeTime1),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p12.Marker = 's';
% end
% hold on
% if LifeTime2 ~= -1
%     p13 = line(LifeTime2,Typ_Pri(LifeTime2),'Color',[.2 .8 1.0],'linewidth',1.5);
%     p13.Marker = 's';
% end

%--------------------------------  Additional  ----------------------------
% x_curve3 = Sol_Total(ii-1,1) : 2*( tt_T(ii+5) + t7_t1(ii-1,1) );
% x_curve4 = Sol_Total(ii-1,3) : 2*( tt_T(ii+5) + t7_t1(ii-1,1) );
% 
% y12 = a1_Previous * x_curve3 + a0_Previous ...
%     + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * x_curve3 + Sol_Total(ii-1,17));
% 
% y13 = a1 * x_curve4 + a0 ...
%     + Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * x_curve4 + Sol_Total(ii-1,26));
% 
% if 2*( tt_T(ii+5) + t7_t1(ii-1,1) ) > MM
%     figure()
%     plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
%     hold on
%     plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
%     hold on
% else
%     figure()
%     plot(Sol_Total(ii-1,1):2*( tt_T(ii+5) + t7_t1(ii-1,1) ), ...
%         Typ_Pri(Sol_Total(ii-1,1):2*( tt_T(ii+5) + t7_t1(ii-1,1) )),'k.');
%     hold on
%     bi = ii + 5;
%     while ( bi < length(tt_T) ) && ( tt_T(bi) < 2*( tt_T(ii+5) + t7_t1(ii-1,1) ))
%         bi = bi + 1;
%     end
%     plot(tt_T(ii-1:bi),b_T(ii-1:bi),'mo','linewidth',2);
%     hold on
% end
% plot(x_curve3,y12,'b',x_curve4,y13,'g');
% 
% % With noise
% y14 = a1_Previous * x_curve3 + a0_Previous ...
%     + Sol_Total(ii-1,15) * sin(Sol_Total(ii-1,16) * x_curve3 + Sol_Total(ii-1,17)) ...
%     + Sol_Total(ii-1,18) * sin(Sol_Total(ii-1,19) * x_curve3 + Sol_Total(ii-1,20));
% 
% y15 = a1 * x_curve4 + a0 ...
%     + Sol_Total(ii-1,24) * sin(Sol_Total(ii-1,25) * x_curve4 + Sol_Total(ii-1,26)) ...
%     + Sol_Total(ii-1,27) * sin(Sol_Total(ii-1,28) * x_curve4 + Sol_Total(ii-1,29));
% 
% if 2*( tt_T(ii+5) + t7_t1(ii-1,1) ) > MM
%     figure()
%     plot(A(Sol_Total(ii-1,1):MM,8), Typ_Pri(Sol_Total(ii-1,1):MM),'k.');
%     hold on
%     plot(tt_T(ii-1:length(tt_T)),b_T(ii-1:length(b_T)),'mo','linewidth',2);
%     hold on
% else
%     figure()
%     plot(Sol_Total(ii-1,1):2*( tt_T(ii+5) + t7_t1(ii-1,1) ), ...
%         Typ_Pri(Sol_Total(ii-1,1):2*( tt_T(ii+5) + t7_t1(ii-1,1) )),'k.');
%     hold on
%     bi = ii + 5;
%     while ( bi < length(tt_T) ) && ( tt_T(bi) < 2*( tt_T(ii+5) + t7_t1(ii-1,1) ))
%         bi = bi + 1;
%     end
%     plot(tt_T(ii-1:bi),b_T(ii-1:bi),'mo','linewidth',2);
%     hold on
% end
% plot(x_curve3,y14,'b',x_curve4,y15,'g');

end

%-----------------------------  Write table to file  ----------------------

CT_1 = Data_DT2(:,1); % Continuous time
DT_1 = Sol_Total(:,1); % Discrete time
ST_1 = Time_Spark(:,1); % Spark Time (Delay)
De_1 = Time_Delay(:,1); % Delay
Extreme_distance_with_spark_1 = ST_1 - DT_1;
Phase_change_time_1 = ST_1 - De_1;
Pr_1 = Sol_Total(:,2); % Typical Price
CT_2 = Data_DT2(:,2);
DT_2 = Sol_Total(:,3);
ST_2 = Time_Spark(:,2);
De_2 = Time_Delay(:,2);
Extreme_distance_with_spark_2 = ST_2 - DT_2;
Phase_change_time_2 = ST_2 - De_2;
Pr_2 = Sol_Total(:,4);
CT_3 = Data_DT2(:,3);
DT_3 = Sol_Total(:,5);
ST_3 = Time_Spark(:,3);
De_3 = Time_Delay(:,3);
Extreme_distance_with_spark_3 = ST_3 - DT_3;
Phase_change_time_3 = ST_3 - De_3;
Pr_3 = Sol_Total(:,6);
CT_4 = Data_DT2(:,4);
DT_4 = Sol_Total(:,7);
ST_4 = Time_Spark(:,4);
De_4 = Time_Delay(:,4);
Extreme_distance_with_spark_4 = ST_4 - DT_4;
Phase_change_time_4 = ST_4 - De_4;
Pr_4 = Sol_Total(:,8);
CT_5 = Data_DT2(:,5);
DT_5 = Sol_Total(:,9);
ST_5 = Time_Spark(:,5);
De_5 = Time_Delay(:,5);
Extreme_distance_with_spark_5 = ST_5 - DT_5;
Phase_change_time_5 = ST_5 - De_5;
Pr_5 = Sol_Total(:,10);
CT_6 = Data_DT2(:,6);
DT_6 = Sol_Total(:,11);
ST_6 = Time_Spark(:,6);
De_6 = Time_Delay(:,6);
Extreme_distance_with_spark_6 = ST_6 - DT_6;
Phase_change_time_6 = ST_6 - De_6;
Pr_6 = Sol_Total(:,12);
CT_7 = Data_DT2(:,7);
DT_7 = Sol_Total(:,13);
ST_7 = Time_Spark(:,7);
De_7 = Time_Delay(:,7);
Extreme_distance_with_spark_7 = ST_7 - DT_7;
Phase_change_time_7 = ST_7 - De_7;
Pr_7 = Sol_Total(:,14);

a_11 = Sol_Total(:,15);
w_11 = Sol_Total(:,16);
Theta_11 = Sol_Total(:,17);
a_12 = Sol_Total(:,18);
w_12 = Sol_Total(:,19);
Theta_12 = Sol_Total(:,20);
lambda_1 = Sol_Total(:,21);
Coeff_1 = Sol_Total(:,22);
HASmid_1 = Sol_Total(:,23);

a_21 = Sol_Total(:,24);
w_21 = Sol_Total(:,25);
Theta_21 = Sol_Total(:,26);
a_22 = Sol_Total(:,27);
w_22 = Sol_Total(:,28);
Theta_22 = Sol_Total(:,29);
lambda_2 = Sol_Total(:,30);
Coeff_2 = Sol_Total(:,31);
HASmid_2 = Sol_Total(:,32);

%-----------------------------------  DSP  --------------------------------
w_11_a = abs(w_11);
w_12_a = abs(w_12);
w_21_a = abs(w_21);
w_22_a = abs(w_22);

for i = 1 : length(w_11_a)
    if w_11_a(i) >= w_12_a(i)
        W_Max_Min_1(i,1) = 200 * ( w_11_a(i) - w_12_a(i) ) / ( w_11_a(i) + w_12_a(i) );
        if W_Max_Min_1(i,1) >= 50
            W_Max_Min_50_1(i,1) = 1;
        else
            W_Max_Min_50_1(i,1) = 0;
        end
    elseif w_12_a(i) > w_11_a(i)
        W_Max_Min_1(i,1) = 200 * ( w_12_a(i) - w_11_a(i) ) / ( w_11_a(i) + w_12_a(i) );
        if W_Max_Min_1(i,1) >= 50
            W_Max_Min_50_1(i,1) = 1;
        else
            W_Max_Min_50_1(i,1) = 0;
        end
    end
    if w_21_a(i) >= w_22_a(i)
        W_Max_Min_2(i,1) = 200 * ( w_21_a(i) - w_22_a(i) ) / ( w_21_a(i) + w_22_a(i) );
        if W_Max_Min_2(i,1) >= 50
            W_Max_Min_50_2(i,1) = 1;
        else
            W_Max_Min_50_2(i,1) = 0;
        end
    elseif w_22_a(i) > w_21_a(i)
        W_Max_Min_2(i,1) = 200 * ( w_22_a(i) - w_21_a(i) ) / ( w_21_a(i) + w_22_a(i) );
        if W_Max_Min_2(i,1) >= 50
            W_Max_Min_50_2(i,1) = 1;
        else
            W_Max_Min_50_2(i,1) = 0;
        end
    end
end
%--------------------------  Stability coefficient  -----------------------
i1 = 1;
i5 = 0; % DSP total counter
i6 = 0; % DSP units counter
SC = 0;
if (Trend_Change(i1,1) == 1) && (Trend_Change(i1+1,1) == 1) % 1-1
    %Stability_coef_Ex(i1,1) = 0;
    Stability_coef_Ex(i1+1,1) = 0; % Stability coefficient
    Stability_coef_Ex(i1+1,4) = 1001; % Additional
    DSP(i1+1,1) = 0; % DSP
    i1 = i1 + 1;
elseif (Trend_Change(i1,1) == 1) && (Trend_Change(i1+1,1) == 0) % 1-0
    SC = SC + 7;
    i1 = i1 + 1;
    i5 = i5 + 1; % DSP total counter
    if W_Max_Min_50_2(i1,1) == 1
        i6 = i6 + 1; % DSP units counter
    end
    if i1 < length(Trend_Change)%**********
        while (Trend_Change(i1+1,1) == 0) && (i1 + 1 < length(Trend_Change))
            SC = SC + 1;
            i1 = i1 + 1;
            i5 = i5 + 1; % DSP total counter
            if W_Max_Min_50_2(i1,1) == 1
                i6 = i6 + 1; % DSP units counter
            end
        end
        i1 = i1 + 1;
    end%**********
    Stability_coef_Ex(i1,1) = SC; % Stability coefficient
    Stability_coef_Ex(i1,2) = 2; % Start candle time
    Stability_coef_Ex(i1,3) = 7 + (i1 - 2); % End candle time
    Stability_coef_Ex(i1,4) = 1002; % Additional
    Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
        - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
    Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
        / Stability_coef_Ex(i1,1);
    
    Stability_coef_Ex(i1,5) = i5; % DSP total counter
    Stability_coef_Ex(i1,6) = i6; % DSP units counter
    DSP(i1,1) = 100 * i6 / i5; % DSP
elseif (Trend_Change(i1,1) == 0) && (Trend_Change(i1+1,1) == 1) % 0-1
    i5 = i5 + 1; % DSP total counter
    if W_Max_Min_50_2(i1,1) == 1
        i6 = i6 + 1; % DSP units counter
    end
    i1 = i1 + 1;
    Stability_coef_Ex(i1,1) = 7; % Stability coefficient
    Stability_coef_Ex(i1,2) = 1; % Start candle time
    Stability_coef_Ex(i1,3) = 7 + (i1 - 2); % End candle time
    Stability_coef_Ex(i1,4) = 1003; % Additional
    Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
        - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
    Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
        / Stability_coef_Ex(i1,1);
    
    Stability_coef_Ex(i1,5) = i5; % DSP total counter
    Stability_coef_Ex(i1,6) = i6; % DSP units counter
    DSP(i1,1) = 100 * i6 / i5; % DSP
elseif (Trend_Change(i1,1) == 0) && (Trend_Change(i1+1,1) == 0) % 0-0
    i5 = i5 + 1; % DSP total counter
    if W_Max_Min_50_2(i1,1) == 1
        i6 = i6 + 1; % DSP units counter
    end
    SC = SC + 8;
    i1 = i1 + 1;
    i5 = i5 + 1; % DSP total counter
    if W_Max_Min_50_2(i1,1) == 1
        i6 = i6 + 1; % DSP units counter
    end
    if i1 < length(Trend_Change)%**********
        while (Trend_Change(i1+1,1) == 0) && (i1 + 1 < length(Trend_Change))
            SC = SC + 1;
            i1 = i1 + 1;
            i5 = i5 + 1; % DSP total counter
            if W_Max_Min_50_2(i1,1) == 1
                i6 = i6 + 1; % DSP units counter
            end
        end
        i1 = i1 + 1;
    end%**********
    Stability_coef_Ex(i1,1) = SC; % Stability coefficient
    Stability_coef_Ex(i1,2) = 1; % Start candle time
    Stability_coef_Ex(i1,3) = 7 + (i1 - 2); % End candle time
    Stability_coef_Ex(i1,4) = 1004; % Additional
    Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
        - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
    Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
        / Stability_coef_Ex(i1,1);
    
    Stability_coef_Ex(i1,5) = i5; % DSP total counter
    Stability_coef_Ex(i1,6) = i6; % DSP units counter
    DSP(i1,1) = 100 * i6 / i5; % DSP
end

while i1 < length(Trend_Change) - 1
    i5 = 0; % DSP total counter
    i6 = 0; % DSP units counter
    SC = 0;
    if (Trend_Change(i1,1) == 1) && (Trend_Change(i1+1,1) == 1) %1-1
        Stability_coef_Ex(i1+1,1) = 0;
        Stability_coef_Ex(i1+1,4) = 2001;
        i1 = i1 + 1;
    elseif (Trend_Change(i1,1) == 1) && (Trend_Change(i1+1,1) == 0) % 1-0
        TSC = i1;
        SC = SC + 7;
        i1 = i1 + 1;
        i5 = i5 + 1; % DSP total counter
        if W_Max_Min_50_2(i1,1) == 1
            i6 = i6 + 1; % DSP units counter
        end
        while (Trend_Change(i1+1,1) == 0) && (i1 + 1 < length(Trend_Change))
            SC = SC + 1;
            i1 = i1 + 1;
            i5 = i5 + 1; % DSP total counter
            if W_Max_Min_50_2(i1,1) == 1
                i6 = i6 + 1; % DSP units counter
            end
        end
        i1 = i1 + 1;
        Stability_coef_Ex(i1,1) = SC;
        Stability_coef_Ex(i1,2) = TSC + 1; % or Stability_coef_Ex(i1,3)=i1
        Stability_coef_Ex(i1,3) = 7 + (i1 - 2); % End candle time
        Stability_coef_Ex(i1,4) = 2002;
        Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
            - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
        Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
            / Stability_coef_Ex(i1,1);
        
        Stability_coef_Ex(i1,5) = i5; % DSP total counter
        Stability_coef_Ex(i1,6) = i6; % DSP units counter
        DSP(i1,1) = 100 * i6 / i5; % DSP
    end
end

if (i1 == 2 && Trend_Change(1,1) == 1 && Trend_Change(2,1) == 1)%************
    Stability_coef_Ex(i1,4) = 3004; % Additional
    Stability_coefficient_time(i1,1) = 0;
    Stability_coefficient_ratio(i1,1) = 0;
elseif (i1 == 2 && Trend_Change(1,1) == 1 && Trend_Change(2,1) == 0)%************
    Stability_coef_Ex(i1,3) = 7 + (i1 - 1); % End candle time
    Stability_coef_Ex(i1,4) = 3005; % Additional
    Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
        - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
    Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
        / Stability_coef_Ex(i1,1);
elseif (i1 == 2 && Trend_Change(1,1) == 0 && Trend_Change(2,1) == 1)%************
    Stability_coef_Ex(i1,3) = 7; % End candle time
    Stability_coef_Ex(i1,4) = 3006; % Additional
    Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
        - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
    Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
        / Stability_coef_Ex(i1,1);
elseif (i1 == 2 && Trend_Change(1,1) == 0 && Trend_Change(2,1) == 0)%************
    Stability_coef_Ex(i1,3) = 7 + (i1 - 1); % End candle time
    Stability_coef_Ex(i1,4) = 3007; % Additional
    Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
        - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
    Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
        / Stability_coef_Ex(i1,1);
elseif i1 == length(Trend_Change) - 1
    if (Trend_Change(i1,1) == 1) && (Trend_Change(i1+1,1) == 0) % 1-0
        TSC = i1;
        i1 = i1 + 1;
        Stability_coef_Ex(i1,1) = 7;
        Stability_coef_Ex(i1,2) = TSC + 1; % or Stability_coef_Ex(i1,3)=i1
        Stability_coef_Ex(i1,3) = 7 + (i1 - 1); % End candle time
        Stability_coef_Ex(i1,4) = 3001;
        Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
            - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
        Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
            / Stability_coef_Ex(i1,1);
        
        i5 = 1; % DSP total counter
        i6 = 0; % DSP units counter
        if W_Max_Min_50_2(i1,1) == 1
            i6 = 1; % DSP units counter
        end
        Stability_coef_Ex(i1,5) = i5; % DSP total counter
        Stability_coef_Ex(i1,6) = i6; % DSP units counter
        DSP(i1,1) = 100 * i6 / i5; % DSP
    elseif (Trend_Change(i1,1) == 1) && (Trend_Change(i1+1,1) == 1) % 1-1
        i1 = i1 + 1;
        Stability_coef_Ex(i1,1) = 0;
        Stability_coef_Ex(i1,4) = 3002;
        Stability_coefficient_time(i1,1) = 0;
        Stability_coefficient_ratio(i1,1) = 0;
        
        DSP(i1,1) = 0; % DSP
    end
elseif i1 == length(Trend_Change)
    if Trend_Change(i1,1) == 0
        Stability_coef_Ex(i1,1) = Stability_coef_Ex(i1,1) + 1;
        Stability_coef_Ex(i1,3) = 7 + (i1 - 1); % End candle time
        Stability_coef_Ex(i1,4) = 3003;
        Stability_coefficient_time(i1,1) = EXPRICE(Stability_coef_Ex(i1,3),3) ....
            - EXPRICE(Stability_coef_Ex(i1,2),3) + 1;
        Stability_coefficient_ratio(i1,1) = Stability_coefficient_time(i1,1) ...
            / Stability_coef_Ex(i1,1);
        
        i5 = i5 + 1; % DSP total counter
        if W_Max_Min_50_2(i1,1) == 1
            i6 = i6 + 1; % DSP units counter
        end
        Stability_coef_Ex(i1,5) = i5; % DSP total counter
        Stability_coef_Ex(i1,6) = i6; % DSP units counter
        DSP(i1,1) = 100 * i6 / i5; % DSP
    end
end

Stability_coefficient = Stability_coef_Ex(:,1);
% DSP1(:,1) = Stability_coef_Ex(:,5); % DSP total counter
% DSP2(:,1) = Stability_coef_Ex(:,6); % DSP units counter

%---------------------------  Three reward to risk  -----------------------
rr1 = 0;
%IgnoreSequence_Sel_R1 = NumberZero(SequenceZeroOne_Sel_R1);
[IgnoreSequence_Sel_R1,MaxZeroCunt_Sel_R1] = NumberZeroCont(SequenceZeroOne_Sel_R1);
if IgnoreSequence_Sel_R1 == 1
    rr1 = rr1 + 1;
    %RR(rr1,1) = RR_TP_SL(TypePosition_Sel_R1, 1);
    RR(rr1,1) = MaxZeroCunt_Sel_R1;
    RR2(rr1,1) = {'Direct_100_100'};
    Additional(rr1,1) = IgnoreSequence_Sel_R1; % Additional
end

%IgnoreSequence_Rev_R1 = NumberZero(SequenceZeroOne_Rev_R1);
[IgnoreSequence_Rev_R1,MaxZeroCunt_Rev_R1] = NumberZeroCont(SequenceZeroOne_Rev_R1);
if IgnoreSequence_Rev_R1 == 1
    rr1 = rr1 + 1;
    %RR(rr1,1) = RR_TP_SL(TypePosition_Rev_R1, 1);
    RR(rr1,1) = MaxZeroCunt_Rev_R1;
    RR2(rr1,1) = {'Reverse_100_100'};
    Additional(rr1,1) = IgnoreSequence_Rev_R1; % Additional
end

%--------------------------------------------------------------------------
%IgnoreSequence_Sel_R2 = NumberZero(SequenceZeroOne_Sel_R2);
[IgnoreSequence_Sel_R2,MaxZeroCunt_Sel_R2] = NumberZeroCont(SequenceZeroOne_Sel_R2);
if IgnoreSequence_Sel_R2 == 1
    rr1 = rr1 + 1;
    %RR(rr1,1) = RR_TP_SL(TypePosition_Sel_R2, 1);
    RR(rr1,1) = MaxZeroCunt_Sel_R2;
    RR2(rr1,1) = {'Direct_50_50'};
    Additional(rr1,1) = IgnoreSequence_Sel_R2; % Additional
end

%IgnoreSequence_Rev_R2 = NumberZero(SequenceZeroOne_Rev_R2);
[IgnoreSequence_Rev_R2,MaxZeroCunt_Rev_R2] = NumberZeroCont(SequenceZeroOne_Rev_R2);
if IgnoreSequence_Rev_R2 == 1
    rr1 = rr1 + 1;
    %RR(rr1,1) = RR_TP_SL(TypePosition_Rev_R2, 1);
    RR(rr1,1) = MaxZeroCunt_Rev_R2;
    RR2(rr1,1) = {'Reverse_50_50'};
    Additional(rr1,1) = IgnoreSequence_Rev_R2; % Additional
end

%--------------------------------------------------------------------------
%IgnoreSequence_Sel_R3 = NumberZero(SequenceZeroOne_Sel_R3);
[IgnoreSequence_Sel_R3,MaxZeroCunt_Sel_R3] = NumberZeroCont(SequenceZeroOne_Sel_R3);
if IgnoreSequence_Sel_R3 == 1
    rr1 = rr1 + 1;
    %RR(rr1,1) = RR_TP_SL(TypePosition_Sel_R3, 3);
    RR(rr1,1) = MaxZeroCunt_Sel_R3;
    RR2(rr1,1) = {'Direct_150_50'};
    Additional(rr1,1) = IgnoreSequence_Sel_R3; % Additional
end

%IgnoreSequence_Rev_R3 = NumberZero(SequenceZeroOne_Rev_R3);
[IgnoreSequence_Rev_R3,MaxZeroCunt_Rev_R3] = NumberZeroCont(SequenceZeroOne_Rev_R3);
if IgnoreSequence_Rev_R3 == 1
    rr1 = rr1 + 1;
    %RR(rr1,1) = RR_TP_SL(TypePosition_Rev_R3, 3);
    RR(rr1,1) = MaxZeroCunt_Rev_R3;
    RR2(rr1,1) = {'Reverse_150_50'};
    Additional(rr1,1) = IgnoreSequence_Rev_R3; % Additional
end
%--------------------------------------------------------------------------
if rr1 > 0
%     RR_max = max(RR);
%     RR_Index = find(RR == RR_max);
    RR_min = min(RR);
    RR_Index = find(RR == RR_min);
    %disp('max(RR * TP - SL) =')
    disp('Minimum Sequence zeros:')
    disp(RR_min)
    disp('The Best reward to risk:')
    disp( RR2(RR_Index,1) )
elseif rr1 ==0
    disp('We do not have optimal reward to risk.')
end
%--------------------------------------------------------------------------

T = table(CT_1, DT_1, ST_1, De_1, Extreme_distance_with_spark_1, Phase_change_time_1, Pr_1, ...
    CT_2, DT_2, ST_2, De_2, Extreme_distance_with_spark_2, Phase_change_time_2, Pr_2, ...
    CT_3, DT_3, ST_3, De_3, Extreme_distance_with_spark_3, Phase_change_time_3, Pr_3, ...
    CT_4, DT_4, ST_4, De_4, Extreme_distance_with_spark_4, Phase_change_time_4, Pr_4, ...
    CT_5, DT_5, ST_5, De_5, Extreme_distance_with_spark_5, Phase_change_time_5, Pr_5, ...
    CT_6, DT_6, ST_6, De_6, Extreme_distance_with_spark_6, Phase_change_time_6, Pr_6, ...
    CT_7, DT_7, ST_7, De_7, Extreme_distance_with_spark_7, Phase_change_time_7, Pr_7, ...
    a_11, w_11, Theta_11, a_12, w_12, Theta_12, lambda_1, Coeff_1, HASmid_1, ...
    a_21, w_21, Theta_21, a_22, w_22, Theta_22, lambda_2, Coeff_2, HASmid_2, ...
    Trend_Num_1, Point_Num_1, Trend_1, Trend_Mac_1, Up_Down_1, Up_Down_Mac_1, ...
    Min_Max_1, W_Max_Min_1, W_Max_Min_50_1, ...
    Trend_Num_2, Point_Num_2, Trend_2, Trend_Mac_2, Up_Down_2, Up_Down_Mac_2, ...
    Min_Max_2, W_Max_Min_2, W_Max_Min_50_2, t7_t1, Trend_Change, ...
    DSP, Stability_coefficient, Stability_coefficient_time, Stability_coefficient_ratio, ...
    First_encounter_CT, First_encounter_DT, First_PEncounter, ...
    Second_encounter_CT, Second_encounter_DT, Second_PEncounter, ...
    Third_encounter_CT, Third_encounter_DT, Third_PEncounter, ...
    First_Forecast_CT, First_Forecast_DT, First_PForecast, ...
    Second_Forecast_CT, Second_Forecast_DT, Second_PForecast, ...
    Third_Forecast_CT, Third_Forecast_DT, Third_PForecast, ...
    XNUD_Encounter, Hunt_Encounter, Deviation_E, Max_PercentE, ...
    XNUD_Forecast, Hunt_Forecast, ...
    Gravity_tE, Gravity_PE, XNUD_EG, Hunt_EG, Deviation_EG, Max_PercentEG, ...
    Gravity_tF, Gravity_PF, XNUD_FG, Hunt_FG, Deviation_FG, Max_PercentFG, ...
    LifeTimeA, LifeTimeB, ...
    TypePosition_E, TypePos_MaxPer_Dev_E, TypePos_E, ...
    TypePosition_F, SellBuy, ...
    TimeOpen, PriceOpen, TimeClose, PriceClose, ...
    New_lineP, New_lineP_D, Hunt_Price, ...
    XNUD_Forecast_Sel_R1, Hunt_Forecast_Sel_R1, SellBuy_Sel_R1, TimeOpen_Sel_R1, ...
    PriceOpen_Sel_R1, TimeClose_Sel_R1, PriceClose_Sel_R1, TypePosition_Sel_R1, ...
    Cut_Number, ...
    XNUD_Forecast_Rev_R1, Hunt_Forecast_Rev_R1, SellBuy_Rev_R1, TimeOpen_Rev_R1, ...
    PriceOpen_Rev_R1, TimeClose_Rev_R1, PriceClose_Rev_R1, TypePosition_Rev_R1, ...
    XNUD_Forecast_Sel_R2, Hunt_Forecast_Sel_R2, SellBuy_Sel_R2, TimeOpen_Sel_R2, ...
    PriceOpen_Sel_R2, TimeClose_Sel_R2, PriceClose_Sel_R2, TypePosition_Sel_R2, ...
    XNUD_Forecast_Rev_R2, Hunt_Forecast_Rev_R2, SellBuy_Rev_R2, TimeOpen_Rev_R2, ...
    PriceOpen_Rev_R2, TimeClose_Rev_R2, PriceClose_Rev_R2, TypePosition_Rev_R2, ...
    XNUD_Forecast_Sel_R3, Hunt_Forecast_Sel_R3, SellBuy_Sel_R3, TimeOpen_Sel_R3, ...
    PriceOpen_Sel_R3, TimeClose_Sel_R3, PriceClose_Sel_R3, TypePosition_Sel_R3, ...
    XNUD_Forecast_Rev_R3, Hunt_Forecast_Rev_R3, SellBuy_Rev_R3, TimeOpen_Rev_R3, ...
    PriceOpen_Rev_R3, TimeClose_Rev_R3, PriceClose_Rev_R3, TypePosition_Rev_R3); % Deviation_F, Max_PercentF
% TypePos_MaxPer_Dev_F, TypePos_F, FinalPrice_F, PiP, NetProfit, PerPosition,
writetable(T,'Final_Result_3.csv')
writetable(T,'Final_Result_3.xlsx')
%type 'Final_Result.csv'

Sequence = table(SequenceZeroOne');
writetable(Sequence,'SequenceZeroOne_3.csv')

Sequence_Sel_R1 = table(SequenceZeroOne_Sel_R1');
writetable(Sequence_Sel_R1,'SequenceZeroOne_Sel_R1_3.csv')
Sequence_Rev_R1 = table(SequenceZeroOne_Rev_R1');
writetable(Sequence_Rev_R1,'SequenceZeroOne_Rev_R1_3.csv')

Sequence_Sel_R2 = table(SequenceZeroOne_Sel_R2');
writetable(Sequence_Sel_R2,'SequenceZeroOne_Sel_R2_3.csv')
Sequence_Rev_R2 = table(SequenceZeroOne_Rev_R2');
writetable(Sequence_Rev_R2,'SequenceZeroOne_Rev_R2_3.csv')

Sequence_Sel_R3 = table(SequenceZeroOne_Sel_R3');
writetable(Sequence_Sel_R3,'SequenceZeroOne_Sel_R3_3.csv')
Sequence_Rev_R3 = table(SequenceZeroOne_Rev_R3');
writetable(Sequence_Rev_R3,'SequenceZeroOne_Rev_R3_3.csv')

end

toc
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% figure()
% plot(x,y,'r.',x,y5,'b')
% hold on
% plot(ttt,bb,'c*',tt,b,'kh','linewidth',2)
% title('System n: r+a_1*sin(w_1t+\theta_1)+a_2*sin(w_2t+\theta_2)')
% legend('Real Data','r+a_1*sin(w_1t+\theta_1)+a_2*sin(w_2t+\theta_2)',...
%     'Mide Points','Six Points')
%--------------------------------------------------------------------------


