function [solution_2,Maximum] = SELECT_6(a4, solution_1, M, S3, S2, b, tt, m_1, r, x, THETA_1) % SELECT_2
W2 = linspace(0, solution_1(2)*a4, M); %w_2

k3 = 0;
for i = 1 : m_1
    S3_1 = S3(b(i), tt(i));
    S4_1 =@(w) S2(tt(i), w);
    for j = 1 : M-1
        S4_2 = S4_1(W2(j+1));
        k3 = k3 + 1;
        Index_8(k3,1) = S3_1 / S4_2;
        Index_8(k3,2) = W2(j+1);
        SUM_2(k3) = 0;
        for k5 = 1 : m_1
            S33(k5) = abs( S3(b(k5), tt(k5)) - Index_8(k3,1) * S2(tt(k5), W2(j+1)) );
            SUM_2(k3) = SUM_2(k3) + S33(k5);
        end
    end
end
Minimum_2 = min(SUM_2);
[Index_9,Index_10] = find(SUM_2==Minimum_2);
solution_2(1,1:2) = [Index_8(Index_10,1), Index_8(Index_10,2)]; % Solution Set

y6 = r + solution_1(1)*sin(solution_1(2)*x+THETA_1) + solution_2(1)*sin(solution_2(2)*x);

Index_14(1,1) = 0;
for j = 1 : m_1-1 % 1:5 => Nomber of six point
    i1 = 0;
    for i = ( tt(j) : tt(j+1) - 2 ) - ( tt(1) - 1 )
        if ( y6(i+1)<y6(i) && y6(i+1)<y6(i+2) ) || ( y6(i+1)>y6(i) && y6(i+1)>y6(i+2) )
            i1 = i1 + 1;
            Index_14(j) = i1;
%             Index_14(j,1) = i1;
%             Index_14(j,i1+1) = i + 1;
%             Index_14(j,i1+1+5) = tt(j) - ( tt(1) - 1 );
%             Index_14(j,i1+1+10) = tt(j+1) - ( tt(1) - 1 ) - 2;
        else
            continue
        end
    end
end
Maximum = max(Index_14);
%Maximum = max(Index_14(:,1));

end

