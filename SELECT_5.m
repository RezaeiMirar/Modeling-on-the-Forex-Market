function [solution,Maximum] = SELECT_5(m_1, M, lambda, ttt, bb, x, r, S1_1, S2) % SELECT_3
fix = 12 * lambda * pi / ( ttt(m_1-1) - ttt(1) );
W1 = linspace(0,fix,M); %w_1

k1 = 0;
for i = 1 : m_1-1
    S1_2 = S1_1(bb(i), ttt(i));
    S2_1 =@(w) S2(ttt(i), w);
    for j = 1 : M-1
        S2_2 = S2_1(W1(j+1));
        k1 = k1 + 1;
        Index_1(k1,1) = S1_2 / S2_2;
        Index_1(k1,2) = W1(j+1);
        SUM(k1) = 0;
        for k = 1 : m_1-1
            S31(k) = abs( S1_1(bb(k), ttt(k)) - Index_1(k1,1) * S2(ttt(k), W1(j+1)) );
            SUM(k1) = SUM(k1) + S31(k);
        end
    end
end
Minimum = min(SUM);
[Index_4,Index_5] = find(SUM==Minimum);
solution(1,1:2) = [Index_1(Index_5,1), Index_1(Index_5,2)]; % Solution Set

y2 = r + solution(1) * sin(solution(2)*x);

Index_2(1,1) = 0;
for j = 1 : m_1-2 % 1:4 => Nomber of mide point
    i1 = 0;
    for i = ( floor(ttt(j)) - ( floor(ttt(1)) - 1 ) ) : ( floor(ttt(j+1)) - ( floor(ttt(1)) - 1 ) - 2 )
        if ( y2(i+1)<y2(i) && y2(i+1)<y2(i+2) ) || ( y2(i+1)>y2(i) && y2(i+1)>y2(i+2) )
            i1 = i1 + 1;
            Index_2(j) = i1;
%             Index_2(j,1) = i1;
%             Index_2(j,i1+1) = i + 1;
%             Index_2(j,i1+1+5) = floor(ttt(j)) - ( floor(ttt(1)) - 1 );
%             Index_2(j,i1+1+10) = floor(ttt(j+1)) - ( floor(ttt(1)) - 1 ) - 2;
        else
            continue
        end
    end
end
Maximum = max(Index_2);
%Maximum = max(Index_2(:,1));

end

