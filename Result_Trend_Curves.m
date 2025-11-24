function Trend_Curves_Cut = Result_Trend_Curves(Ext1_Cut_1_time, Ext2_Cut_1_time, ...
    Ext1_Cut_1_MinMax, Ext2_Cut_1_MinMax, YForecast, S_E_1, S_E_2)

i2 = 0;
if Ext1_Cut_1_time == Ext2_Cut_1_time
    if Ext1_Cut_1_MinMax == 1 % Min
        Trend_Curve1_Cut1 = 2; % Increase
    elseif Ext1_Cut_1_MinMax == 2 % Max
        Trend_Curve1_Cut1 = 1; % Decrease
    end
    if Ext2_Cut_1_MinMax == 1 % Min
        Trend_Curve2_Cut1 = 2; % Increase
    elseif Ext2_Cut_1_MinMax == 2 % Max
        Trend_Curve2_Cut1 = 1; % Decrease
    end
elseif Ext1_Cut_1_time > Ext2_Cut_1_time
    for i = 1 : size(S_E_1,1)
        if (S_E_1(i,1) < YForecast) && (S_E_1(i,1) > Ext2_Cut_1_time)
            i2 = i2 + 1;
            SEIndex(i2,1) = S_E_1(i,1);
            SEIndex(i2,2) = S_E_1(i,3);
        end
    end
    if i2 == 0
        Trend_Curve1_Cut1 = 0;
    end
    if Ext2_Cut_1_MinMax == 1 % Min
        Trend_Curve2_Cut1 = 2; % Increase
    elseif Ext2_Cut_1_MinMax == 2 % Max
        Trend_Curve2_Cut1 = 1; % Decrease
    end
elseif Ext1_Cut_1_time < Ext2_Cut_1_time
    for i = 1 : size(S_E_2,1)
        if (S_E_2(i,1) < YForecast) && (S_E_2(i,1) > Ext1_Cut_1_time)
            i2 = i2 + 1;
            SEIndex(i2,1) = S_E_2(i,1);
            SEIndex(i2,2) = S_E_2(i,3);
        end
    end
    if i2 == 0
        Trend_Curve2_Cut1 = 0;
    end
    if Ext1_Cut_1_MinMax == 1 % Min
        Trend_Curve1_Cut1 = 2; % Increase
    elseif Ext1_Cut_1_MinMax == 2 % Max
        Trend_Curve1_Cut1 = 1; % Decrease
    end
end
%-------------------------------  One Extreme  ----------------------------
if (i2 == 1) && (Ext1_Cut_1_time >= Ext2_Cut_1_time) && (SEIndex(1,2) == 2)
    s1 = Ext1_Cut_1_time - Ext2_Cut_1_time + 1; % Increase
    s2 = YForecast - Ext1_Cut_1_time; % Decrease
    if s1 >= 2 * s2
        Trend_Curve1_Cut1 = 2; % Increase
    elseif s2 >= 2 * s1
        Trend_Curve1_Cut1 = 1; % Decrease
    else
        Trend_Curve1_Cut1 = 0; % Range
    end
elseif (i2 == 1) && (Ext1_Cut_1_time >= Ext2_Cut_1_time) && (SEIndex(1,2) == 1)
    s1 = Ext1_Cut_1_time - Ext2_Cut_1_time + 1; % Decrease
    s2 = YForecast - Ext1_Cut_1_time; % Increase
    if s1 >= 2 * s2
        Trend_Curve1_Cut1 = 1; % Decrease
    elseif s2 >= 2 * s1
        Trend_Curve1_Cut1 = 2; % Increase
    else
        Trend_Curve1_Cut1 = 0; % Range
    end
elseif (i2 == 1) && (Ext2_Cut_1_time >= Ext1_Cut_1_time) && (SEIndex(1,2) == 2)
    s1 = Ext2_Cut_1_time - Ext1_Cut_1_time + 1; % Increase
    s2 = YForecast - Ext2_Cut_1_time; % Decrease
    if s1 >= 2 * s2
        Trend_Curve2_Cut1 = 2; % Increase
    elseif s2 >= 2 * s1
        Trend_Curve2_Cut1 = 1; % Decrease
    else
        Trend_Curve2_Cut1 = 0; % Range
    end
elseif (i2 == 1) && (Ext2_Cut_1_time >= Ext1_Cut_1_time) && (SEIndex(1,2) == 1)
    s1 = Ext2_Cut_1_time - Ext1_Cut_1_time + 1; % Decrease
    s2 = YForecast - Ext2_Cut_1_time; % Increase
    if s1 >= 2 * s2
        Trend_Curve2_Cut1 = 1; % Decrease
    elseif s2 >= 2 * s1
        Trend_Curve2_Cut1 = 2; % Increase
    else
        Trend_Curve2_Cut1 = 0; % Range
    end
%-------------------------------  Tow Extreme  ----------------------------
elseif (i2 == 2) && (Ext1_Cut_1_time >= Ext2_Cut_1_time) && (SEIndex(1,2) == 2)
    s1 = SEIndex(1,1) - Ext2_Cut_1_time + 1; % Increase
    s2 = SEIndex(2,1) - SEIndex(1,1); % Decrease
    s3 = YForecast - SEIndex(2,1); % Increase
    if (s1 + s3) >= 2 * s2
        Trend_Curve1_Cut1 = 2; % Increase
    elseif s2 >= 2 * (s1 + s3)
        Trend_Curve1_Cut1 = 1; % Decrease
    else
        Trend_Curve1_Cut1 = 0; % Range
    end
elseif (i2 == 2) && (Ext1_Cut_1_time >= Ext2_Cut_1_time) && (SEIndex(1,2) == 1)
    s1 = SEIndex(1,1) - Ext2_Cut_1_time + 1; % Decrease
    s2 = SEIndex(2,1) - SEIndex(1,1); % Increase
    s3 = YForecast - SEIndex(2,1); % Decrease
    if (s1 + s3) >= 2 * s2
        Trend_Curve1_Cut1 = 1; % Decrease
    elseif s2 >= 2 * (s1 + s3)
        Trend_Curve1_Cut1 = 2; % Increase
    else
        Trend_Curve1_Cut1 = 0; % Range
    end
elseif (i2 == 2) && (Ext2_Cut_1_time >= Ext1_Cut_1_time) && (SEIndex(1,2) == 2)
    s1 = SEIndex(1,1) - Ext1_Cut_1_time + 1; % Increase
    s2 = SEIndex(2,1) - SEIndex(1,1); % Decrease
    s3 = YForecast - SEIndex(2,1); % Increase
    if (s1 + s3) >= 2 * s2
        Trend_Curve2_Cut1 = 2; % Increase
    elseif s2 >= 2 * (s1 + s3)
        Trend_Curve2_Cut1 = 1; % Decrease
    else
        Trend_Curve2_Cut1 = 0; % Range
    end
elseif (i2 == 2) && (Ext2_Cut_1_time >= Ext1_Cut_1_time) && (SEIndex(1,2) == 1)
    s1 = SEIndex(1,1) - Ext1_Cut_1_time + 1; % Decrease
    s2 = SEIndex(2,1) - SEIndex(1,1); % Increase
    s3 = YForecast - SEIndex(2,1); % Decrease
    if (s1 + s3) >= 2 * s2
        Trend_Curve2_Cut1 = 1; % Decrease
    elseif s2 >= 2 * (s1 + s3)
        Trend_Curve2_Cut1 = 2; % Increase
    else
        Trend_Curve2_Cut1 = 0; % Range
    end
end
%--------------------------------------------------------------------------

if i2 > 2 % Range
    Trend_Curves_Cut = 0;
elseif Trend_Curve1_Cut1 == 0
    Trend_Curves_Cut = 0;
elseif Trend_Curve2_Cut1 == 0
    Trend_Curves_Cut = 0;
elseif Trend_Curve1_Cut1 ~= Trend_Curve2_Cut1
    Trend_Curves_Cut = 0;
elseif Trend_Curve1_Cut1 == Trend_Curve2_Cut1
    Trend_Curves_Cut = 1;
end

end

