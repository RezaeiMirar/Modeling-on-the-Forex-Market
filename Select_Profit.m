function [XNUD_Forecast, Hunt_Forecast, SellBuy, TimeOpen, PriceOpen, TimeClose, ...
    PriceClose, TypePosition_F, NSZO] ...
    = Select_Profit(Typ_Pri, YForecast, PForecast, Spark, NSZO, b_T_6, b_T_7, ...
    TimeClose, ii, Digit, MM)

if (b_T_6 < b_T_7) && (b_T_7 < PForecast) % b_T(7) is max
    XNUD_Forecast = {'XU'};
elseif (b_T_6 < b_T_7) && (b_T_7 > PForecast) % b_T(7) is max
    XNUD_Forecast = {'XD'};
elseif (b_T_6 > b_T_7) && (b_T_7 < PForecast) % b_T(7) is min
    XNUD_Forecast = {'NU'};
elseif (b_T_6 > b_T_7) && (b_T_7 > PForecast) % b_T(7) is min
    XNUD_Forecast = {'ND'};
elseif (b_T_6 < b_T_7) && (b_T_7 == PForecast) % b_T(7) is equal to PForecast(1,1).
    XNUD_Forecast = {'XE'};
elseif (b_T_6 > b_T_7) && (b_T_7 == PForecast) % b_T(7) is equal to PForecast(1,1).
    XNUD_Forecast = {'NE'};
end

TimeOpen = Spark;
PriceOpen = Typ_Pri(Spark);

NumberTime = 0;
for w = 1 : ii - 2
    if TimeOpen < TimeClose(ii-1-w,1)
        NumberTime = NumberTime + 1;
    end
end

if NumberTime == 0
    if Typ_Pri(Spark) < PForecast % Hunt Up
        Hunt_Forecast = {'Hunt_U'};
        SellBuy = {'Buy'};
        TypePosition_F = -1; % No mode
        New_lineP = PForecast + 4 * 10^(-(Digit-1));
        New_lineP_D = Typ_Pri(Spark) - abs(PForecast - Typ_Pri(Spark)) + 4 * 10^(-(Digit-1));
        
        if YForecast <= MM + 1 %******
            for i = Spark : YForecast - 1
                if (((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                        || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))) ...
                        && (((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                        || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D)))
                    TypePosition_F = 2; % Deviation
                    NSZO = NSZO + 1;
                    %SequenceZeroOne = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose = New_lineP_D;
                    break
                elseif ((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                        || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))
                    TypePosition_F = 1; % Max_Percent
                    NSZO = NSZO + 1;
                    %SequenceZeroOne = 1;
                    TimeClose(ii-1,1) = i;
                    PriceClose = New_lineP;
                    break
                elseif ((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                        || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D))
                    TypePosition_F = 2; % Deviation
                    NSZO = NSZO + 1;
                    %SequenceZeroOne = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose = New_lineP_D;
                    break
                end
            end
            i = i + 1;
            if TypePosition_F == -1 % No mode
                if Typ_Pri(YForecast) >= Typ_Pri(Spark) + 4 * 10^(-(Digit-1))
                    TypePosition_F = 3; % No mode & in profit
                    TimeClose(ii-1,1) = YForecast;
                    PriceClose = Typ_Pri(YForecast);
                else
                    while (Typ_Pri(i) < Typ_Pri(Spark) + 4 * 10^(-(Digit-1))) ...
                            && (Typ_Pri(i) > New_lineP_D) && (i < MM)
                        i = i + 1;
                    end
                    
                    if i > MM
                        TypePosition_F = 100; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose = -1;
                    elseif Typ_Pri(i) >= Typ_Pri(Spark) + 4 * 10^(-(Digit-1))
                        TypePosition_F = 32; % No mode & Risk Free
                        TimeClose(ii-1,1) = i;
                        PriceClose = Typ_Pri(Spark) + 4 * 10^(-(Digit-1));
                    elseif Typ_Pri(i) <= New_lineP_D
                        TypePosition_F = 33; % No mode & Deviation
                        NSZO = NSZO + 1;
                        %SequenceZeroOne = 0;
                        TimeClose(ii-1,1) = i;
                        PriceClose = New_lineP_D;
                    else
                        TypePosition_F = 50; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose = -1;
                    end
                end
            end
        else
            Hunt_Forecast = {'There is no real data'};
            SellBuy = {'There is no real data'};
            
            TypePosition_F = -1;
            
            TimeClose(ii-1,1) = -1;
            PriceClose = -1;
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark) > PForecast % Hunt Down
        Hunt_Forecast = {'Hunt_D'};
        SellBuy = {'Sell'};
        TypePosition_F = -1; % No mode
        New_lineP = PForecast - 4 * 10^(-(Digit-1));
        New_lineP_D = Typ_Pri(Spark) + abs(PForecast - Typ_Pri(Spark)) - 4 * 10^(-(Digit-1));
        
        if YForecast <= MM + 1 %******
            for i = Spark : YForecast - 1
                if (((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                        || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))) ...
                        && (((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                        || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D)))
                    TypePosition_F = 2; % Deviation
                    NSZO = NSZO + 1;
                    %SequenceZeroOne = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose = New_lineP_D;
                    break
                elseif ((Typ_Pri(i) >= New_lineP) && (Typ_Pri(i+1) < New_lineP)) ...
                        || ((Typ_Pri(i) <= New_lineP) && (Typ_Pri(i+1) > New_lineP))
                    TypePosition_F = 1; % Max_Percent
                    NSZO = NSZO + 1;
                    %SequenceZeroOne = 1;
                    TimeClose(ii-1,1) = i;
                    PriceClose = New_lineP;
                    break
                elseif ((Typ_Pri(i) <= New_lineP_D) && (Typ_Pri(i+1) > New_lineP_D)) ...
                        || ((Typ_Pri(i) >= New_lineP_D) && (Typ_Pri(i+1) < New_lineP_D))
                    TypePosition_F = 2; % Deviation
                    NSZO = NSZO + 1;
                    %SequenceZeroOne = 0;
                    TimeClose(ii-1,1) = i;
                    PriceClose = New_lineP_D;
                    break
                end
            end
            i = i + 1;
            if TypePosition_F == -1 % No mode
                if Typ_Pri(YForecast) <= Typ_Pri(Spark) - 4 * 10^(-(Digit-1))
                    TypePosition_F = 3; % No mode & in profit
                    TimeClose(ii-1,1) = YForecast;
                    PriceClose = Typ_Pri(YForecast);
                else
                    while (Typ_Pri(i) > Typ_Pri(Spark) - 4 * 10^(-(Digit-1))) ...
                            && (Typ_Pri(i) < New_lineP_D) && (i < MM)
                        i = i + 1;
                    end
                    
                    if i > MM
                        TypePosition_F = 100; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose = -1;
                    elseif Typ_Pri(i) <= Typ_Pri(Spark) - 4 * 10^(-(Digit-1))
                        TypePosition_F = 32; % No mode & Risk Free
                        TimeClose(ii-1,1) = i;
                        PriceClose = Typ_Pri(Spark) - 4 * 10^(-(Digit-1));
                    elseif Typ_Pri(i) >= New_lineP_D
                        TypePosition_F = 33; % No mode & Deviation
                        NSZO = NSZO + 1;
                        %SequenceZeroOne = 0;
                        TimeClose(ii-1,1) = i;
                        PriceClose = New_lineP_D;
                    else
                        TypePosition_F = 100; % No mode & ignore
                        TimeClose(ii-1,1) = -1;
                        PriceClose = -1;
                    end
                end
            end
        else
            Hunt_Forecast = {'There is no real data'};
            SellBuy = {'There is no real data'};
            
            TypePosition_F = -1;
            
            TimeClose(ii-1,1) = -1;
            PriceClose = -1;
        end
        %----------------------------------------------
    elseif Typ_Pri(Spark) == PForecast % Hunt Equal
        Hunt_Forecast = {'Hunt_E'};
        SellBuy = {'No'};
%         Max_PercentF(ii-1,1) = -1;
%         Deviation_F(ii-1,1) = -1;
        
        TypePosition_F = -1;
%         TypePos_MaxPer_Dev_F(ii-1,1) = -1;
%         TypePos_F(ii-1,1) = {'NO'};
        
        %FinalPrice_F(ii-1,1) = Typ_Pri(YForecast);
%        PiP(ii-1,1) = -1;
        %NetProfit(ii-1,1) = -1;
        %PerPosition(ii-1,1) = -1;
        
        TimeClose(ii-1,1) = -1;
        PriceClose = -1;
        
%         New_lineP = -1;
%         New_lineP_D = -1;
    end
else
    Hunt_Forecast = {'Position at the same time'};
    SellBuy = {'Position at the same time'};
    
    TypePosition_F = -1;
%         TypePos_MaxPer_Dev_F(ii-1,1) = -1;
%         TypePos_F(ii-1,1) = {'Position at the same time'};
        
        %FinalPrice_F(ii-1,1) = Typ_Pri(YForecast);
        %PiP(ii-1,1) = -1;
        %NetProfit(ii-1,1) = -1;
        %PerPosition(ii-1,1) = -1;
    
    TimeClose(ii-1,1) = -1;
    PriceClose = -1;
    
%         New_lineP = -1;
%         New_lineP_D = -1;
end

end

