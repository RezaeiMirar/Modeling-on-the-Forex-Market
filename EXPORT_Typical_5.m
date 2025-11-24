function [Ext_Price, Spark, Delay] = EXPORT_Typical_5(M, Typ_Pri, t, EMA_11, EMA_22)

Period1 = 12;
Period2 = 26;

%--------------------------------------------------------------------------
Alpha1 = 1 - 2 / (Period1 + 1);
Beta1 = 2 / (Period1 + 1);

EMA_1(1) = EMA_11; % 12 candels
for i = 2 : M - Period1 + 1
    EMA_1(i) = ( Alpha1 * EMA_1(i-1) + Beta1 * Typ_Pri(Period1+i-1) )/ (Alpha1 + Beta1);
end

%--------------------------------------------------------------------------
Alpha2 = 1 - 2 / (Period2 + 1);
Beta2 = 2 / (Period2 + 1);

EMA_2(1) = EMA_22; % 26 candels
for i = 2 : M - Period2 + 1
    EMA_2(i) = ( Alpha2 * EMA_2(i-1) + Beta2 * Typ_Pri(Period2+i-1) )/ (Alpha2 + Beta2);
end
SIZE_1 = length(EMA_2);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>  The MACD line formula  >>>>>>>>>>>>>>>>>>>>>>
% MACD line = 12-day EMA – 26-day EMA

% EMA_1(15) - EMA_2(1) = 1.079386 - 1.079422;
% Period2 - Period1 + 1 = 26 - 12 + 1 = 15
for i = 1 : SIZE_1
    Macd(i) = EMA_1(Period2-Period1+i) - EMA_2(i); % Start from candle 15 in EMA_1. 
end
SIZE_2 = length(Macd);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% figure()
% plot(t,Typ_Pri,'k-',t(Period1:M),EMA_1,'m-',t(Period2:M),EMA_2,'b-')
% legend('Price',['EMA Period: ',num2str(Period1)],['EMA Period: ',num2str(Period2)])
% %--------------------------------------------------------------------------
% figure()
% bar(t(Period2:M),Macd)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  MACD sign count  >>>>>>>>>>>>>>>>>>>>>>>>>>

SIGN = sign(Macd);
c1 = 0;
if SIGN(1) == 1
    cont(1,1) = CONT(c1, SIZE_2, SIGN); % +
    cont(1,2) = SIGN(1);
elseif SIGN(1) == -1
    cont(1,1) = CONT(c1, SIZE_2, SIGN); % -
    cont(1,2) = SIGN(1);
end
sum = cont(1,1);

if SIGN(SIZE_2) == 1
    ww = CONT_2(SIZE_2, SIGN);
elseif SIGN(SIZE_2) == -1
    ww = CONT_2(SIZE_2, SIGN);
end

i1 = 2;
while sum < SIZE_2 - ww
    if SIGN( sum + 1 ) == 1
        cont(i1,1) = CONT(sum, SIZE_2, SIGN);
        cont(i1,2) = SIGN( sum + 1 );
    elseif SIGN( sum + 1 ) == -1
        cont(i1,1) = CONT(sum, SIZE_2, SIGN);
        cont(i1,2) = SIGN( sum + 1 );
    end
    sum = sum + cont(i1,1);
    i1 = i1 + 1;
end
cont(i1,1) = ww;
cont(i1,2) = SIGN(SIZE_2);
disp(cont)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Change of MACD sign  >>>>>>>>>>>>>>>>>>>>>>>
j1 = 1;
i2 = 1;
if cont(i2,1) < 40  && cont(i2+1,1) >= 40
    Index(j1,1) = cont(i2,1)+cont(i2+1,1);
    Index(j1,2) = cont(i2+1,2);
    Index(j1,3) = cont(i2,1)+cont(i2+1,1);
    Index(j1,4) = i2 + 1;
    Index(j1,5) = 1001;
    %Index(j1,7) = Index(j1,3); % Spark Time
elseif cont(i2,1) >= 40
    Index(j1,1) = cont(i2,1);
    Index(j1,2) = cont(i2,2);
    Index(j1,3) = cont(i2,1);
    Index(j1,4) = i2;
    Index(j1,5) = 1002;
    %Index(j1,7) = Index(j1,3); % Spark Time
elseif cont(i2,1) < 40  && cont(i2+1,1) < 40
    s1 = cont(1,1); s2 = 0; s3 = cont(1,1); Diff = 0;
    while ( cont(i2+1,1) < 40 ) && i2 < size(cont,1)-1 && ~(Diff>=40 || Diff<=-40)
        if cont(i2+1,2) == cont(1,2)
            s1 = s1 + cont(i2+1,1);
        elseif cont(i2+1,2) ~= cont(1,2)
            s2 = s2 + cont(i2+1,1);
        end
        Diff = s1 - s2;
        DIFF(i2,1) = Diff;
        s3 = s3 + cont(i2+1,1);
        i2 = i2 + 1;
    end
    if Diff >= 40
        Index(j1,1) = s3;
        Index(j1,2) = cont(1,2);
        Index(j1,3) = s3;
        Index(j1,4) = i2;
        Index(j1,5) = 1003;
        Index(j1,6) = 1;
        %Index(j1,7) = Index(j1,3); % Spark Time
    elseif Diff <= -40
        Index(j1,1) = s3;
        Index(j1,2) = -cont(1,2);
        Index(j1,3) = s3;
        Index(j1,4) = i2;
        Index(j1,5) = 1004;
        Index(j1,6) = 1;
        %Index(j1,7) = Index(j1,3); % Spark Time
    else
        Index(j1,1) = s3 + cont(i2+1,1);
        Index(j1,2) = cont(i2+1,2);
        Index(j1,3) = s3 + cont(i2+1,1);
        Index(j1,4) = i2 + 1;
        Index(j1,5) = 1005;
        %Index(j1,7) = Index(j1,3); % Spark Time
    end
end
disp('Index =')
disp(Index)

j1 = j1 + 1;
i2 = i2 + 1;
while i2 < length(cont(:,1))-1
if cont(i2,1) >= 40 && cont(i2,2) ~= Index(j1-1,2)
    Index(j1,2) = cont(i2,2);
    Index(j1,3) = Index(j1-1,3) + cont(i2,1);
    Index(j1,1) = Index(j1,3) - Index(j1-1,3);
    Index(j1,4) = i2;
    Index(j1,5) = 2001;
    Index(j1,7) = Index(j1-1,3) + 40;
    i2 = i2 + 1;
    j1 = j1 + 1;
elseif cont(i2,1) >= 40  && cont(i2+1,1) < 40
    s1 = 0; s2 = 0; s3 = 0; Diff = 0;
    while ( ~(cont(i2+1,1) >= 40 && cont(i2+1,2)~=Index(j1-1,2)) ) ...
            && i2 < size(cont,1)-1
        if cont(i2+1,2) == Index(j1-1,2)
            s1 = s1 + cont(i2+1,1);
        elseif cont(i2+1,2) ~= Index(j1-1,2)
            s2 = s2 + cont(i2+1,1);
        end
        Diff = s1 - s2;
        DIFF(i2,1) = Diff;
        s3 = s3 + cont(i2+1,1);
        if Diff <= -40
            i2 = i2 + 2;
            break
        else
            i2 = i2 + 1;
        end
    end
    if Diff <= -40
        Index(j1,1) = s3;
        Index(j1,2) = -Index(j1-1,2);
        Index(j1,3) = Index(j1-1,3) + s3;
        Index(j1,4) = i2-1;
        Index(j1,5) = 2002; % Additional
        Index(j1,6) = 1; % Additional
        Index(j1,7) = Index(j1,3) - ( abs(Diff) - 40 );
        j1 = j1 + 1;
    else
        Index(j1-1,1) = Index(j1-1,1) + s3;
        Index(j1-1,3) = Index(j1-1,3) + s3;
        Index(j1-1,4) = i2;
        Index(j1-1,5) = 2003;
     end
elseif cont(i2+1,1) >= 40 && cont(i2+1,2) ~= Index(j1-1,2)
    Index(j1,2) = cont(i2+1,2);
    Index(j1,3) = Index(j1-1,3) + cont(i2+1,1);
    Index(j1,1) = Index(j1,3) - Index(j1-1,3);
    Index(j1,4) = i2 + 1;
    Index(j1,5) = 2004;
    Index(j1,7) = Index(j1-1,3) + 40;
    i2 = i2 + 1;
    j1 = j1 + 1;
elseif cont(i2-1,1) >= 40  && cont(i2,1) < 40
    s1 = 0; s2 = 0; s3 = 0; Diff = 0;
    while ( ~(cont(i2,1) >= 40 && cont(i2,2)~=Index(j1-1,2)) ) ...
            && i2 <= size(cont,1)-1
        if cont(i2,2) == Index(j1-1,2)
            s1 = s1 + cont(i2,1);
        elseif cont(i2,2) ~= Index(j1-1,2)
            s2 = s2 + cont(i2,1);
        end
        Diff = s1 - s2;
        DIFF(i2,1) = Diff;
        s3 = s3 + cont(i2,1);
        if Diff <= -40
            i2 = i2 + 1;
            break
        else
            i2 = i2 + 1;
        end
    end
    if Diff <= -40
        Index(j1,1) = s3;
        Index(j1,2) = -Index(j1-1,2);
        Index(j1,3) = Index(j1-1,3) + s3;
        Index(j1,4) = i2-1;
        Index(j1,5) = 2005;
        Index(j1,6) = 1;
        Index(j1,7) = Index(j1,3) - ( abs(Diff) - 40 );
        j1 = j1 + 1;
    else
        Index(j1-1,1) = Index(j1-1,1) + s3;
        Index(j1-1,3) = Index(j1-1,3) + s3;
        Index(j1-1,4) = i2-1;
        Index(j1-1,5) = 2006;
        i2 = i2 - 1;
    end
elseif cont(i2-1,1) < 40  && cont(i2,1) < 40
    s1 = 0; s2 = 0; s3 = 0; Diff = 0;
    while ( ~(cont(i2,1) >= 40 && cont(i2,2)~=Index(j1-1,2)) ) ...
            && i2 <= size(cont,1)-1
        if cont(i2,2) == Index(j1-1,2)
            s1 = s1 + cont(i2,1);
        elseif cont(i2,2) ~= Index(j1-1,2)
            s2 = s2 + cont(i2,1);
        end
        Diff = s1 - s2;
        DIFF(i2,1) = Diff;
        s3 = s3 + cont(i2,1);
        if Diff <= -40
            i2 = i2 + 1;
            break
        else
            i2 = i2 + 1;
        end
    end
    if Diff <= -40
        Index(j1,1) = s3;
        Index(j1,2) = -Index(j1-1,2);
        Index(j1,3) = Index(j1-1,3) + s3;
        Index(j1,4) = i2-1;
        Index(j1,5) = 2007;
        Index(j1,6) = 1;
        Index(j1,7) = Index(j1,3) - ( abs(Diff) - 40 );
        j1 = j1 + 1;
    else
        Index(j1-1,1) = Index(j1-1,1) + s3;
        Index(j1-1,3) = Index(j1-1,3) + s3;
        Index(j1-1,4) = i2-1;
        Index(j1-1,5) = 2008;
        i2 = i2 - 1;
    end

end

end

disp('Index =')
disp(Index)
if size(cont,1) == Index(j1-1,4) + 2 % i2 == Index(j1-1,4) + 1
    if cont(i2,1) < 40 && cont(i2+1,1) < 40
        Index(j1-1,1) = Index(j1-1,1) + cont(i2,1) + cont(i2+1,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2,1) + cont(i2+1,1);
        Index(j1-1,4) = i2+1;
        Index(j1-1,5) = 5007;
    elseif cont(i2,1) >= 40 && cont(i2+1,1) < 40 && cont(i2,2) == Index(j1-1,2)
        Index(j1-1,1) = Index(j1-1,1) + cont(i2,1) + cont(i2+1,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2,1) + cont(i2+1,1);
        Index(j1-1,4) = i2+1;
        Index(j1-1,5) = 5008;
    elseif cont(i2,1) < 40 && cont(i2+1,1) >= 40 && cont(i2+1,2) == Index(j1-1,2)
        Index(j1-1,1) = Index(j1-1,1) + cont(i2,1) + cont(i2+1,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2,1) + cont(i2+1,1);
        Index(j1-1,4) = i2+1;
        Index(j1-1,5) = 5009;
    elseif cont(i2,1) >= 40 && cont(i2+1,1) < 40 && cont(i2,2) == -Index(j1-1,2)
        Index(j1,1) = cont(i2,1) + cont(i2+1,1);
        Index(j1,2) = cont(i2,2);
        Index(j1,3) = Index(j1-1,3) + cont(i2,1) + cont(i2+1,1);
        Index(j1,4) = i2+1;
        Index(j1,5) = 5010;
        Index(j1,7) = Index(j1-1,3) + 40;
    elseif cont(i2,1) < 40 && cont(i2+1,1) >= 40 && cont(i2+1,2) == -Index(j1-1,2)
        Index(j1-1,1) = Index(j1-1,1) + cont(i2,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2,1);
        Index(j1-1,4) = i2;
        Index(j1-1,5) = 5011;
        
        Index(j1,1) = cont(i2+1,1);
        Index(j1,2) = cont(i2+1,2);
        Index(j1,3) = Index(j1-1,3) + cont(i2+1,1);
        Index(j1,4) = i2+1;
        Index(j1,5) = 5012;
        Index(j1,7) = Index(j1-1,3) + 40;
    elseif cont(i2,1) >= 40 && cont(i2+1,1) >= 40 && cont(i2,2) == Index(j1-1,2)
        Index(j1-1,1) = Index(j1-1,1) + cont(i2,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2,1);
        Index(j1-1,4) = i2;
        Index(j1-1,5) = 5015;
        
        Index(j1,1) = cont(i2+1,1);
        Index(j1,2) = cont(i2+1,2);
        Index(j1,3) = Index(j1-1,3) + cont(i2+1,1);
        Index(j1,4) = i2+1;
        Index(j1,5) = 5016;
        Index(j1,7) = Index(j1-1,3) + 40;
    elseif cont(i2,1) >= 40 && cont(i2+1,1) >= 40 && cont(i2,2) == -Index(j1-1,2)
        Index(j1,1) = cont(i2,1);
        Index(j1,2) = cont(i2,2);
        Index(j1,3) = Index(j1-1,3) + cont(i2,1);
        Index(j1,4) = i2;
        Index(j1,5) = 5013;
        Index(j1,7) = Index(j1-1,3) + 40;
        
        Index(j1+1,1) = cont(i2+1,1);
        Index(j1+1,2) = cont(i2+1,2);
        Index(j1+1,3) = Index(j1,3) + cont(i2+1,1);
        Index(j1+1,4) = i2+1;
        Index(j1+1,5) = 5014;
        Index(j1+1,7) = Index(j1,3) + 40;
    end
elseif i2 == size(cont,1) - 1
    if cont(i2+1,1) < 40
        Index(j1-1,1) = Index(j1-1,1) + cont(i2+1,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2+1,1);
        Index(j1-1,4) = i2+1;
        Index(j1-1,5) = 5001;
    elseif cont(i2+1,1) >= 40 && cont(i2+1,2) == Index(j1-1,2)
        Index(j1-1,1) = Index(j1-1,1) + cont(i2+1,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2+1,1);
        Index(j1-1,4) = i2+1;
        Index(j1-1,5) = 5002;
    elseif cont(i2+1,1) >= 40 && cont(i2+1,2) == -Index(j1-1,2)
        Index(j1,1) = cont(i2+1,1);
        Index(j1,2) = cont(i2+1,2);
        Index(j1,3) = Index(j1-1,3) + cont(i2+1,1);
        Index(j1,4) = i2+1;
        Index(j1,5) = 5003;
        Index(j1,7) = Index(j1-1,3) + 40;
    end
elseif i2 == size(cont,1)
    if cont(i2,1) < 40
        Index(j1-1,1) = Index(j1-1,1) + cont(i2,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2,1);
        Index(j1-1,4) = i2;
        Index(j1-1,5) = 5004;
    elseif cont(i2,1) >= 40 && cont(i2,2) == Index(j1-1,2)
        Index(j1-1,1) = Index(j1-1,1) + cont(i2,1);
        Index(j1-1,3) = Index(j1-1,3) + cont(i2,1);
        Index(j1-1,4) = i2;
        Index(j1-1,5) = 5005;
    elseif cont(i2,1) >= 40 && cont(i2,2) == -Index(j1-1,2)
        Index(j1,1) = cont(i2,1);
        Index(j1,2) = cont(i2,2);
        Index(j1,3) = Index(j1-1,3) + cont(i2,1);
        Index(j1,4) = i2;
        Index(j1,5) = 5006;
        Index(j1,7) = Index(j1-1,3) + 40;
    end
end

disp('Index =')
disp(Index)

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Extreme prices  >>>>>>>>>>>>>>>>>>>>>>>>
if Index(1,2) == 1
    Ext_Price(1,1) = max( Typ_Pri(1:Index(1,3)+Period2-1) );
    ind = find( Typ_Pri(1:Index(1,3)+Period2-1) == Ext_Price(1,1) );
    Ext_Price(1,2) = ind(1);
    Ext_Price(1,3) = Ext_Price(1,2);
    if length(ind) > 1
        Ext_Price(1,4) = 1; % Additional
    else
        Ext_Price(1,4) = 0; % Additional
    end
else
    Ext_Price(1,1) = min( Typ_Pri(1:Index(1,3)+Period2-1) );
    ind = find( Typ_Pri(1:Index(1,3)+Period2-1) == Ext_Price(1,1) );
    Ext_Price(1,2) = ind(1);
    Ext_Price(1,3) = Ext_Price(1,2);
    if length(ind) > 1
        Ext_Price(1,4) = 1; % Additional
    else
        Ext_Price(1,4) = 0; % Additional
    end
end
for i = 2 : size(Index,1)
    if Index(i,2) == 1
        Ext_Price(i,1) = max( Typ_Pri(Index(i-1,3)+Period2:Index(i,3)+Period2-1) );
        ind = find( Typ_Pri(Index(i-1,3)+Period2:Index(i,3)+Period2-1) == Ext_Price(i,1) );
        Ext_Price(i,2) = ind(1);
        Ext_Price(i,3) = Index(i-1,3) + Period2 + Ext_Price(i,2) - 1;
        if length(ind) > 1
            Ext_Price(i,4) = 1; % Additional
        else
            Ext_Price(i,4) = 0; % Additional
        end
    else
        Ext_Price(i,1) = min( Typ_Pri(Index(i-1,3)+Period2:Index(i,3)+Period2-1) );
        ind = find( Typ_Pri(Index(i-1,3)+Period2:Index(i,3)+Period2-1) == Ext_Price(i,1) );
        Ext_Price(i,2) = ind(1);
        Ext_Price(i,3) = Index(i-1,3) + Period2 + Ext_Price(i,2) - 1;
        if length(ind) > 1
            Ext_Price(i,4) = 1; % Additional
        else
            Ext_Price(i,4) = 0; % Additional
        end
    end
end
disp(' Extreme prices')
disp(Ext_Price)

%--------------------------------------------------------------------------
i3 = 0;
SHOW = NaN(1); % Additional
for i = 1 : length(DIFF)
    if DIFF(i,1) <= -40
        i3 = i3 + 1;
        SHOW(i3,1) = i;
        SHOW(i3,2) = DIFF(i,1);
    else
        continue
    end
end
disp(' Index   Negetive')
disp(SHOW)

for i = 1 : size(Index,1) - 1
    Spark(i,1) = Index(i+1,7) + Period2 - 1;
    Delay(i,1) = Index(i+1,7) - Index(i,3);
end
Spark(size(Index,1),1) = 0;
Delay(size(Index,1),1) = 0;

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% figure()
% subplot(2,1,1)
% plot(t,Typ_Pri,'k-',t(Period1:M),EMA_1,'m-',t(Period2:M),EMA_2,'b-')
% axis tight
% Y_1 = ylim; % Y = [min(A(1:M,4)), max(A(1:M,4))];
% hold on
% plot(Ext_Price(:,3),Ext_Price(:,1),'go','linewidth',2)
% legend('Price',['EMA Period: ',num2str(Period1)],['EMA Period: ',num2str(Period2)])
% 
% for i = 1 : size(Index,1) % Additional
%     hold on
%     X_1 = [Ext_Price(i,3), Ext_Price(i,3)];
%     line(X_1,Y_1,'Color',[.6 .1 .9])
%     X_2 = [Index(i,3)+Period2-1, Index(i,3)+Period2-1];
%     line(X_2,Y_1,'Color',[.3 .9 .9])
% end
% 
% subplot(2,1,2)
% bar(t,[zeros(1,25),Macd]) % bar(t(Period2:M),Macd)
% axis tight
% Y_2 = ylim;
% for i = 1 : size(Index,1) % Additional
%     hold on
%     X_1 = [Ext_Price(i,3), Ext_Price(i,3)];
%     line(X_1,Y_2,'Color',[.6 .1 .9])
%     X_2 = [Index(i,3)+Period2-1, Index(i,3)+Period2-1];
%     line(X_2,Y_2,'Color',[.3 .9 .9])
% end
% %--------------------------------------------------------------------------
% figure()
% plot(t,Typ_Pri,'k')
% hold on
% plot(Ext_Price(:,3),Ext_Price(:,1),'mo','linewidth',2)
% legend('Price','Extreme')

end

