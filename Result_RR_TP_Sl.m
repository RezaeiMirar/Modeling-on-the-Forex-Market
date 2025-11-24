function Result = Result_RR_TP_Sl(Sequence_Sel_R1, Sequence_Rev_R1, Sequence_Sel_R2, ...
    Sequence_Rev_R2, Sequence_Sel_R3, Sequence_Rev_R3, ...
    TypePosition_Sel_R1, TypePosition_Rev_R1, TypePosition_Sel_R2, TypePosition_Rev_R2, ...
    TypePosition_Sel_R3, TypePosition_Rev_R3)

IgnoreSequence_Sel_R1 = NumberZero(Sequence_Sel_R1);
IgnoreSequence_Rev_R1 = NumberZero(Sequence_Rev_R1);

IgnoreSequence_Sel_R2 = NumberZero(Sequence_Sel_R2);
IgnoreSequence_Rev_R2 = NumberZero(Sequence_Rev_R2);

IgnoreSequence_Sel_R3 = NumberZero(Sequence_Sel_R3);
IgnoreSequence_Rev_R3 = NumberZero(Sequence_Rev_R3);

if (IgnoreSequence_Sel_R1 == 1) && (IgnoreSequence_Rev_R1 == 1) && ...
        (IgnoreSequence_Sel_R2 == 1) && (IgnoreSequence_Rev_R2 == 1) && ...
        (IgnoreSequence_Sel_R3 == 1) && (IgnoreSequence_Rev_R3 == 1)
    RR(1) = RR_TP_SL(TypePosition_Sel_R1, 1);
    RR(2) = RR_TP_SL(TypePosition_Rev_R1, 1);
    RR(3) = RR_TP_SL(TypePosition_Sel_R2, 1);
    RR(4) = RR_TP_SL(TypePosition_Rev_R2, 1);
    RR(5) = RR_TP_SL(TypePosition_Sel_R3, 1);
    RR(6) = RR_TP_SL(TypePosition_Rev_R3, 1);
    RR_max = max(RR);
    RR_Index = find(RR == RR_max);
    if RR_Index == 1
        Result = {'Direct_100_100'};
    elseif RR_Index == 2
        Result = {'Reverse_100_100'};
    elseif RR_Index == 3
        Result = {'Direct_50_50'};
    elseif RR_Index == 4
        Result = {'Reverse_50_50'};
    elseif RR_Index == 5
        Result = {'Direct_150_50'};
    elseif RR_Index == 6
        Result = {'Reverse_150_50'};
    end
end




end
