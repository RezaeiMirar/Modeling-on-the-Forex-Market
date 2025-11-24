function RR = RR_TP_SL(Sequence, Reward_Risk)

Num_Profit = 0; % TP
Num_Loss = 0; % SL
for i = 1 : length(Sequence)
   if Sequence(i) == 1
       Num_Profit = Num_Profit + 1;
   elseif (Sequence(i) == 2) || (Sequence(i) == 33)
       Num_Loss = Num_Loss + 1;
   end
end

RR = Reward_Risk * Num_Profit - Num_Loss;

end

% clear all
% clc
% d=[-1,-1,2,33,1,2,1,-1];
% Num_Profit = 0;
% Num_Loss = 0;
% for i = 1 : length(d)
%    if d(i) == 1
%        Num_Profit = Num_Profit + 1;
%    elseif (d(i) == 2) || (d(i) == 33)
%        Num_Loss = Num_Loss + 1;
%    end
% end
% disp(Num_Profit)
% disp(Num_Loss)
% [Num_Profit, Num_Loss] = RR_TP_SL(d, 1);
% disp(Num_Profit)
% disp(Num_Loss)

