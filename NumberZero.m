function IgnoreSequence = NumberZero(Sequence)

IgnoreSequence = 1;
num = 0;
for i = 1 : length(Sequence)
    if Sequence(i) == 0
        num = num + 1;
    else
        num = 0;
    end
    
    if num > 5
        IgnoreSequence = 0;
        break
    end
end
end

% clear all
% clc
% a=[0,0,0,0,0,0,1,0,0,0,0,0,0,1,0];
% IgnoreSequence = 1;
% num = 0;
% for i = 1 : length(a)
%     b(i)=a(i)
%     if a(i) == 0
%         num = num + 1
%     else
%         num = 0;
%     end
%     
%     if num > 5
%         IgnoreSequence = 0
%         break
%     end
% end
% disp(IgnoreSequence)


