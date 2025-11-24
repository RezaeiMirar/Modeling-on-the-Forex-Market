function [IgnoreSequence,MaxZeroCunt] = NumberZeroCont(Sequence)
IgnoreSequence = 1;
num = 0;
i = 1;
i2 = 1;
d(i2) = 0;

if Sequence(1) == 1
    while Sequence(i) == 1
        i = i + 1;
    end
    while i <= length(Sequence)
        if Sequence(i) == 0
            num = num + 1;
            d(i2) = num;
            i = i + 1;
        else
            num = 0;
            while (i <= length(Sequence)) && (Sequence(i) == 1 )
                i = i + 1;
            end
            i2 = i2 +1;
        end
        if num > 5
            IgnoreSequence = 0;
            break
        end
    end
else
    while i <= length(Sequence)
        if Sequence(i) == 0
            num = num + 1;
            d(i2) = num;
            i = i + 1;
        else
            num = 0;
            while (i <= length(Sequence)) && (Sequence(i) == 1 )
                i = i + 1;
            end
            i2 = i2 +1;
        end
        if num > 5
            IgnoreSequence = 0;
            break
        end
    end
end
d
MaxZeroCunt = max(d);

end


% clear all
% clc
% %a=[0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0];
% %a=[0,0,0,0,1,0,0,0,0,0,1,0];
% %a=[0,0,1,1,1,0,0,0,0,0,1,0];
% %a=[1,1,1,1,0,0,0,0,1,0,0,0,0,0,1,0];
% %a=[1,0,0,0,0,1,0,0,0,0,0,1,0];
% a=[1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0];
% IgnoreSequence = 1;
% num = 0;
% i = 1;
% i2 = 1;
% d(i2) = 0;
% 
% if a(1) == 1
%     while a(i) == 1
%         %b(i)=a(i)
%         i = i + 1;
%     end
%     while i <= length(a)
%         %b(i)=a(i)
%         if a(i) == 0
%             num = num + 1;
%             d(i2) = num;
%             i = i + 1;
%         else
%             num = 0;
%             while a(i) == 1
%                 %b(i)=a(i)
%                 i = i + 1;
%             end
%             i2 = i2 +1;
%         end
%         if num > 5
%             IgnoreSequence = 0;
%             break
%         end
%     end
% else
%     %i=1
%     while i <= length(a)
%         %b(i)=a(i)
%         if a(i) == 0
%             num = num + 1;
%             d(i2) = num;
%             i = i + 1;
%         else
%             num = 0;
%             while a(i) == 1
%                 %b(i)=a(i)
%                 i = i + 1;
%             end
%             i2 = i2 +1;
%         end
%         if num > 5
%             IgnoreSequence = 0;
%             break
%         end
%     end
% end
% disp(IgnoreSequence)
% disp(max(d))



