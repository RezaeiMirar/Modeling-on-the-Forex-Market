function c = CONT_2(SIZE_2, SIGN)
c = 0;
for i = SIZE_2 : -1 : 1
    if SIGN(i) ~= SIGN(i-1)
        c = c + 1;
        break
    else
        c = c + 1;
    end
end

end
