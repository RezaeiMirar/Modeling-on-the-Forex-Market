function c = CONT(N, SIZE_2, SIGN)
c = 0;
for i = N + 1 : SIZE_2
    if SIGN(i) ~= SIGN(i+1)
        c = c + 1;
        break
    else
        c = c + 1;
    end
end

end
