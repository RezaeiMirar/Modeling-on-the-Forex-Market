function c = CONT_R_2(Slo_1, m_1)
c = 0;
for i = 1 : m_1
    if Slo_1(i) == 1
        c = c + 1;
    else
        continue
    end
end

end
