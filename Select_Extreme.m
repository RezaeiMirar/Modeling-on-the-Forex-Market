function S_E = Select_Extreme(y, Index_Begin, Index_End) % Select extreme for curve
i1 = 0;
for i = ( Index_Begin : Index_End - 2 ) - ( Index_Begin - 1 )
    if ( y(i+1)<y(i) && y(i+1)<y(i+2) )
        i1 = i1 + 1;
        S_E(i1,1) = i + 1 + ( Index_Begin - 1 );
        S_E(i1,2) = y(i+1);
        S_E(i1,3) = 1; % min
    elseif ( y(i+1)>y(i) && y(i+1)>y(i+2) )
        i1 = i1 + 1;
        S_E(i1,1) = i + 1 + ( Index_Begin - 1 );
        S_E(i1,2) = y(i+1);
        S_E(i1,3) = 2; % max
    else
        continue
    end
end
if i1 == 0
    i1 = i1 + 1;
    S_E(i1,1) = 0;
    S_E(i1,2) = 0;
    S_E(i1,3) = 0; % max or min
end

end
