function Q_ED_II = Q_ED_II_assign(i,Q_flow2,Q_ED_II,s_time,h)

for j = 1:numel(Q_flow2)
    %If the value of j+1 is larger than the Q_ED_I array the value
    %assignment can stop
    if (j+i>=round(s_time/h)+1)
        return;
    end
    Q_ED_II(i+j) = Q_flow2(j);
    j = j+1;   
end
