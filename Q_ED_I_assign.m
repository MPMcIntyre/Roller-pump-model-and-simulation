function Q_ED_I = Q_ED_I_assign(i,Q_flow,Q_ED_I,s_time,h)

for j = 1:numel(Q_flow)
    %If the value of j+1 is larger than the Q_ED_I array the value
    %assignment can stop
    if (j+i>=round(s_time/h)+1)
        return;
    end
    
    Q_ED_I(i+j) = Q_flow(j);
    j = j+1;
end
