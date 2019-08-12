function red_data = compute_reduced_data(m,d)

for i=1:size(d.AII1,2)
    red_data.ANII1{i} = d.RB1'*d.AII1{i}*d.RB1;
    red_data.ANGI1{i} = d.RBG1'*d.AGI1{i}*d.RB1;
    red_data.ANIG1{i} = d.RB1'*d.AIG1{i}*d.RBG1;
    red_data.ANGG1{i} = d.RBG1'*d.AGG1{i}*d.RBG1;
end



for i=1:size(d.AII2,2)
    red_data.ANII2{i} = d.RB2'*d.AII2{i}*d.RB2;
    red_data.ANGI2{i} = d.RBG2'*d.AGI2{i}*d.RB2;
    red_data.ANIG2{i} = d.RB2'*d.AIG2{i}*d.RBG2;
    red_data.ANGG2{i} = d.RBG2'*d.AGG2{i}*d.RBG2;
end

red_data.LNGG1 = d.RBG1'*d.LGG1*d.RBG1;
red_data.LNGG2 = d.RBG2'*d.LGG2*d.RBG2;

red_data.LNGG1_cross = d.RBG1'*d.LGG1*d.RBG2;
red_data.LNGG2_cross = d.RBG2'*d.LGG2*d.RBG1;

red_data.W1 = d.RB1'*d.W1*d.RB1;
red_data.W2 = d.RB2'*d.W2*d.RB2;

%% Now right hand side
for i=1:size(d.r1,2)
    red_data.rN1{i} = d.RB1'*d.r1{i};
    red_data.rNG1{i} = d.RBG1'*d.rG1{i};
end
for i=1:size(d.r2,2)
    red_data.rN2{i} = d.RB2'*d.r2{i};
    red_data.rNG2{i} = d.RBG2'*d.rG2{i};
end

end