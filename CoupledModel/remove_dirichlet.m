function A=remove_dirichlet(A,mode,I)
if mode == 0
    A(I,:) = [];
    A(:,I)= [];
else
    for i = 1:length(A)
        A{1,i}(I,:) = [];
        A{1,i}(:,I) = [];
    end
end
end