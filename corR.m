function cor = corR(original, signal)
cor = {};
for i = 1: size(signal,2)
    cor{i} = corrcoef(original, signal(:,i));
    cor{i} = sum(cor{1,i}(:,1))-1;
end
cor = cell2mat(cor);
end