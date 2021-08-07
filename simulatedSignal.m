%% Signal Simulation
% Ref: An improved complementary ensemble empirical mode decomposition with
% adaptive noise and its application to rolling element bearing fault diagnosis

clc; clear;

for n = 1:1000
    if n >=1 && n <=500
        x1{n} = 0;
    elseif n >= 501 && n <= 750
        x1{n} = sin(2*pi*0.3*(n - 501));
    else
        x1{n} = 0;
    end
    x2{n} = cos(2*pi*0.05*(n - 1));
    x3{n} = (n/1000) + (n/1000).^2;
end

y = (cell2mat(x1)+cell2mat(x2)+cell2mat(x3))';
csvwrite('simulatedSignal.csv',y);
X = {cell2mat(x1),cell2mat(x2),cell2mat(x3), y};

a = cell2mat(x1)';
b = cell2mat(x2)';
c = cell2mat(x3)';

sig = [a, b, c, y];
