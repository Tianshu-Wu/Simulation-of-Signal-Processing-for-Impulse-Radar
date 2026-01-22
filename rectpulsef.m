%% 矩形窗函数
function y = rectpulsef(t,Tw)
y = zeros(size(t));
index = find(abs(t)<=Tw/2);
y(index) = 1;
end