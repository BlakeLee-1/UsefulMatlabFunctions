function test()
x = -pi:0.1:pi;
y = cos(x.*4)+1.5;

%polarplot(x(1:floor(length(x)/2)),y(1:floor(length(x)/2)))
polarplot(x,y)

