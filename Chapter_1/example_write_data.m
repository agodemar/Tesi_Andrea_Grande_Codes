x = linspace(0,2*pi,20);
y = sin(x);

data_xy = [x', y'];

save("data_xy.txt", 'data_xy', '-ascii')