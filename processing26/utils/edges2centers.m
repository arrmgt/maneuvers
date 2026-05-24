function CENTERS=edges2centers(EDGES);
windowSize = 2; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
CENTERS=filter(b,a,EDGES);
CENTERS=CENTERS(2:end);