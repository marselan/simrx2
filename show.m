
clear;
m = Inf;
M = -Inf;
for i=1:360
%figure
fp = fopen(strcat('detector',num2str(i)), 'r');
I = fread(fp, [501 501], 'uint32');
fclose(fp);
m = min(m, min(min(I)));
M = max(M, max(max(I)));
%imshow(I, []);
end
m
M






