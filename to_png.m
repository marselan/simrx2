images = glob("image*");
for i = 1:length(images)
    m = fopen(strcat("results/image", num2str(i)), "r");
    v = fread(m, Inf, "uint16");
    v = reshape(v, 501, 501);
    v = v';
    M = max(max(v));
    v = (v/M) * 255;
    w = uint8(v);
    imwrite(w, strcat('results/png/image', num2str(i), '.png'));
endfor
