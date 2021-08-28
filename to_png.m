for i = 1:1
    m = fopen(strcat("result/image", num2str(i)), "r");
    v = fread(m, Inf, "uint16");
    v = reshape(v, 501, 501);
    v = v';
    M = max(max(v));
    v = (v/M) * 255;
    w = uint8(v);
    imwrite(w, strcat('result/png/image', num2str(i), '.png'));
endfor
