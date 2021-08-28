fd = fopen("image1", "r");
img = fread(fd, Inf, "uint16");
img = reshape(img, 501, 501);
img = img';
M = max(max(img));
img = (img/M) * 255;    
w = img(350:430, 100:400);
b = img(200:300, 200:300);
bm = mean(vec(b));
wm = mean(vec(w));
mu = -log(bm/wm) / 6;
tmu = 1.0217;
abs_err = abs((mu -  tmu) / tmu) * 100;
if (abs_err <= 5.0)
    printf("\nOK - Silicon relative error is less than 5%%\n");
else
    printf("\nERROR - Silicon relative error is greater than 5%%\n");
    exit;
endif

