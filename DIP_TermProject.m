clc;
clear;
close all;

img_s = imresize(imread('source.jpg'),[300 600]);
img_t = imresize(imread('target.jpg'),[400 700]);
[xs,ys] = size(img_s(:,:,1));
[xt,yt] = size(img_t(:,:,1));
lab_s = rgb2lab(img_s);
lab_t = rgb2lab(img_t);

ls = lab_s(:,:,1);
as = lab_s(:,:,2);
bs = lab_s(:,:,3);
lt = lab_t(:,:,1);
at = lab_t(:,:,2);
bt = lab_t(:,:,3);

%%%%% l
max_ls = round(max(ls(:))); min_ls = round(min(ls(:)));
max_lt = round(max(lt(:))); min_lt = round(min(lt(:)));
V = 1;
Bls = round((max_ls - min_ls)/V);
Blt = round((max_lt - min_lt)/V);
hls = zeros([1 Bls]);
hlt = zeros([1 Blt]);
vs = zeros([1 Bls]);
vt = zeros([1 Blt]);
Pls = zeros([xs ys]);
Plt = zeros([xt yt]);
%histogram
for i = 1:Bls
    for n = 1:xs
        for m = 1:ys
            if i == round(((ls(n,m) - min_ls)/V)+1)
                Pls(n,m) = 1;
                hls(i) = hls(i) + Pls(n,m);
            else
                Pls(n,m) = 0;
                hls(i) = hls(i) + Pls(n,m);
            end
        end
    end
end
for i = 1:Blt
    for n = 1:xt
        for m = 1:yt
            if i == round(((lt(n,m) - min_lt)/V)+1)
                Plt(n,m) = 1;
                hlt(i) = hlt(i) + Plt(n,m);
            else
                Plt(n,m) = 0;
                hlt(i) = hlt(i) + Plt(n,m);
            end
        end
    end
end

%scale
k = 2;
syms t;
expr = log2(t);
smax_l = round(eval(subs(expr, t, Bls/10)));
Blk = round(Bls*2^(k-smax_l));


hls_k = hls;
hlt_k = hlt;
% % gradient
d_hlt_k = zeros([1 size(hlt_k,2)-1]);
d2_hlt_k = zeros([1 size(d_hlt_k,2)-1]);
for i=1:size(hlt_k,2)-1
    d_hlt_k(i) = hlt_k(i) - hlt_k(i+1);
end
for i=1:size(d_hlt_k,2)-1
    d2_hlt_k(i) = d_hlt_k(i) - d_hlt_k(i+1);
end
Rmin_k = zeros([1 size(hlt_k,2)]);
for i=1:size(d2_hlt_k,2)
    if (d_hlt_k(i)*d_hlt_k(i+1) < 0) && (d2_hlt_k(i) < 0)
        Rmin_k(i+1) = 1;
    end
end
count = 0;
for i=1:size(Rmin_k,2)
    if Rmin_k(1,i) == 1
        count = count + 1;
    end
end
j = zeros([1,count]);
z = 1;
for i=1:size(Rmin_k,2)
    if Rmin_k(i) == 1
        j(z) = i;
        z = z + 1;
    end
end
a = j(end-1);
b = j(end)-1;

sub_hls_k = hls_k(a:b);
meanls_k = mean(sub_hls_k);
stdls_k = std(sub_hls_k);

sub_hlt_k = hlt_k(a:b);
meanlt_k = mean(sub_hlt_k);
stdlt_k = std(sub_hlt_k);

wsk = k/smax_l;
hlo_k = zeros([1,size(hls_k,2)]);

for i = 1:size(hls_k,2)
    hlo_k(i) = (hls_k(i) - wsk*meanls_k)*(((1-wsk)*stdlt_k)/(wsk*stdls_k)) + (1-wsk)*meanlt_k;
end

Cs = zeros([1,size(hls_k,2)]);
Co = zeros([1,size(hls_k,2)]);
sum_cs = 0;
sum_co = 0;
for j=1:size(hls_k,2)
    for i=1:j
        sum_cs = sum_cs + hls_k(i);
        sum_co = sum_co + hlo_k(i);
    end
    Cs(j) = sum_cs;
    Co(j) = sum_co;
    sum_cs = 0;
    sum_co = 0;
end

Ilo = zeros([xs,ys]);

for n = 1:xs
    for m = 1:ys
        index_cs = round(((ls(n,m) - min_ls+1)/V));
        if index_cs > size(Cs,2)
            index_cs =  size(Cs,2);
        else
            index_cs = index_cs;
        end
        [~, index_co] = min(abs(Co - Cs(index_cs)));
        Ilo(n,m) = index_co ; 
    end
end


Ilo = 4.5*(Ilo + min_ls) ;

max_as = round(max(as(:))); min_as = round(min(as(:)));
max_at = round(max(at(:))); min_at = round(min(at(:)));

Bas = round((max_as - min_as)/V);
Bat = round((max_at - min_at)/V);
has = zeros([1 Bas]);
hat = zeros([1 Bat]);
vs = zeros([1 Bas]);
vt = zeros([1 Bat]);
Pas = zeros([xs ys]);
Pat = zeros([xt yt]);
%histogram
for i = 1:Bas
    v(i) = min_as + (i-1)*V;
    for n = 1:xs
        for m = 1:ys
            if i == round(((as(n,m) - min_as)/V)+1)
                Pas(n,m) = 1;
                has(i) = has(i) + 1;
            end
        end
    end
end
for i = 1:Bat
    for n = 1:xt
        for m = 1:yt
            if i == round(((at(n,m) - min_at)/V)+1)
                Pat(n,m) = 1;
                hat(i) = hat(i) + 1;

            end
        end
    end
end


%Bk
syms t;
expr = log2(t);
smax_a = round(eval(subs(expr, t, Bas/10)));
Bak = round(Bas*2^(k-smax_a));

                                                                                                                      
has_k = has;
hat_k = hat;

d_hat_k = zeros([1 size(hat_k,2)-1]);
d2_hat_k = zeros([1 size(d_hat_k,2)-1]);
for i=1:size(hat_k,2)-1
    d_hat_k(i) = hat_k(i) - hat_k(i+1);
end
for i=1:size(d_hat_k,2)-1
    d2_hat_k(i) = d_hat_k(i) - d_hat_k(i+1);
end
Rmina_k = zeros([1 size(hat_k,2)]);
for i=1:size(d2_hat_k,2)
    if (d_hat_k(i)*d_hat_k(i+1) < 0) && (d2_hat_k(i) < 0)
        Rmina_k(i+1) = 1;
    end
end
count = 0;
for i=1:size(Rmina_k,2)
    if Rmina_k(i) == 1
        count = count + 1;
    end
end
j = zeros([1,count]);
z = 1;
for i=1:size(Rmina_k,2)
    if Rmina_k(i) == 1
        j(z) = i;
        z = z + 1;
    end
end
a = j(1);
b = j(2)-1;

sub_has_k = has_k(a:b);
meanas_k = mean(sub_has_k);
stdas_k = std(sub_has_k);

sub_hat_k = hat_k(a:b);
meanat_k = mean(sub_hat_k);
stdat_k = std(sub_hat_k);

wask = k/smax_a;
hao_k = zeros([1,size(has_k,2)]);

for i = 1:size(has_k,2)
    hao_k(i) = (has_k(i) - wask*meanas_k)*(((1-wask)*stdat_k)/(wask*stdas_k)) + (1-wask)*meanat_k;
end

Cs = zeros([1,size(has_k,2)]);
Co = zeros([1,size(has_k,2)]);
sum_cs = 0;
sum_co = 0;
for j=1:size(has_k,2)
    for i=1:j
        sum_cs = sum_cs + has_k(i);
        sum_co = sum_co + hao_k(i);
    end
    Cs(j) = sum_cs;
    Co(j) = sum_co;
    sum_cs = 0;
    sum_co = 0;
end

Iao = zeros([xs,ys]);

for n = 1:xs
    for m = 1:ys
        index_cs = round(((as(n,m) - min_as+1)/V));
        if index_cs > 124
            index_cs = 124;
        else
            index_cs = index_cs;
        end
        [~, index_co] = min(abs(Co - Cs(index_cs)));
        Iao(n,m) = index_co ;
            
         
    end
end
Iao = 0.5*(Iao + min_as);

max_bs = round(max(bs(:))); min_bs = round(min(bs(:)));
max_bt = round(max(bt(:))); min_bt = round(min(bt(:)));

Bbs = round((max_bs - min_bs)/V);
Bbt = round((max_bt - min_bt)/V);
hbs = zeros([1 Bbs]);
hbt = zeros([1 Bbt]);
vs = zeros([1 Bbs]);
vt = zeros([1 Bbt]);
Pbs = zeros([xs ys]);
Pbt = zeros([xt yt]);
%histogram
for i = 1:Bbs
    v(i) = min_bs + (i-1)*V;
    for n = 1:xs
        for m = 1:ys
            if i == round(((bs(n,m) - min_bs)/V)+1)
                Pbs(n,m) = 1;
                hbs(i) = hbs(i) + 1;
            end
        end
    end
end
for i = 1:Bbt

    for n = 1:xt
        for m = 1:yt
            if i == round(((bt(n,m) - min_bt)/V)+1)
                Pbt(n,m) = 1;
                hbt(i) = hbt(i) + 1;

            end
        end
    end
end


%Bk
syms t;
expr = log2(t);
smax_b = round(eval(subs(expr, t, Bbs/10)));
Bbk = round(Bbs*2^(k-smax_b));

                                                                                                                             
hbs_k = hbs;
hbt_k = hbt;

% % gradient
d_hbt_k = zeros([1 size(hbt_k,2)-1]);
d2_hbt_k = zeros([1 size(d_hbt_k,2)-1]);
for i=1:size(hbt_k,2)-1
    d_hbt_k(i) = hbt_k(i) - hbt_k(i+1);
end
for i=1:size(d_hbt_k,2)-1
    d2_hbt_k(i) = d_hbt_k(i) - d_hbt_k(i+1);
end
Rminb_k = zeros([1 size(hbt_k,2)]);
for i=1:size(d2_hbt_k,2)
    if (d_hbt_k(i)*d_hbt_k(i+1) < 0) && (d2_hbt_k(i) < 0)
        Rminb_k(i+1) = 1;
    end
end
count = 0;
for i=1:size(Rminb_k,2)
    if Rminb_k(i) == 1
        count = count + 1;
    end
end
j = zeros([1,count]);
z = 1;
for i=1:size(Rminb_k,2)
    if Rminb_k(i) == 1
        j(z) = i;
        z = z + 1;
    end
end
a = j(1);
b = j(2)-1;
% % ls mean std
sub_hbs_k = hbs_k(a:b);
meanbs_k = mean(sub_hbs_k);
stdbs_k = std(sub_hbs_k);
% % lt mean std
sub_hbt_k = hat_k(a:b);
meanbt_k = mean(sub_hbt_k);
stdbt_k = std(sub_hbt_k);
% % his_o
wbsk = k/smax_b;
hbo_k = zeros([1,size(hbs_k,2)]);

for i = 1:size(hbs_k,2)
    hbo_k(i) = (hbs_k(i) - wbsk*meanbs_k)*(((1-wbsk)*stdbt_k)/(wbsk*stdbs_k)) + (1-wbsk)*meanbt_k;
end
% % Cs Co
Cs = zeros([1,size(hbs_k,2)]);
Co = zeros([1,size(hbs_k,2)]);
sum_cs = 0;
sum_co = 0;
for j=1:size(hbs_k,2)
    for i=1:j
        sum_cs = sum_cs + hbs_k(i);
        sum_co = sum_co + hbo_k(i);
    end
    Cs(j) = sum_cs;
    Co(j) = sum_co;
    sum_cs = 0;
    sum_co = 0;
end

% % L Io
Ibo = zeros([xs,ys]);

for n = 1:xs
    for m = 1:ys
        index_cs = round(((bs(n,m) - min_bs+1)/V));
        if index_cs > 97
            index_cs = 97;
        else
            index_cs = index_cs;
        end
        [~, index_co] = min(abs(Co - Cs(index_cs)));
        Ibo(n,m) = index_co ;
            
         
    end
end
Ibo = 2.5*(Ibo + min_bs - 20);
lab_match = cat(3,Ilo,Iao,Ibo);
ct_output = lab2rgb(lab_match);
figure(1)
subplot(1,3,1);
imshow(img_s);title('Source')
subplot(1,3,2);
imshow(img_t);title('Target')
subplot(1,3,3);
imshow(ct_output);title('Result')

figure(2)
subplot(2,3,1);
plot(hls);title('ls')
subplot(2,3,2);
plot(has);title('as')
subplot(2,3,3);
plot(hbs);title('bs')
subplot(2,3,4);
plot(hlt);title('lt')
subplot(2,3,5);
plot(hat);title('at')
subplot(2,3,6);
plot(hat);title('bt')