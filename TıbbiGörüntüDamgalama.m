
alfa=0.01;
%Y1=imread('D:\image-000001.jpg');
dcmfile=dicomread('D:\proje\image-000242.dcm');
%figure(1);
%imshow(dcmfile,[]); title('Damgalanacak týbbi görüntü');
dcmImage=uint8(255*mat2gray(dcmfile));
imwrite(dcmImage, 'WM_dcmImage.png', 'png');




figure(1);
imshow(dcmImage,[]);title('Damgalanacak Týbbi Görüntü');


[LL1, HL1, LH1, HH1]= dwt2(double(dcmImage),'haar');
[LL2, HL2, LH2, HH2]= dwt2(LL1,'haar');
p=size(LL2);

[U2, S2, V2]= svd(LL2);
q=size(S2);

I_w1=imread('D:\proje\qr_kod.png');
I_w1=I_w1(:,:,1);
I1_w1=imresize(I_w1,p);

figure(2); imshow(I_w1); title('Eklenecek QR Kod');

[Ud, Sd, Vd] = svd(double(I1_w1));

Sd= alfa*Sd;
Smark= S2 + Sd;

LL2_1= U2*Smark*V2';

LL1_1= idwt2(LL2_1, HL2, LH2, HH2,'haar');
Y_1= idwt2(LL1_1, HL1, LH1, HH1, 'haar');

figure(3); imshow(uint8(Y_1)); title('Damgalanmýþ Týbbi Görüntü');


z=size(Y_1);

[LL1_wm, HL1_wm, LH1_wm, HH1_wm] = dwt2 (double(Y_1),'haar');
[LL2_wm, HL2_wm, LH2_wm, HH2_wm] = dwt2(LL1_wm, 'haar');

[U_wm, S_wm, V_wm] = svd(LL2_wm);
Swyeni= (S_wm - S2)/alfa;
WMy =  Ud*Swyeni*Vd';
%WMy_1=imresize(WMy,z);
WMy2= imresize(WMy, [300 300]);
figure(4); imshow(uint8(WMy2)); title('Çýkarýlan QR Kod');

cikti= imsubtract( uint8(Y_1), dcmImage );

figure(5);
imshow(cikti); title('Ýki resim arasýndaki fark');


J = imnoise(uint8(Y_1),'salt & pepper',0.00002);
figure(6);
imshow(J); title('Tuz karabiber filtresi');

[LL1_j, HL1_j, LH1_j, HH1_j] = dwt2 (double(J),'haar');
[LL2_j, HL2_j, LH2_j, HH2_j] = dwt2(LL1_j, 'haar');

[U_j, S_j, V_j] = svd(LL2_j);
Swj= (S_j - S2)/alfa;
WMj =  Ud*Swj*Vd';
WMt= imresize(WMj, [300 300]);
figure(7);
imshow(uint8(WMt)); title('Tuz karabiber filtresinden çýkan damga');

D= rot90(uint8(Y_1));
figure(8); imshow(uint8(D)); title('Damgalanmýþ verinin 90 derece döndürülmüþ hali');





[LL1_d, HL1_d, LH1_d, HH1_d] = dwt2 (double(D),'haar');
[LL2_d, HL2_d, LH2_d, HH2_d] = dwt2(LL1_d, 'haar');

[U_d, S_d, V_d] = svd(LL2_d);
Swd= (S_d - S2)/alfa;
WMd =  Ud*Swd*Vd';

WMd2= imresize(WMd, [300 300]);

figure(9);
imshow(uint8(WMd2)); title('Damgalanmýþ ve 90 derece döndürülmüþ veriden çýkan damga');





%Gauss filtresi 
gaussFiltresi= imnoise(uint8(Y_1),'gaussian',0, 0.00002); 
figure(10);
imshow(gaussFiltresi,[]); title('Gauss Filtresi');




[LL1_g, HL1_g, LH1_g, HH1_g] = dwt2 (gaussFiltresi,'haar');
[LL2_g, HL2_g, LH2_g, HH2_g] = dwt2(LL1_g, 'haar');

[U_g, S_g, V_g] = svd(LL2_g);
Swg= (S_g - S2)/alfa;
WMg =  Ud*Swg*Vd';

WMg2= imresize(WMg, [300 300]);
figure(11);
imshow(uint8(WMg2)); title('Gauss filtresinden çýkan damga');


%Keskinleþtirme
hUnsharpFiltresi = fspecial('unsharp',0.00002);
unsharpFiltresi = imfilter(uint8(Y_1), hUnsharpFiltresi);


figure(12);
imshow(unsharpFiltresi);
title('Keskinleþtirme Filtresi');



[LL1_k, HL1_k, LH1_k, HH1_k] = dwt2 (unsharpFiltresi,'haar');
[LL2_k, HL2_k, LH2_k, HH2_k] = dwt2(LL1_k, 'haar');

[U_k, S_k, V_k] = svd(LL2_k);
Swk= (S_k - S2)/alfa;
WMk =  Ud*Swk*Vd';

WMk2= imresize(WMk, [300 300]);
figure(13);
imshow(uint8(WMk2)); title('Keskinleþtirme filtresinden çýkan damga');


%Bulanýklaþtýrma
hBlurringFiltresi = fspecial('disk',0.02);
blurred = imfilter(uint8(Y_1),hBlurringFiltresi);


figure(14);
imshow(blurred);
title('Bulanýklaþtýrma Filtresi');

[LL1_b, HL1_b, LH1_b, HH1_b] = dwt2 (blurred,'haar');
[LL2_b, HL2_b, LH2_b, HH2_b] = dwt2(LL1_b, 'haar');

[U_b, S_b, V_b] = svd(LL2_b);
Swb= (S_b - S2)/alfa;
WMb =  Ud*Swb*Vd';

WMb2= imresize(WMb, [300 300]);
figure(15);
imshow(uint8(WMb2)); title('Bulanýklaþtýrma filtresinden çýkan damga');


[peaksnr, snr]= psnr(uint8(gaussFiltresi), uint8(dcmImage));
[peaksnr1, snr1] = psnr(uint8(WMg2), uint8(I_w1));
err = immse(uint8(gaussFiltresi), uint8(dcmImage));
err2 = immse(uint8(WMg2), uint8(I_w1));
[ssimval, ssimmap] = ssim(uint8(gaussFiltresi), uint8(dcmImage));
[ssimval2, ssimmap2] = ssim(uint8(WMg2), uint8(I_w1));



R=corrcoef(double(LL2_wm),double(LL2));



fprintf('\n Ortalama karesel hata görüntü% 0.4f \n',err);
fprintf('\n Ortalama karesel hata damga% 0.4f \n',err2);






