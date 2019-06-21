function[Motion20]=s20mentai()
% 正二十面体の不変変換生成

GR = (1+sqrt(5))/2; % 黄金比
Motion20 = zeros(3,3,120); % 正二十面体の不変変換

S20 = [GR GR 1 0 0 1 -GR -1 0 0 -1 -GR;
       -1 1 0 -GR -GR 0 -1 0 GR GR 0 1;
       0 0 GR 1 -1 -GR 0 -GR -1 1 GR 0];
   
% scatter3(S20(1,:),S20(2,:),S20(3,:))

inv_123 = inv(emat(S20,1,2,3));

for i=1:6
    Motion20(:,:,20*i-19) = emat(S20,1,2,3) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-18) = emat(S20,1,3,4) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-17) = emat(S20,1,4,5) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-16) = emat(S20,1,5,6) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-15) = emat(S20,1,6,2) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-14) = emat(S20,5,6,8) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-13) = emat(S20,6,8,9) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-12) = emat(S20,6,9,2) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-11) = emat(S20,9,2,10) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-10) = emat(S20,10,2,3) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-9) = emat(S20,10,11,3) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-8) = emat(S20,11,4,3) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-7) = emat(S20,11,4,7) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-6) = emat(S20,4,5,7) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-5) = emat(S20,5,7,8) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-4) = emat(S20,12,7,8) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-3) = emat(S20,12,8,9) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-2) = emat(S20,12,9,10) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i-1) = emat(S20,12,10,11) * permatrix(i,3,3) * inv_123;
    Motion20(:,:,20*i) = emat(S20,12,11,7) * permatrix(i,3,3) * inv_123;
end

