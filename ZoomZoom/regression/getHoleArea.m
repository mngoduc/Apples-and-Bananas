function [area_bw perim sliver_area sliver_found] = getHoleArea(r)

axis equal;
axis([-r, r,-r,r]);
frame = getframe;
[img,~] = frame2im(frame);
BW = im2bw(img, 0.5);

%Get hole area
area_ratio = bwarea(BW)/(numel(BW));
area_extra = ((2*r)^2)-(pi*r^2);
area_bw = (area_ratio * ((r*2)^2)) - area_extra;

%Get hole perimeter

P = regionprops(BW,'perimeter');
corr_factor = (2*r)/length(BW);

perim = cat(1, P.Perimeter);
perim = sum(perim);
perim = perim*corr_factor;

%Search for slivers

BW2=imfill(BW,'holes');

noholes=bwarea(BW2); % white part without holes
withholes=bwarea(BW); % white part with holes
sliver_bw = noholes-withholes; % total area in pixels of the holes alone
sliver_area_ratio = sliver_bw/(numel(BW));
sliver_area = (sliver_area_ratio * ((r*2)^2));

threshold = 10;

if sliver_bw > threshold
    sliver_found = 1;
else
    sliver_found = 0;
end

end

