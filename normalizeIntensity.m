% normalize the Intensity of different exp. batch
function norIMG = normalizeIntensity(img, bitwiseMask, BWedge)
% img's Attribute is double.

img0 = img.*bitwiseMask;
valueEdge = img0*BWedge;
sumEdge = sum(valueEdge(:));
averageEdge = sumEdge/sum(bitwiseMask(:));
norIMG = img0/averageEdge;