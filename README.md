# ECSpress

input data dimension: fluorescence(time,color,x,y)
Microscopy : 
Objective :
Sample : Bright background with black spots (penetrating vessels), some line (capillaries), ,and noise (scattered pepper like)

## Mask generation & motion correction
Step: 
 1. Denoise the image using a gaussian blur to correctly label each penetrating vessels 
 2. detect most invariable in their size (supposedly penetrating vein)
 3. calculate the center coordinate of each section of penetrating vessels
 4. 
 5. align the images while cutting the edge as 

## Single penetrating vessel analysis
Step
1. Detect penetrating vessels at c (PVs)
2. Calculate the average fluorescence of the outer boundary of a PV, save distance from PV (initial: 0)
3. Expand the outer boundary by one pixel (distance: +1) and do the same thing as step 2
4. Repeat 3 until 

## Multiple penetrating vessels
