#ifndef IFT_FILTERING_H_
#define IFT_FILTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftSort.h"
#include "iftMatrix.h"
#include "iftKernel.h"

iftImage  *iftLinearFilter(const iftImage *img, iftKernel *K);
iftImage  *iftLinearFilterInRegion(iftImage *img, iftImage *mask, iftKernel *K);
iftImage  *iftFastLinearFilter(iftImage *img, iftKernel *K, char crop);
iftImage  *iftCroppedFastLinearFilterByMatrixMult(iftImage *img, iftKernel *K);
iftFImage *iftFLinearFilter(iftFImage *img, iftKernel *K);
iftFImage  *iftFLinearFilterInRegion(iftFImage *img, iftImage *mask, iftKernel *K);
iftFImage *iftFastFLinearFilter(iftFImage *img, iftKernel *K, char crop);
iftMImage *iftMLinearFilter(iftMImage *img, iftMKernel *K);
iftMImage *iftMMLinearFilter(iftMImage *img, iftMMKernel *k_bank);
iftImage  *iftMedianFilter(iftImage *img, iftAdjRel *A);
iftMImage *iftMMedianFilter(iftMImage *img, iftAdjRel *A);
iftImage  *iftModaFilter(iftImage *img, iftAdjRel *A);
iftImage  *iftSobelGradientMagnitude(iftImage *img);
iftMatrix *iftImageToMatrix(iftImage *img, iftFastAdjRel *F, char crop);
iftMatrix *iftImageAdjToMatrix(iftImage *img, iftAdjRel *A);
iftImage  *iftMatrixToImage(iftMatrix *M, int xsize, int ysize, int zsize);
iftMatrix *iftFImageToMatrix(iftFImage *img, iftFastAdjRel *F, char crop);
iftFImage  *iftMatrixToFImage(iftMatrix *M, int xsize, int ysize, int zsize);
iftMatrix *iftMImageToMatrix(iftMImage *img, iftFastAdjRel *F);
iftMImage *iftMatrixToMImage(iftMatrix *M, int xsize, int ysize, int zsize, int nbands);
iftMatrix *iftKernelToMatrix(iftKernel *K);
iftMatrix *iftMKernelToMatrix(iftMKernel *K);
iftMatrix *iftMMKernelToMatrix(iftMMKernel *k_bank);
iftImage  *iftSmoothImage(iftImage *img, iftAdjRel *A, float sigma);
iftImage  *iftSmoothImageInRegion(iftImage *img, iftImage *mask, iftAdjRel *A, float sigma);

void iftFastBilateralFilter2DAux(const float * __restrict__ _input,
                                 float       * __restrict__ _output,
                                 int   width,
                                 int   height,
                                 int 	channels,
                                 int   s_sigma,
                                 float r_sigma);

iftImage *iftFastBilateralFilter2D(iftImage *img, int s_sigma, float r_sigma);
iftMImage *iftFastMBilateralFilter2D(iftMImage *img, int s_sigma, float r_sigma);

iftImage *iftNormalizeImage(iftImage *img, iftAdjRel *A, int Imax);
iftImage *iftAlphaPooling(iftImage *img, iftAdjRel *A, int stride, float alpha);



#ifdef __cplusplus
}
#endif

#endif
