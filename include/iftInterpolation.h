#ifndef IFT_INTERPOLATION_H_
#define IFT_INTERPOLATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftColor.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "iftMImage.h"
#include "iftDataSet.h"
#include "iftAdjacency.h"
#include "iftPlane.h"
#include "iftMatrix.h"
#include "iftSegmentation.h"
#include "iftRepresentation.h"
#include "iftMSPS.h"

// Returns the coordinate on the border of the image for the given
// dimension that is closest to the original coordinate
static inline int iftBorderInterpolateCoord(int dim_size, int p) {
  return (p < 0) ? 0 : ((p >= dim_size)? dim_size - 1 : p);
}

// Interpolates a voxel outside the float image's domain to a border
// voxel
static inline iftVoxel iftBorderInterpolateVoxel(iftImage *img, iftVoxel v) {
  v.x = iftBorderInterpolateCoord(img->xsize, v.x);
  v.y = iftBorderInterpolateCoord(img->ysize, v.y);
  v.z = iftBorderInterpolateCoord(img->zsize, v.z);

  return v;
}

// Interpolates a voxel outside the float image's domain to a border
// voxel
static inline iftVoxel iftFBorderInterpolateVoxel(iftFImage *img, iftVoxel v) {
  v.x = iftBorderInterpolateCoord(img->xsize, v.x);
  v.y = iftBorderInterpolateCoord(img->ysize, v.y);
  v.z = iftBorderInterpolateCoord(img->zsize, v.z);

  return v;
}


// Interpolates a voxel outside the multi-band image's domain to a border
// voxel
static inline iftVoxel iftMBorderInterpolateVoxel(iftMImage *img, iftVoxel v) {
  v.x = iftBorderInterpolateCoord(img->xsize, v.x);
  v.y = iftBorderInterpolateCoord(img->ysize, v.y);
  v.z = iftBorderInterpolateCoord(img->zsize, v.z);

  return v;
}

iftPlane *iftFindBestCutPlane(iftImage *weight, iftPoint pos, int xviewsize, int yviewsize);
iftPlane *iftFindBestObjectCutPlane(iftImage *obj, iftImage *weight);

/**
 * @brief Reslice a 3D image by following a specified order among axes
 * x, y, and z, and with increments dx, dy, dz from the corresponding
 * octant, which is found from dx, dy, and dz values. The input 3D
 * image starts from octant (0,0,0) and its voxels are accessed from
 * that octant by following the order XYZ with increments
 * dx=dy=dz=1. In order to reslice it from octant
 * (xsize-1,ysize-1,zsize-1), for instance, the increments must be
 * specified as dx=dy=dz=-1. If we also wish to change the voxel
 * access order to z, y, and x, then the axis order must be ZYX.
 * @author Alexandre Falcao 
 * @date Feb, 23rd 2016
 * @ingroup Interpolation
 * @param  img: A 3D input image.
 * @param  dx : increment in {-1,1} along the axis x.
 * @param  dy : increment in {-1,1} along the axis y.
 * @param  dz : increment in {-1,1} along the axis z.
 * @param  axis_order : new order among axes from {XYZ, XZY, YXZ, YZX, ZXY, ZYX}
 * @return Image that results from the reslicing of the input image
 * with increments dx, dy, dz from the corresponding octant and with a
 * given axis order.
 */

iftImage  *iftResliceImage(iftImage *img, iftPlane *pl, int xsize, int ysize, int zsize);
iftImage *iftResliceImageSimple(const iftImage *img, int dx, int dy, int dz, iftAxisOrder axis_order);
iftImage  *iftGetSlice(iftImage *img, iftPlane *pl, int xsize, int ysize);
iftImage  *iftInterp2D(const iftImage *img, float sx, float sy);
iftImage  *iftInterp(const iftImage *img, float sx, float sy, float sz);
iftImage* iftResizeImage(const iftImage *img, int xsize, int ysize, int zsize);
iftFImage  *iftFGetSlice(iftFImage *img, iftPlane *pl, int xsize, int ysize);
iftFImage *iftFInterp2D(iftFImage *img, float sx, float sy);
iftMImage *iftMInterp2D(iftMImage *mimg, float sx, float sy);
iftMImage *iftMInterp(iftMImage *img, float sx, float sy, float sz);
iftMImage* iftResizeMImage(iftMImage* img, int xsize, int ysize, int zsize);
iftFImage *iftFInterp(iftFImage *img, float sx, float sy, float sz);
iftImage  *iftShapeBasedInterp2D(iftImage *label, float sx, float sy);
iftImage  *iftShapeBasedInterp(iftImage *label, float sx, float sy, float sz);
iftImage  *iftResliceOnPrincipalAxis(iftImage *img, iftImage *bin);

#ifdef __cplusplus
}
#endif

#endif

