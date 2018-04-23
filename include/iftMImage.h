#ifndef IFT_MIMAGE_H_
#define IFT_MIMAGE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "iftColor.h"
#include "iftAdjacency.h"


#define iftMGetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define iftMGetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define iftMGetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))
#define iftMGetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])
#define iftMDiagonalSize(s) (iftRound(sqrtf(s->xsize*s->xsize + s->ysize*s->ysize + s->zsize*s->zsize)))

/* Multiband image: iftMImage *img; img->band[b].val[p]  */

typedef struct ift_mimage {
  iftBand *band; /* a list of images, one for each band */
  int      xsize,ysize,zsize; /* image dimensions of each band */
  float    dx,dy,dz;  /* voxel size */
  int     *tby, *tbz; /* LUT to speed up index to/from coordinate conversions */
  int      n,m; /* number of voxels and number of bands */
} iftMImage;

  int         iftMXSize(iftMImage *img);
  int         iftMYSize(iftMImage *img);
  int         iftMZSize(iftMImage *img);
  iftVoxel    iftMGetVoxelCoord(iftMImage *img, int p);
  iftMImage  *iftCreateMImage(int xsize,int ysize,int zsize, int nbands);
  void        iftDestroyMImage(iftMImage **img);
  iftMImage  *iftCopyMImage(iftMImage *img);
  char        iftMValidVoxel(iftMImage *img, iftVoxel v);
  void        iftMCopyVoxelSize(iftMImage *img1, iftMImage *img2);
/**
 * @brief Copies the voxel size from an iftImage to an iftMImage
 *
 * @author Thiago Vallin Spina
 *
 * @param img1 An iftImage.
 * @param img2 An iftMImage.
 */
  void        iftMCopyVoxelSizeFromImage(iftImage *img1, iftMImage *img2);
/**
 * @brief Copies the voxel size from an iftMImage to an iftImage
 *
 * @author Thiago Vallin Spina
 *
 * @param img1 An iftMImage.
 * @param img2 An iftImage.
 */
  void        iftMCopyVoxelSizeToImage(iftMImage *img1, iftImage *img2);
  iftMImage  *iftImageToMImage(const iftImage *img, char color_space); /* See options for color_space in iftColor.h */
  iftImage   *iftMImageToImage(iftMImage *img, int Imax, int band);
  char        iftIs3DMImage(iftMImage *img);
  
  iftMImage   *iftReadMImage(const char *filename);
  void  	     iftWriteMImage(iftMImage *mimg, const char *filename);
  void        iftWriteMImageBands(iftMImage *mimg, char *base_filename);
  
  iftMImage  *iftMAddFrame(iftMImage *img, int bx, int by, int bz, float value);
  iftMImage  *iftMRemFrame(iftMImage *fimg, int bx, int by, int bz);
  void        iftSetMImage(iftMImage *img, float value);
  iftImage   *iftEuclMImageBasins(iftMImage *img, iftAdjRel *A);
  iftImage   *iftMImageBasins(iftMImage *img, iftAdjRel *A);
  iftImage   *iftBorderProbImage(iftMImage *img);
  iftImage   *iftRegionProbImage(iftMImage *img);
  iftImage   *iftUniformProbImage(iftMImage *img);
  void        iftMultMImageByScalar(iftMImage *Z, float scalar);
  
  iftMImage  *iftGradientVector(iftImage *img, iftImage *mask, iftAdjRel *A);
  
  /* Voxel sampling methods that return a binary image with 0/1 value,
     by using the border information to avoid outliers. */

  iftImage *iftSelectNonBorderVoxels(iftMImage *img, iftImage *mask, int nsamples);

/**
 * @brief Does a grid sampling in the multi-image
 * @param img A multi-image
 * @param mask The mask that defines the region in img to be sampled
 * @param nsamples Desired number of samples
 * @return A mask image with the selected samples
 *
 */
  iftImage *iftGridSampling(iftMImage *img, iftImage *mask, int nsamples);

  iftImage *iftAltMixedSampling(iftMImage *img, iftImage *mask, int nsamples);

  iftImage *iftSelectNSamplesFromMask(iftMImage *img, iftImage *mask1, int nsamples);

  // If band is negative then the function will search for the maximum value among all bands
  float iftMMaximumValue(iftMImage *img, int band);
  // If band is negative then the function will search for the minimum value among all bands
  float iftMMinimumValue(iftMImage *img, int band);

  /**
  * @brief Computes an iftImageTiles record given a multidimensional image and a number of tiles for each image dimension.
  *
  * The dimensions of the image tiles are equally sized to approximately divide the image into regular patches.
  * In practice, this function uses iftComputeBoundingBoxImageTiles with the entire image as a bounding box.
  *
  * @author Thiago V. Spina
  * @date Feb 26, 2015
  *
  * @param mimg    A multidimensional image that must be tiled..
  * @param ntiles_x The number of tiles in the X axis.
  * @param ntiles_y The number of tiles in the Y axis.
  * @param ntiles_z The number of tiles in the Z axis.
  * @return The allocated record.
  *
  * @warning The number of tiles may actually be less than ntiles_x*ntiles_y*ntiles_z since we prioritize the fact that
  * the patches must have approximately the same size. REFER TO tiles->ntiles_x/y/z when accessing the actual number of
  * computed tiles.
  *
  * @sa iftBoundingBoxImageTilesByEquallyDividingAxes
  */
  iftImageTiles *iftMImageTilesByEquallyDividingAxes(const iftMImage *mimg, int ntiles_x, int ntiles_y, int ntiles_z);

  /**
  * @brief Extends multidimensional image by adding the voxel coordinates. They can be added with (normalize_voxel_coord=1) or without (normalize_voxel_coord=0) normalization. 
  *
  *
  * @author Alexandre Falcao
  * @date Nov 10, 2016
  *
  * @param img    An input multidimensional image to be extended with the voxel coordinates.
  * @param normalize A flag that indicates whether or not the voxel coordinates must be normalized for addition.
  * @return the extended multidimensional image. 
  *
  */

  iftMImage *iftExtendMImageByVoxelCoord(iftMImage *img, char normalize);

/**
 * @brief Extends a MImage with only one band with its Local Binary Pattern (LPB)
 */
iftMImage *iftExtendMImageByLBP(iftMImage *img, iftAdjRel *A, char normalize);

  /**
  * @brief Extends a multidmensional image by the grayscale values of the bands of adjacent pixels and its voxel coordenates.They can be added with (normalize_voxel_coord=1) or without (normalize_voxel_coord=0) normalization.
  *
  * @author Adan Echemendia
  * @date Nov 20, 2016
  *
  * @param img The input multi-band image.
  * @param A The adjacency relation corresponding to the patch around each voxel.
  * @param normalize_voxel_coord A flag that indicates whether or not the voxel coordinates must be normalized for addition.
  * @return A extended multidimensional image.
  *
  */
  iftMImage *iftExtendMImageByAdjacencyAndVoxelCoord(iftMImage *mimg, iftAdjRel *A, char normalize_voxel_coord);

iftImage *iftMImageTilesToLabelImage(iftMImage *mimg, int n_tiles);
#ifdef __cplusplus
}
#endif

#endif


