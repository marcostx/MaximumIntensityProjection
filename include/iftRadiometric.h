#ifndef IFT_RADIOMETRIC_H_
#define IFT_RADIOMETRIC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "iftColor.h"
#include "iftGQueue.h"
#include "iftHistogram.h"

/* Radiometric Transformations */

iftImage *iftLinearStretch(iftImage *img, double f1, double f2, double g1, double g2);
iftFImage *iftFLinearStretch(iftFImage *img, double f1, double f2, double g1, double g2);
iftImage *iftNormalize(iftImage *img, double minval, double maxval);


/**
 * @brief Normalize an Image to the range [minval, maxval] ignoring voxel outliers with high brightness.
 *
 * It computes the normalized accumulated histogram to figure out the percentage of voxels from 0 until 
 * each voxel v. When attaining a given percentage perc, it will ignore the brightness of the 
 * remaining voxels (outliers)
 * 
 * @param  img    Input Image
 * @param  minval Inferior bound of the normalization.
 * @param  maxval Inferior bound of the normalization.
 * @param  perc   Percentage of true voxels to be considered in normalization (Suggestion: 0.98)
 * @return        Normalized Image.
 *
 * @author Samuka Martins
 * @date Jun 20, 2017
 */
iftImage *iftNormalizeWithNoOutliers(const iftImage *img, int minval, int maxval, float perc);
iftImage *iftWindowAndLevel(iftImage *img, int width, int level, int maxval);
iftImage *iftGaussianStretch(iftImage *img, double mean, double stdev, int maxval);
iftImage *iftExponenStretch(iftImage *img, double f1, double f2, double g1, double g2);

/**
 * @brief Equalizes (uniformly) the histogram of an image by the sorting method.
 *
 * Input image will be equalized uniformly in the range of [0, max_value].
 * 
 * 1) Sort all pixels from <img> in ascending order of brightness; \n
 * 2) Split the sorted pixels in (max_val + 1) buckets, k = 0, 1, ..., max_val, with the same number of voxels; \n
 * 3) Equalized Image <out>, such that out[p] = the bucket k where the voxel p has its brightness img[p] mapped to. \n\n
 *
 * POSSIBLE PROBLEM:
 * If the passed max_val <= 0, it will be considered the maximum value of the image.
 * Noises with high brightness values can result in a wrong equalization, stretching the brightness range incorrectly.
 * 
 * @param  img     Image to be equalized.
 * @param  max_val Maximum value of the output image. If <= 0, it is considered the maximum value of the image. 
 * @return         Equalized Image.
 *
 * @author Samuka
 * @date Jan 3, 2017
 */
iftImage *iftEqualize(const iftImage *img, int max_val);


iftImage *iftMatchHistogram(iftImage *img, iftImage *ref);


/**
 * @brief Matches the histogram of an image <img> with a reference image <ref> according to
 * the Traditional Algorithm.
 *
 * This traditional algorithm works as follows:
 * 1) builds the normalized accumulated histograms for the images (himg and href)
 * 2) For each pixel p in img
 *    tries to find the bucket of href with value exactly equal to himg[img->val[p]]
 *    this will be value of p in the matched output image. 
 *
 * OBS: Equal images (img == ref) does not necessarily return an exact copy of ref.
 * A set of voxels can have values slightly different from the original ones, because it gets
 * the first ref img's hist bucket with the same value of the img's hist.
 * Because it is an accumulated histogram, different sequential buckets can have the same value.
 * 
 * @param  img Image whose histogram will be matched.
 * @param  ref Reference image.
 * @return     Resulting Image with Matched Histogram.
 *
 * @author Samuka
 * @date Jan 3, 2017
 */
iftImage *iftTradMatchHistogram(const iftImage *img, const iftImage *ref);

  /**
 * @brief Selects the highest voxel intensity among a percentage of voxels in the image. 
 *
 * @param  img        Input image 
 * @param  voxel_perc Percentage of voxels in (0.0,1.0)
 * @return maximum value among the given percentage of voxels
 *
 * @author Samuka and Falcao
 * @date Jul 4th, 2017
 */
  int iftMaximumValueAmongPercOfVoxels(const iftImage *img, const float voxel_perc);


/* Radiometric Resolution */
int 	  iftRadiometricResolution(iftImage *img);


#ifdef __cplusplus
}
#endif

#endif
