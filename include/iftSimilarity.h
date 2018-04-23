/**
 * @file iftSimilarity.h
 * @brief Functions to compute Similarities or Distances between images.
 * @author Samuel Martins
 * @ingroup ImageScore
 */

#ifndef IFT_SIMILARITY_H
#define IFT_SIMILARITY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftFile.h"
#include "iftImage.h"
#include "iftMatrix.h"
#include "iftMetrics.h"
#include "iftSet.h"
#include "iftSeeds.h"
#include "iftString.h"


/**
 * @brief Metrics used between Label Images.
 * @author Samuel Martins
 * @date Mar 23, 2016
 * @ingroup ImageMetrics
 */
typedef enum {
    IFT_DICE,
    IFT_ASSD
} iftLabelImgMetric;



/**
 * @brief Computes the Similarity between two Binary Images using Dice.
 * @author Samuel Martins
 * @date Sep 1st, 2015
 * @ingroup ImageMetrics
 *
 * Computes the Dice Similarity between two Binary Images using Dice.\n
 * The binary image must have only the value 0 and another value (not necessarily 0 and 1).\n
 *
 * @param bin_source First Binary Image.
 * @param bin_target Second Binary Image. Commonly the ground truth.
 * @return The dice similarity between bin_source and bin_target.
 *
 * @note The dice similarity goes from 0 (no similarity) and 1 (perfect similarity).
 * 
 * @warning The function does not check if the input images are really Binary Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
double iftDiceSimilarity(const iftImage *bin_source, const iftImage *bin_target);


/**
 * @brief Computes the similarity between two Label Images using Dice.
 * @author Samuel martins
 * @date Feb 18, 2016
 * @ingroup ImageMetrics
 *
 * The (Binary) Dice Similarity is computed for each object (1..N).\n
 * The result is stored in the ith position of a double array, that indeed is its the object index (label value).\n
 * Index 0 holds the average for all objects.\n\n
 * 
 * E.g:\n
 * output[0] = Average Dice; \n
 * output[1] = Dice of the object 1 (its voxel value is 1)\n
 * output[N] = Dice of the object N (its voxel value is N)\n
 * 
 * @param label_img_source Source Label Image.
 * @param label_img_target Target Label Image. Commonly the ground truth.
 * @param n_objects Number of objects of the input label images.
 * @return A double array with the dices. Index 0 has the average dice and ith index has the dice of
 * the ith object.
 *
 * @note The dice similarity goes from 0 (no similarity) and 1 (perfect similarity).
 * 
 * @warning The function does not check if the input images are really Label Images with <b>n_objects</b> objects,
 * whose labels goes from 0 to n_objects.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
iftDblArray *iftDiceSimilarityMultiLabel(const iftImage *label_source,
                                         const iftImage *label_target, int n_objects);


/**
 * @brief Computes the Average Surface Distance between two Binary Images.
 * @author Samuel Martins
 * @date Feb 21, 2016
 * @ingroup ImageMetrics
 *
 * Computes the Mean Euclidean Distance between the boundary spels from the Source Binary Image to the boundary
 * spels of the Target Binary Image.\n
 * For each boundary spel from the Source Binary Image, it computes the minimum euclidean distance to a boundary spel
 * of the Target Binary Image. Then, the average minimum distance is computed and returned.
 * EDT (implemented with IFT) is used to compute the euclidean distances.\n
 * 
 * @param bin_source Source Binary Image.
 * @param bin_target Target Binary Image.
 * @return The Mean Euclidean Distance.
 * 
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).
 * @note It is not symmetric because it only computes the distance <b>one-way</b> (from the sourceto the target binary image).
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Binary Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
double iftASD(const iftImage *bin_source, const iftImage *bin_target);


/**
 * @brief Computes the Average Surface Distance between two Label Images with label from 0 to N.
 * @author Samuel Martins
 * @date Feb 21, 2016
 * @ingroup ImageMetrics
 *
 * The (Binary) ASD is computed for each object (1..N).\n
 * The result is stored in the ith position of a double array, that indeed is its the object index (label value).\n
 * Index 0 holds the average for all objects.\n\n
 * 
 * E.g:\n
 * output[0] = Average ASD; \n
 * output[1] = ASD of the object 1 (its voxel value is 1)\n
 * output[N] = ASD of the object N (its voxel value is N)\n
 *
 * @param label_source Source Label Image.
 * @param label_target Target Label Image.
 * @param n_objects Number of Objects of the two Label Images.
 * @return A double array with all ASD results. Index 0 has the average ASD and ith index has the ASD of
 * the ith object.
 *
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).\n
 * @note It is not symmetric because it only computes the distance <b>one-way</b> (from the sourceto the target binary image).
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Label Images with <b>n_objects</b> objects,
 * whose labels goes from 0 to n_objects.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
iftDblArray *iftASDMultiLabel(const iftImage *label_source, const iftImage *label_target, int n_objects);


/**
 * @brief Computes the Average Symmetric Surface Distance (ASSD) between two Binary Images.
 * @author Samuel Martins
 * @date Sep 14th, 2015
 * @ingroup ImageMetrics
 *
 * Computes the Average Distance from Source Binary Image to the Target Binary Image and vice versa (Symmetric).\n
 * The binary image must have only the value 0 and another value.
 * 
 * @param bin_source Source Binary Image.
 * @param bin_target Target Binary Image.
 * @return The resulting ASSD.
 * 
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).
 * @note The binary image must have only the value 0 and another value.
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Binary Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
double iftASSD(const iftImage *bin_source, const iftImage *bin_target);


/**
 * @brief Computes the Average Symmetric Surface Distance (ASSD) between two Label Images.
 * @author Samuel Martins
 * @date Feb 20, 2016
 * @ingroup ImageMetrics
 *
 * The (Binary) ASSD is computed for each object (1..N).\n
 * The result is stored in the ith position of a double array, that indeed is its the object index (label value).\n
 * Index 0 holds the average for all objects.\n\n
 * 
 * E.g:\n
 * output[0] = Average ASSD; \n
 * output[1] = ASSD of the object 1 (its voxel value is 1)\n
 * output[N] = ASSD of the object N (its voxel value is N)\n
 *
 * @param label_source Source Label Image.
 * @param label_target Target Label Image.
 * @param n_objects Number of Objects of the two Label Images.
 * @return A double array with the ASD. Index 0 has the average ASD and ith index has the ASD of
 * the ith object.
 * 
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Label Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
iftDblArray *iftASSDMultiLabel(const iftImage *label_source, const iftImage *label_target, int n_objects);


/**
 * // TODO: Refactor this function and create a function to Multilabel
 * @brief Returns the Sum of the Gradient of a given border points.
 * @author Samuel Martins
 * @date Mar 8, 2016 
 */
double iftBorderGradMatchingScore(iftSet *borders, iftImage *grad);

/**
 * // TODO: Refactor this function and create a function to Multilabel
 * @brief Returns the sum of the maximum gradient values on the region of influence of a given set of border points.
 *
 * The region of influence is given by the euclidean distance transform using the set of border voxels as roots.
 *
 * @date Mar 11, 2016
 * @author Thiago Vallin Spina
 *
 * @param borders The set of border voxels.
 * @param grad The original image gradient.
 * @param A The adjacency relation to be used by the EDT.
 * @param max_dist The maximum distance of voxels to the border set to be considered for computing the region of influence.
 */
double iftBorderRegionOfInfluenceGradMatchingScore(iftSet *borders, iftImage *grad, iftAdjRel *A, double max_dist);


/**
 * @brief Computes the Score Matrix (Similarity or Distance Matrix) from a set of Label Images.
 * @author Samuel Martins
 * @date Oct 23, 2015
 * @ingroup ImageMetrics
 *
 * Computes the Score Matrix (Similarity or Distance Matrix) from a set of Label Images based on a given Image Metric.
 * Each row i=p and column j=p from the resulting Score Matrix corresponds to the image img_files->files[p].
 * The resulting Score Matrix can also be a Similarity Matrix or a Distance Matrix, according to the similarity function used.
 *
 * @param metric Image Metric Function used to compute the scores between the Label Images.
 * @param img_files Array with the pathnames from the Label Images to be computed.
 * @param n_objects Number of Objects that the label images have.
 * @param extra_params Extra parameters that the Image Metric Functions may have.
 * @return The resulting Score Matrix.
 * 
 * @exception The number of objects is <= 0.
 * @exception The metric uses extra parameters and this one is NULL.
 */
iftMatrix *iftComputeScoreMatrixFromLabelImgFiles(iftLabelImgMetric metric,
                                                  const iftFileSet *img_files,
                                                  int n_objects, const iftDict *extra_params);




double          iftShannonEntropy(const iftImage *image);
double          iftJointEntropy(const iftImage *image1, const iftImage *image2);
double          iftNormalizedMutualInformation(const iftImage *image1, const iftImage *image2);
float           iftMeanDistBetweenBoundaries(iftImage *bin1, iftImage * bin2);

#ifdef __cplusplus
}
#endif

#endif


