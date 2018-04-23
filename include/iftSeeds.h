#ifndef IFT_SEEDS_H_
#define IFT_SEEDS_H_

#ifdef __cplusplus
extern "C" {
#endif


#include "iftCommon.h"
#include "iftImage.h"
#include "iftAdjacency.h"
#include "iftSet.h"
#include "iftLabeledSet.h"
#include "iftGQueue.h"
#include "iftLIFO.h"
#include "iftFIFO.h"
#include "iftBMap.h"
#include "iftMathMorph.h"
#include "iftDataSet.h"
#include "iftRegion.h"
#include "iftMImage.h"

iftLabeledSet *iftReadSeeds(char *filename, iftImage *img);
iftLabeledSet *iftReadSeedsComplete(const char *filename, const iftImage *img);
iftLabeledSet *iftReadSeeds2D(const char *filename, const iftImage *img);
iftSet 	      *iftExtractRemovalMarkers(iftLabeledSet **s);
iftLabeledSet *iftLabelObjBorderSet(iftImage *bin, iftAdjRel *A);
iftLabeledSet *iftImageBorderLabeledSet(iftImage *img);
iftLabeledSet *iftLabelCompSet(iftImage *bin, iftAdjRel *A);
iftLabeledSet *iftFuzzyModelToLabeledSet(iftImage *model);
iftLabeledSet *iftMAdjustSeedCoordinates(iftLabeledSet *S, iftMImage *input, iftMImage *output);
iftLabeledSet *iftAdjustSeedCoordinates(iftLabeledSet *Sin, iftImage *orig, iftMImage *output);

iftSet        *iftImageBorderSet(iftImage *img);
iftLabeledSet *iftMultiObjectBorderLabeledSet(iftImage *img, iftAdjRel *A);
iftImage *iftLabeledBorders(iftImage *label, iftAdjRel *A);

/**
 * @brief Gets the voxel's index from all object voxels (val > 0), returning them into a set
 * @param  label_img Label Image.
 * @return           Set with the voxel's indices from all objects.
 */
iftSet *iftObjectToSet(const iftImage *label_img);


/**
 * @brief Finds and returns an iftSet with the coordinates of all border spels of a label image for a given adjacency <A>.
 * @param label Image label.
 * @param A Adjacency Relation.
 * @return The iftSet with the coordinates of all border spels of the label image <label>.
 */
iftSet *iftObjectBorderSet(const iftImage *label, iftAdjRel *Ain);

/**
 * @brief Finds and returns an iftBMap with the coordinates of all border spels of a label image for a given adjacency <A>.
 * @param label Image label.
 * @param A Adjacency Relation.
 * @param n_border_spels Number of border pixels/voxels of the label image that will be returned by reference.
 * @return The iftBMap with the coordinates of all border spels of the label image <label>.
 */
iftBMap *iftObjectBorderBitMap(const iftImage *label, const iftAdjRel *A, int *n_border_spels);


/**
 * @brief Stores all pixels from a Bin Image (except the background 0) to an input Set.
 * @author Samuel Martins
 * @dataset Jul 12, 2016
 * 
 * @param  img Input Image.
 * @return     The resulting Set.
 *
 * @note The input Set S may already have elements.
 */
void iftBinImageToSet(const iftImage *img, iftSet **S);

/**
 * @brief   Copy all labeled region voxels to a Labeled Set.
 * @author  Alexandre Falcao
 * @dataset Feb 07th, 2017
 * 
 * @param  label Labeled region image.
 * @return The resulting labeled set.
 *
 * @note The input labeled set S may already have elements.
 */
  
void iftRegionsToLabeledSet(const iftImage *label, iftLabeledSet **S);
  


iftSet        *iftEndPoints(iftImage *skel, iftAdjRel *A);
iftSet        *iftFindPathOnSkeleton(iftImage *skel, iftAdjRel *A, int src, int dst);
iftSet        *iftSkeletonPoints(iftImage *skel);
iftImage      *iftObjectBorders(iftImage *label, iftAdjRel *A);
iftImage      *iftLabelObjBorders(iftImage *bin, iftAdjRel *A);
iftImage      *iftEasyLabelComp(iftImage *bin, iftAdjRel *A);
iftImage      *iftLabelComp(iftImage *bin, iftAdjRel *A);
iftImage      *iftSelectLargestComp(iftImage *bin, iftAdjRel *A);
iftImage      *iftSelectSmallestComp(iftImage *bin, iftAdjRel *A);
iftImage      *iftSelectKLargestComp(iftImage *bin, iftAdjRel *A, int K);
iftImage      *iftSelectKSmallestComp(iftImage *bin, iftAdjRel *A, int K);
iftImage      *iftComponentArea(iftImage *bin, iftAdjRel *A);
iftImage      *iftSelectCompAboveArea(iftImage *bin, iftAdjRel *A, int thres);
iftImage      *iftSelectCompBelowArea(iftImage *bin, iftAdjRel *A, int thres);
iftImage      *iftSelectCompInAreaInterval(iftImage *bin, iftAdjRel *A, int thres_min, int thres_max);
iftImage      *iftRegionalMaxima(iftImage *img);
iftImage      *iftRegionalMinima(iftImage *img);
iftImage      *iftRegionalMaximaInRegion(iftImage *img, iftImage *mask);
iftImage      *iftRegionalMinimaInRegion(iftImage *img, iftImage *mask);
iftImage      *iftRootVoxels(iftImage *pred);
iftImage      *iftLeafVoxels(iftImage *pred, iftAdjRel *A);


  /* Assume that label contains regions with values in
     {1,2,..,n}. That is, label equal to zero is considered background
     and it does not define any region */

iftImage      *iftRegionArea(iftImage *label);
iftImage      *iftSelectLargestRegion(iftImage *label);
iftImage      *iftSelectSmallestRegion(iftImage *label);
iftImage      *iftSelectKLargestRegions(iftImage *label, int K);
iftImage      *iftSelectKLargestRegionsAndPropagateTheirLabels(iftImage *label, iftAdjRel *A, int K);
iftImage      *iftSelectKSmallestRegions(iftImage *label, int K);
iftImage      *iftSelectRegionsAboveArea(iftImage *label, int thres);
iftImage      *iftSelectRegionsAboveAreaAndPropagateTheirLabels(iftImage *label, int thres);
iftImage      *iftSelectRegionsBelowArea(iftImage *label, int thres);
/**
 * @brief Selects the regions with area size within the given interval.
 *
 * The regions are relabeled automatically in the output image.
 *
 * @param in_label The input label.
 * @param min_thres The minimum area size.
 * @param max_thres The maximum area size.
 * @return The relabeled image with removed regions.
 *
 */

iftImage      *iftSelectRegionsInAreaInterval(iftImage *label, int min_area, int max_area);
/**
 * @brief Selects the regions with area size within the given interval.
 *
 * The regions are relabeled automatically in-place, as opposed to iftSelectRegionsInAreaInterval.
 *
 * @param in_label The input label.
 * @param min_thres The minimum area size.
 * @param max_thres The maximum area size.
 * @sa iftSelectRegionsInAreaInterval
 *
 */
void           iftSelectRegionsInAreaIntervalInplace(iftImage *label, int min_area, int max_area);

  /* ---------------------------- */

iftImage      *iftLabelContPixel(iftImage *bin);
char           iftValidContPoint(iftImage *bin, iftAdjRel *A, iftAdjRel *L, iftAdjRel *R, int p);

void iftWriteSeeds2D(char* filename, iftLabeledSet* seed, iftImage* image);
void iftWriteSeeds(const char *filename, iftLabeledSet* seed, const iftImage* image);
void iftWriteSeedsComplete(const char *filename, iftLabeledSet *seed, const iftImage *image, ...);

iftImage	  *iftSeedImageFromLabeledSet(iftLabeledSet* labeled_set, iftImage *image);

iftLabeledSet *iftLabeledSetFromSeedImage(iftImage* seed_image);
iftLabeledSet *iftLabeledSetFromSeedImageMarkersAndHandicap(iftImage* seed_image, iftImage *marker, iftImage *handicap);

/**
 * @brief This function computes the segmentation error components and relabels them to ensure that each one is given a
 * unique id (refactored from iftBorderMarkersForPixelSegmentation)
 *
 * @author Thiago Vallin Spina
 *
 * @param gt_image Ground truth image
 * @param label Segmentation result. It may be NULL, in which case we return an error component image with one component
 * per label
 * @param adj_relabeling Adjacency relation considered for relabeling the error components
 */
iftImage *iftRelabelSegmentationErrorComponents(iftImage *gt_image, iftImage *label, iftAdjRel *adj_relabeling);
iftLabeledSet *iftBorderMarkersForPixelSegmentation(iftImage *grad_image, iftImage *gt_image, float border_distance);
iftLabeledSet *iftGeodesicMarkersForSegmentation(iftImage *gt_image, iftImage *label);
iftLabeledSet* iftBorderMarkersForSuperpixelSegmentation(iftImage* label_image,iftImage* gt_image, iftDataSet* dataset);

//Pops the first "nelem" elements from "lset" with label "label"
iftLabeledSet* iftGetSeeds(iftLabeledSet* S, int nelem, int label);
iftLabeledSet* iftGetMisclassifiedSeeds(iftLabeledSet* S, int nelem, int label, iftImage* gt_image, iftImage* cl_image);

// Binary Segmentation only
int iftCheckNewSeeds(int *nelem, int length);

int iftSelectCircularRobotSeeds(iftImage *seed_image, iftBMap *used_seeds, iftImage *gt, double dist_border,
                                double max_marker_radius, double min_marker_radius, iftAdjRel *distance_border,
                                int center_seed);
//int iftMarkersFromMisclassifiedSeeds(iftImage* seed_image, iftLabeledSet* all_seeds, iftBMap* used_seeds, int nseeds,iftImage* gt_image, iftImage* cl_image, int dist_border, int max_marker_radius, int min_marker_radius);
int iftMarkersFromMisclassifiedSeeds(iftImage *seed_image, iftLabeledSet *all_seeds, iftBMap *used_seeds, int nseeds,
                                     int number_of_labels, iftImage *gt, iftImage *label, int dist_border,
                                     int max_marker_radius, int min_marker_radius);

void iftWriteSeedsOnImage(iftImage* image, iftLabeledSet* seed);

int iftRootVoxel(iftImage *pred, int p);

iftImage  *iftFastLabelComp(iftImage *bin, iftAdjRel *A);

iftSet *iftBinaryMaskToSet(iftImage *mask);

#ifdef __cplusplus
} 
#endif

#endif
