#ifndef IFT_REPRESENTATION_H_
#define IFT_REPRESENTATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftAdjacency.h"
#include "iftGQueue.h"
#include "iftImage.h"
#include "iftRadiometric.h"
#include "iftHistogram.h"
#include "iftFImage.h"
#include "iftSeeds.h"
#include "iftSegmentation.h"

  iftFImage      *iftGeodesicDistTrans(iftSet *S, iftImage *mask, iftAdjRel *A);
  iftImage       *iftDistTrans(iftImage *bin, iftAdjRel *A, char side);
  iftImage       *iftBorderDistTrans(iftImage *label, iftAdjRel *A);
  iftImage       *iftShellDistTrans(iftImage *bin, iftAdjRel *A, char side, float max_dist);
  void            iftDistTransRootMap(iftImage *bin, iftAdjRel *A, char side, iftImage **dist, iftImage **root);
  iftFImage      *iftSignedDistTrans(iftImage *bin, iftAdjRel *A);
  iftFImage      *iftShellSignedDistTrans(iftImage *bin, iftAdjRel *A,float max_dist);
  iftFImage      *iftMSSkel(iftImage *bin); /* multiscale skeletons by geodesic length */

  /* compute surface skeleton by thresholding the geodesic skeleton
     and then select a given number of components */

  iftImage *iftSurfaceSkeleton(iftFImage *skel, float thres, int number_of_components); 

  // This function computes the multiscale skeleton and returns the distance transformed
  // used for computing it. Both maps can be used to compute a Medial Axis Transform.
  iftFImage 	 *iftMSSkel2DDistMap(iftImage *bin, iftAdjRel *A, char side, iftImage **dist);
  // This function should be used when the euclidean distance to the skeleton is irrelevant.
  // It is a wrapper for function iftMSSkel2DDistMap
  iftFImage      *iftMSSkel2D(iftImage *bin, iftAdjRel *A, char side);
  iftImage       *iftIntMSSkel2D(iftImage *bin, iftAdjRel *A, char side);
  // This function computes the border of all labels in the label image and then uses the
  // given set of border pixels to compute a distance transform from them.
  // These functions are a superset of the binary case.
  void  	  iftMultiLabelDistTransFromBorders(iftImage *label, iftAdjRel *A, char side, iftImage **dist, iftImage **root);

  /**
   * @brief Computes the distance transform for a multi-label image from a given set S.
   * S is usually composed of the boundary pixels for all labels. This function generalizes
   * iftDistTrans and iftShellDistTrans since it accepts a maximum distance parameter.
   *
   * @author Thiago Vallin Spina
   *
   * @param S An iftSet representing the object's border or any set of points from which the EDT will be computed.
   * @param label A label image used to determine the EDT's size.
   * @param A An adjacency relation used to compute the EDT.
   * @param side Determines whether the EDT will be computed inside, outside, or on both sides of the union of labels.
   * @param max_dist The maximum distance until which the EDT will be computed.
   * @param dist Returns the computed EDT.
   * @param root Returns the root value in S for each spel.
   */
  void            iftMultiLabelShellDistTransFromSet(iftSet *S, iftImage *label, iftAdjRel *A, char side, double max_dist, iftImage **dist, iftImage **root);


  void 		  iftMultiLabelDistTransFromSet(iftSet *S, iftImage *label, iftAdjRel *A,	char side, iftImage **dist, iftImage **root);

  /**
   * @brief Computes the distance transform with root, label, and path propagation.
   *
   * @author Alexandre Falcao
   *
   * @param bin:   input binary image with one or multiple objects.
   * @param A:     adjacency relation for the EDT. 
   * @param side:  side of the boundary in which the EDT will be computed --- IFT_INTERIOR, IFT_EXTERIOR or IFT_BOTH.
   * @param root:  output root map (the closest boundary point in the image). 
   * @param label: output label map (the number of the closest boundary in the image). 
   * @param dist:  output distance map (the distance to the closest boundary point).
   * @param pred:  output predecessor map (the predecessor in the optimum path from the closest boundary point).
   */
  void            iftDistTransAttributes(iftImage *bin, iftAdjRel *A, char side, iftImage **root, iftImage **label, iftImage **dist, iftImage **pred);


  /**
   * @brief Computes the number of descendants each node of the forest
   * has in a target set, as represented by a binary image.
   *
   * @author Alexandre Falcao
   *
   * @param pred:  input spanning forest.
   * @param target_bin: binary image that represents the target set.
   */

  iftImage *iftNumberOfDescendantsInTarget(iftImage *pred, iftImage *target_bin);
  /**
   * @brief Computes the number of ascendants each node has in the forest.
   *
   * @author Alexandre Falcao
   *
   * @param pred:  input spanning forest.
   */

  iftImage *iftNumberOfAscendants(iftImage *pred);

  /**
   * @brief Computes the number of descendants each node has in the forest.
   *
   * @author Alexandre Falcao
   *
   * @param pred:  input spanning forest.
   */

  iftImage *iftNumberOfDescendants(iftImage *pred);

/**
 * @brief Signed version of function iftMultiLabelDistTransFromBorders.
 *
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016 
 *
 * @param label Label image.
 * @param A Adjacency relation to establish local neighborhood.
 * @param dist Signed euclidean distance transform result.
 * @param root The root map. (May be NULL if not desired).
 *
 */
  iftFImage     *iftMultiLabelSignedDistTransFromBorders(iftImage *label, iftAdjRel *A, char side, iftImage **root);
  void           iftLabelRootPropagation(iftImage *bin, iftAdjRel *A, char side, iftImage **root, iftImage **label, iftImage **dist);
  iftImage       *iftRootPropagation(iftImage *bin, iftAdjRel *A, char side, float max_dist);
  iftFImage      *iftLabelContPixelByGeoLen(iftImage *bin);
  iftImage       *iftLiftImage(iftImage *img);
  iftImage       *iftDropImage(iftImage *bin);
  iftFImage      *iftIntegralImage(iftImage *img);
  float           iftGetIntegralValueInRegion(iftFImage *integ, iftVoxel *v, int npts);
  iftImage       *iftMarkGeometricCenters(iftImage *bin);
  iftImage       *iftComponentSizes(iftImage *bin, iftAdjRel *A);
  iftImage       *iftBorderSize(iftImage *bin);

/**
 * @brief Runs the medial axis transform in 2D, such that the output
 * image is a 2D skeleton image with the squared Euclidean distance
 * from the skeleton points to the objects' boundaries.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  binary image with one or multiple objects. 
 * @param  scale threshold for binarization of the multiscale skeleton.
 * @param  object side in which the skeleton will be computed: IFT_INTERIOR, IFT_EXTERIOR, or IFT_BOTH.
 * @return The skeleton image whose values are the squared EDT values for shape reconstruction.
 * 
 */

  iftImage       *iftMedialAxisTrans2D(iftImage *bin, float scale_thres, char side);

/**
 * @brief Reconstructs the shape from the medial axis image. It is a
 * shape smoothing filter.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  medial axis image
 * @param  reconstruction value
 * @return Binary image with the reconstructed shapes. 
 * 
 */

  iftImage       *iftShapeReconstruction(iftImage *medial_axis, int value);


/**
 * @brief Terminal points of a binary skeleton.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  Binary skeleton.
 * @return Binary image with the terminal points of the input skeleton.
 * 
 */

  iftImage       *iftTerminalPoints2D(iftImage *skel);

/**
 * @brief  Extract the set of terminal points of a binary skeleton.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  Binary skeleton.
 * @return Set of terminal points of the input skeleton.
 * 
 */
  
  iftSet         *iftTerminalPointSet2D(iftImage *skel);

  
#ifdef __cplusplus
}
#endif

#endif
