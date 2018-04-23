#ifndef IFT_IMAGEFOREST_H_
#define IFT_IMAGEFOREST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftAdjacency.h"
#include "iftGQueue.h"
#include "iftSet.h"
#include "iftLabeledSet.h"
#include "iftSeeds.h"
#include "iftBMap.h"

typedef struct ift_imageforest {
  iftImage  *pathval;
  iftImage  *label;
  iftImage  *root;
  iftImage  *pred;
  iftImage  *marker;
  iftBMap   *processed;
  iftGQueue *Q;
  iftImage  *img;
  iftAdjRel *A;
} iftImageForest;

  iftImageForest  *iftCreateImageForest(iftImage *img, iftAdjRel *A);
  void             iftResetImageForest(iftImageForest *fst);
  void             iftDestroyImageForest(iftImageForest **fst);
  void             iftCompRemovalWithoutFrontier(iftImageForest *fst, iftLabeledSet *seed);
  iftSet          *iftCompRemoval(iftImageForest *fst, iftLabeledSet *seed);
  iftSet 	  *iftTreeRemoval(iftImageForest *fst, iftSet *trees_for_removal);
  iftSet 	  *iftMarkerRemoval(iftImageForest *fst, iftSet *removal_markers);
iftImage          *iftSwitchTreeLabels(iftImageForest *fst, iftLabeledSet *seed, iftImage *label);
  void             iftRelabelBorderTreesAsBkg(iftImageForest *fst);
  void 		   iftPropagateLabelMarkerAndRootToSubtree(iftImageForest *fst, iftAdjRel *A, iftImage *new_label, iftImage *new_marker, int r);
  int 		   iftPathLength(iftImage *pred, int p);

#ifdef __cplusplus
}
#endif

#endif


