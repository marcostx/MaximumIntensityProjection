#ifndef IFT_LABELEDSET_H_
#define IFT_LABELEDSET_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftSet.h"

typedef struct ift_labeledset {
  int elem;
  int label;
  int marker;
  int handicap;
  struct ift_labeledset *next;
} iftLabeledSet;

void iftInsertLabeledSet(iftLabeledSet **S, int elem, int label);
void iftInsertLabeledSetMarkerAndHandicap(iftLabeledSet **S, int elem, int label, int marker, int handicap);
int  iftRemoveLabeledSet(iftLabeledSet **S, int *label);
void iftRemoveLabeledSetElem(iftLabeledSet **S, int elem);
void iftDestroyLabeledSet(iftLabeledSet **S);

/**
 * @brief Inserts a Set with label <label> into a Labeled Set. All elements are removed from Set and it
 * is assigned to NULL.
 * 
 * @param S     Set to be inserted.
 * @param label Label of the elements from set.
 * @param T     Reference of the Labeled Set
 *
 * @author samuka
 * @date Jan 2, 2017
 */
void iftInsertSetIntoLabeledSet(iftSet **S, int label, iftLabeledSet **T);


/**
 * @brief Translates a Labeled Set by a displacement vector.
 * @param  S        Labeled Set to be translated.
 * @param  img      Image where the labeled set is.
 * @param  disp_vec Displacement vector to translate the labeled set.
 * @return          Translated Labeled Set.
 */
iftLabeledSet *iftTranslateLabeledSet(const iftLabeledSet *S, const iftImage *img, iftVector disp_vec);




iftLabeledSet* iftCopyLabeledSet(iftLabeledSet *s);
iftLabeledSet* iftCopyOrderedLabeledSet(iftLabeledSet *s);
  void iftPrintLabeledSet(iftLabeledSet *S);

/**
 * @brief Copies the labeled set in reverse order
 *
 * @author Thiago Vallin Spina
 *
 * @param s The input seed set.
 * @return The reversed seed set.
 */
iftLabeledSet* iftReverseLabeledSet(iftLabeledSet *s);

void iftConcatLabeledSet(iftLabeledSet **S1,iftLabeledSet **S2);
void iftRemoveSubsetLabeledSet(iftLabeledSet **S1,iftLabeledSet **S2);

int iftLabeledSetSize(iftLabeledSet *s);

iftSet* iftLabeledSetToSet(iftLabeledSet *S, int lb);


/**
 * @brief Copies the elements in the labeled set ot and iftSet* regardless of their labels.
 *
 * @author Thiago Vallin Spina
 *
 * @param S labeled set
 * @return An iftSet* with all elements in <S>
 */
iftSet* iftLabeledSetElemsToSet(iftLabeledSet *S);
iftLabeledSet* iftCopyLabels(iftLabeledSet *S, int lb);

iftSet* iftLabeledSetMarkersToSet(iftLabeledSet *S, int marker);

/**
 * @brief Copies the set of markers into another labeled set.
 *
 * @author Thiago Vallin Spina
 * @date Jan 25, 2016
 *
 * @param S The original set of seeds.
 * @param marker The marker ID that should be copied.
 *
 * @return The filtered labeled set.
 */
iftLabeledSet* iftCopyLabeledSetMarkers(iftLabeledSet *S, int marker);

/**
 * @brief Verifies if an element belongs to the labeled seed set.
 *
 * @param S The labeled seed set.
 * @param elem The desired element.
 * @return 1 if found, 0 otherwise
 * @sa iftSetHasElement
 */
int iftLabeledSetHasElement(iftLabeledSet *S, int elem);

#ifdef __cplusplus
}
#endif

#endif

