#ifndef IFT_SET_H_
#define IFT_SET_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"

/**
 * Set of integer elements. Functions as a single linked list
 */
typedef struct ift_set {
  int elem;
  struct ift_set *next;
} iftSet;


/**
 * @brief Print elements in the set
 * @param S A linked list
 */
void    iftPrintSet(iftSet *S);

/**
 * @brief Inserts an element in the begining of the linked list
 * @param S A linked list
 * @param elem The element to be inserted
 */
void    iftInsertSet(iftSet **S, int elem);

/**
 * @brief Removes and returns the first element of the linked list
 * @param S A linked list
 * @return The first element
 */
int     iftRemoveSet(iftSet **S);

void    iftRemoveSetElem(iftSet **S, int elem);
void    iftDestroySet(iftSet **S);
iftSet* iftSetUnion(iftSet *S1,iftSet *S2);
iftSet* iftSetConcat(iftSet *S1,iftSet *S2);
char    iftUnionSetElem(iftSet **S, int elem);
void    iftInvertSet(iftSet **S);

/**
 * @brief Returns the size of an iftSET.
 */
int 	iftSetSize(iftSet* S);
iftSet* iftSetCopy(iftSet* S);
int     iftSetHasElement(iftSet *S, int elem);

iftSet* iftCopyOrderedSet(iftSet *S);

/**
 * @brief Converts an iftSet to a 1D-array.
 *
 * @author Samuel Martins
 * @date September 14, 2015
 * 
 * @param S iftSet to be copied.
 * @param n_elems Number of elements to be returned by reference.
 * @return The resulting integer 1D-array.
 */
iftIntArray *iftSetToArray(iftSet *S, int *n_elems);

/**
 * @brief Converts a 1D-array to a iftSet.
 *
 * @author Samuel Martins
 * @date September 14, 2015
 * 
 * @param array The integer 1D-array to be copied.
 * @param n_elems Array size.
 * @return The resulting iftSet.
 */
iftSet *iftArrayToSet(iftIntArray *array);



#ifdef __cplusplus
}
#endif

#endif

