#ifndef IFT_DHEAP_H_
#define IFT_DHEAP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"

#define MINVALUE   0 /* define heap to remove node with minimum value */
#define MAXVALUE   1 /* define heap to remove node with maximum value */

#define iftDad(i) ((i - 1) / 2)
#define iftLeftSon(i) (2 * i + 1)
#define iftRightSon(i) (2 * i + 2)
#define iftSetRemovalPolicyDHeap(a,b) a->removal_policy = b 

typedef struct ift_dheap {
  double *value;
  char  *color;
  int   *node;
  int   *pos;
  int    last;
  int    n;
  char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
} iftDHeap;

iftDHeap *iftCreateDHeap(int n, double *value);
void      iftDestroyDHeap(iftDHeap **H);
char      iftFullDHeap(iftDHeap *H);
char      iftEmptyDHeap(iftDHeap *H);
char      iftInsertDHeap(iftDHeap *H, int pixel);
int       iftRemoveDHeap(iftDHeap *H);
void      iftRemoveDHeapElem(iftDHeap *H, int pixel);
void      iftGoUpDHeap(iftDHeap *H, int i);
void      iftGoDownDHeap(iftDHeap *H, int i);
void      iftResetDHeap(iftDHeap *H);

#ifdef __cplusplus
}
#endif

#endif



