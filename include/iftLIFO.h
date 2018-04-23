#ifndef IFT_LIFO_H_
#define IFT_LIFO_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"

typedef struct ift_lifo {
  int  *LIFO;
  char *color;
  int   last;
  int   n;
} iftLIFO;

iftLIFO *iftCreateLIFO(int n);
void     iftDestroyLIFO(iftLIFO **L);
char     iftInsertLIFO(iftLIFO *L, int node);
int      iftRemoveLIFO(iftLIFO *L);
char     iftFullLIFO(iftLIFO *L);
char     iftEmptyLIFO(iftLIFO *L);
void     iftResetLIFO(iftLIFO *L);

#ifdef __cplusplus
}
#endif

#endif
