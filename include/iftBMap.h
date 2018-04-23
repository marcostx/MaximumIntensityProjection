#ifndef IFT_BMAP_H_
#define IFT_BMAP_H_

#ifdef __cplusplus
extern "C" {
#endif


#include "iftCommon.h"

#define   iftBMapValue(b,n) ((b->val[n>>3]&(1<<(n&0x07)))!=0)
#define   iftBMapSet0(b,n) b->val[n>>3]&=((~0)^(1<<(n&0x07)));
#define   iftBMapSet1(b,n) b->val[n>>3]|=(1<<(n&0x07));
#define   iftBMapToggle(b,n) b->val[n>>3]^=(1<<(n&0x07));

typedef struct ift_bitmap {
  char *val;
  int   n, nbytes;
} iftBMap;

iftBMap  *iftCreateBMap(int n);
void      iftDestroyBMap(iftBMap **b);
void      iftFillBMap(iftBMap *b, int value);
iftBMap  *iftCopyBMap(iftBMap *src);


/**
 * @brief Read a Bit Map from File.
 * @author Samuel Martins
 * @date Apr 12, 2016
 * @ingroup BitMap
 *
 * A string is read from the file and is assigned at the char* in the bit map structure. \n
 * Thus, we have one byte more allocated due to the '\0', but the number of bytes from the structure
 * ignores it.
 */
iftBMap *iftReadBMap(const char *path);


/**
 * @brief Write a Bit Map.
 * @author Samuel Martins
 * @date Apr 12, 2016
 * @ingroup BitMap
 * 
 * In order to avoid to write multiple chars in file, we build a string and write it. \n
 * For that, we need to allocate one byte more for the '\0', but nothing changes because of this.
 */
void iftWriteBMap(const iftBMap *bmap, const char *path);

#ifdef __cplusplus
}
#endif

#endif
