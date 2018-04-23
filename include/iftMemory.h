/**
 * @brief Memory Management functions.
 * @author Peixinho
 * @date Oct 8, 2015
 * @ingroup Memory
 */

#ifndef IFT_IFTMEMORY_H
#define IFT_IFTMEMORY_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#if !defined(__APPLE__)
    #include <malloc.h>
#endif


#include "iftBasicDataTypes.h"

/**
 * @brief Array of uchar values.
 * @author Thiago Vallin Spina
 * @date Mar 13, 2016
 * @ingroup Memory
 */
typedef struct ift_uchar_array {
    /** Number of elements. */
    size_t n;
    /** Array of uchar values. */
    uchar *val;
} iftUCharArray;


/**
 * @brief Array of integer values.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
typedef struct ift_int_array {
    /** Number of elements. */
    size_t n;
    /** Array of integer values. */
    int *val;
} iftIntArray;


/**
 * @brief Array of unsigned long values.
 * @author Thiago Vallin spina
 * @date Mar 1, 2016
 * @ingroup Memory
 */
typedef struct ift_ulong_array {
    /** Number of elements. */
    size_t n;
    /** Array of integer values. */
    ulong *val;
} iftULongArray;


/**
 * @brief Array of long values.
 * @author Thiago Vallin spina
 * @date Mar 1, 2016
 * @ingroup Memory
 */
typedef struct ift_long_array {
    /** Number of elements. */
    size_t n;
    /** Array of integer values. */
    long *val;
} iftLongArray;


/**
 * @brief Array of float values.
 * @author Thiago Vallin Spina
 * @date Mar 4, 2016
 * @ingroup Memory
 */
typedef struct ift_flt_array {
    /** Number of elements. */
    size_t n;
    /** Array of float values. */
    float *val;
} iftFloatArray;


/**
 * @brief Array of double values.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
typedef struct ift_dbl_array {
    /** Number of elements. */
    size_t n;
    /** Array of double values. */
    double *val;
} iftDblArray;


/**
 * @brief Array of strings.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
typedef struct ift_str_array {
    /** Number of elements. */
    size_t n;
    /** Array of strings. */
    char **val;
} iftStrArray;

/**
 * @brief Invoke garbage collector.
 * @author Peixinho
 * @ingroup Memory
 * @date Nov, 2016
 */
void iftGC();

/**
 * @brief The amount of memory currently being used by this process, in bytes.
 * @author Peixinho
 * @ingroup Memory
 */
size_t iftMemoryUsed();
size_t iftGetPhysicalSystemMemory();

/**
 * @brief Shows the current number of allocated objects, for debugging purpose. See also iftAllocObjectsCount().
 * @author Peixinho
 * @ingroup Memory
 */
void iftPrintAllocObjectsCount();

/**
 * @brief Prints a memory size, with an appropriate metric (KB, MB, GB).
 * @param mem Memory size in bytes
 * @author Peixinho
 * @date Jun, 2017
 */
void iftPrintMemoryFormatted(size_t mem);

/**
 * @brief Count the current number of allocated objects, for debugging purpose.
 * @author Peixinho
 * @ingroup Memory
 */
size_t iftAllocObjectsCount();

/**
 * @brief Allocate n memory blocks of size sz. See also iftFree().
 * @ingroup Memory
 * @author Peixinho
 */
void* iftAlloc(size_t n, size_t sz);

/**
 * @brief Releases a memory block previously allocated. See also iftAlloc().
 * @ingroup Memory
 * @author Peixinho
 */
void iftFree(void* mem);

/**
 * @brief Allocate reusing memory n blocks of size sz. See also iftFree().
 * @ingroup Memory
 * @author Peixinho
 */
void* iftRealloc(void* p, size_t size);

/**
 * @brief  Verify memory deallocation. 
 * @author Alexandre Falcao
 * @date   Feb 23, 2016
 * @ingroup Memory
 */
bool iftVerifyMemory(size_t MemDinInicial, size_t MemDinFinal);

/**
 * @brief Swap the content of two variables.
 * @warning In pointer, this function only swaps the pointers themselves, not their memory content.
 * @author Peixinho
 * @date Mar, 2016
 */
#define iftSwap(x, y) do { typeof(x) _IFT_SWAP_ = x; x = y; y = _IFT_SWAP_; } while (0)

/**
 * @brief Swap the content of two strings.
 * @author Peixinho
 * @date Mar, 2016
 */
void iftSwapString(char *a, char *b);


/**
 * @brief Shuffles an integer array by the modern version Fisher-Yates shuffle algorithm.
 * @details See https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * @param array Integer array
 * @param n Size of the array
 */
void iftShuffleIntArray(int* array, int n);

/**
 * @brief Allocates a Char Array (string) with <b>n</b> positions
 * @ingroup Memory
 * @{
 */
bool *iftAllocBoolArray(size_t n);
char *iftAllocCharArray(size_t n);
char *iftAllocString(size_t n);
uchar *iftAllocUCharArray(size_t n);
short *iftAllocShortArray(size_t n);
ushort *iftAllocUShortArray(size_t n);
int *iftAllocIntArray(size_t n);
uint *iftAllocUIntArray(size_t n);
long *iftAllocLongIntArray(size_t n);
#ifndef  __cplusplus
long long *iftAllocLongLongIntArray(size_t n);
#endif
ulong *iftAllocULongArray(size_t n);
#ifndef  __cplusplus
ullong *iftAllocULLongArray(size_t n);
#endif
float *iftAllocFloatArray(size_t n);
double *iftAllocDoubleArray(size_t n);
long double *iftAllocLongDoubleArray(size_t n);
iftComplex *iftAllocComplexArray(size_t n);

// Aligned memory allocation
uchar *iftAllocAlignedUCharArray(size_t n, size_t alignment);
int *iftAllocAlignedIntArray(size_t n, size_t alignment);
float *iftAllocAlignedFloatArray(size_t n, size_t alignment);
double *iftAllocAlignedDoubleArray(size_t n, size_t alignment);
/** @} */



/**
 * @brief Creates an iftUCharArray.
 * @author Thiago Vallin Spina
 * @date Mar 13, 2016
 * @ingroup Memory
 */
iftUCharArray *iftCreateUCharArray(size_t n);


/**
 * @brief Destroys an iftUCharArray.
 * @author Thiago Vallin Spina
 * @date Mar 13, 2016
 * @ingroup Memory
 */
void iftDestroyUCharArray(iftUCharArray **iarr);


/**
 * @brief Creates an iftIntArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
iftIntArray *iftCreateIntArray(size_t n);


/**
 * @brief Destroys an iftIntArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyIntArray(iftIntArray **iarr);


/**
 * @brief Creates an iftULongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
iftULongArray *iftCreateULongArray(size_t n);


/**
 * @brief Destroys an iftULongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
void iftDestroyULongArray(iftULongArray **iarr);


/**
 * @brief Creates an iftLongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
iftLongArray *iftCreateLongArray(size_t n);


/**
 * @brief Destroys an iftLongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
void iftDestroyLongArray(iftLongArray **iarr);


/**
 * @brief Creates an iftFloatArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
iftFloatArray *iftCreateFloatArray(size_t n);


/**
 * @brief Destroys an iftDblArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyFloatArray(iftFloatArray **darr);


/**
 * @brief Creates an iftDblArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
iftDblArray *iftCreateDblArray(size_t n);

/**
 * Copies a Double array, from data.
 * @author Peixinho
 * @date April, 2016
 */
iftDblArray *iftCopyDblArray(const double* array, size_t n);


/**
 * @brief Destroys an iftDblArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyDblArray(iftDblArray **darr);


/**
 * @brief Creates an Array of Strings with at most 2048 characters.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
iftStrArray *iftCreateStrArray(size_t n);


/**
 * @brief Destroys an iftStrArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyStrArray(iftStrArray **sarr);

/**
 * @brief Sets all values in <array_dst> to <val>.
 *
 * @param array_dst The array.
 * @param nelems The number of elements.
 * @param val value to set array
 * @ingroup Memory
 */
void iftSetFloatArray(float *array_dst, int nelems, float val);

/**
 * @brief Copies the float <array_src> to <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination array.
 * @param array_src The source array.
 * @param nelems The number of elements.
 * @ingroup Memory
 */
void iftCopyFloatArray(float *array_dst, float *array_src, int nelems);
/**
 * @brief Copies the double <array_src> to <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination array.
 * @param array_src The source array.
 * @param nelems The number of elements.
 * @ingroup Memory
 */
void iftCopyDoubleArray(double *array_dst, double *array_src, int nelems);
/**
 * @brief Copies the int <array_src> to <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination array.
 * @param array_src The source array.
 * @param nelems The number of elements.
 * @ingroup Memory
 */
void iftCopyIntArray(int *array_dst, const int *array_src, int nelems);
/**
 * @brief Concatenates two int arrays.
 *
 * @param array1 The first array.
 * @param n1 The number of elements in the first array.
 * @param array2 The first array.
 * @param n2 The number of elements in the second array.
 * @param nelems Returns the number of elements in the new array.
 * @return The concatenated arrays.
 */
int *iftConcatIntArray(int *array1, int n1, int *array2, int n2, int *nelems);


/**
 * @author Peixinho
 * @date Nov, 2016
 * @brief Creates a double array in the specified range
 * @param begin Starting value
 * @param end Ending value (not inclusive)
 * @param inc Increment ammount
 * @return
 */
iftDblArray *iftRange(double begin, double end, double inc);


/**
 * @brief Creates an integer array in the specified range.
 *
 * @note Inpired in range function from python.
 * 
 * @param begin Starting value
 * @param end Ending value (inclusive)
 * @param inc Increment ammount
 * @return Range Integer Array
 * 
 * @author Samuka
 * @date Jan 3, 2017
 */
iftIntArray *iftIntRange(int begin, int end, int inc);

/**
 * @brief Creates a string array with specific values.
 * @author Peixinho
 * @param n Number of values to search.
 * @param ... Sequence of values to search.
 */
iftStrArray *iftStrValues(int n, ...);

/**
 * @brief Creates a double array with specific values.
 * @warning Since this function uses undefined number of arguments, you HAVE to pass double values as parameter, not integers.
 * @author Peixinho
 * @param n Number of values to search.
 * @param ... Sequence of values to search.
 */
iftDblArray *iftDblValues(int n, ...);

/**
 * @author Peixinho
 * @date Nov, 2016
 * @brief Creates a double array in the specified geometric range
 * @param begin Starting value
 * @param end Ending value (not inclusive)
 * @param inc Geommetric increment ammount
 * @return
 */
iftDblArray *iftGeomRange(double begin, double end, double mul);

void iftCopyData(void* destination,const void* source,size_t n, size_t elementSize);


void iftCopySizeTArray(size_t *array_dst, const size_t *array_src, size_t nelems);

#endif //IFT_IFTMEMORY_H
