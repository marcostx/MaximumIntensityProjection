/**
 * @file iftCommonDataTypes.h
 * @brief Definition of Common Data Types.
 * @author Samuel Martins
 * @date Feb 29, 2016
 * @ingroup DataTypes
 *
 * @note Programs:
 * * @note iftTestGVal.c = It shows how to use the Generic Value (iftGVal)
 */

#ifndef IFT_COMMON_DATATYPES_H
#define IFT_COMMON_DATATYPES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ctype.h>
#include <limits.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iftDialog.h>

/**
 * @brief Limits for common data types.
 * @ingroup DataTypes
 * @{
 */
#define IFT_INFINITY_INT       INT_MAX
#define IFT_INFINITY_INT_NEG   INT_MIN
#define IFT_INFINITY_LONG      LONG_MAX
#define IFT_INFINITY_LONG_NEG  LONG_MIN
#define IFT_INFINITY_ULONG     ULONG_MAX
#define IFT_INFINITY_FLT       FLT_MAX
#define IFT_INFINITY_FLT_NEG  -FLT_MAX
#define IFT_INFINITY_DBL       DBL_MAX
#define IFT_INFINITY_DBL_NEG  -DBL_MAX
#define IFT_INFINITY_LDBL      LDBL_MAX
#define IFT_INFINITY_LDBL_NEG -LDBL_MAX
/**
 * @}
 */

/**
 * @brief This value should be used when allocating static strings in stack. It helps to prevent buffer overflows.
 * @author Thiago Vallin Spina
 * @date Feb 23, 2016
 * @ingroup DataTypes
 */
#define IFT_STR_DEFAULT_SIZE BUFSIZ

/**
 * @brief A "NIL" value for unsigned long.
 * @author Samuel Martins
 * @date Feb 29, 2016
 * @ingroup DataTypes
 */
#define IFT_ULONG_NIL IFT_INFINITY_ULONG

typedef struct ift_dict iftDict;
typedef struct ift_str_array iftStrArray;
typedef struct ift_dbl_array iftDblArray;

/**
 * @brief Short names for data types.
 * @ingroup DataTypes
 * @{
 */
typedef struct timeval timer;
typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;
#ifndef  __cplusplus
typedef long long llong;
typedef unsigned long long ullong;
#endif
/**
 * @}
 */


/**
 * @brief Flags to identify Several C Data Types.
 * @author Samuel Martins
 * @date Jan 31, 2016
 * @ingroup DataTypes
 */
typedef enum ift_cdata_type {
    IFT_UNTYPED, IFT_BOOL_TYPE, IFT_CHAR_TYPE, IFT_UCHAR_TYPE, IFT_STR_TYPE,
    IFT_INT_TYPE, IFT_UINT_TYPE, IFT_FLT_TYPE, IFT_LONG_TYPE, IFT_ULONG_TYPE,
    IFT_DBL_TYPE, IFT_DBL_ARRAY_TYPE, IFT_STR_ARRAY_TYPE, IFT_DICT_TYPE, IFT_PTR_TYPE
} iftCDataType;


/**
 * @brief Basic datatypes useful for Multi-Band kernels and images.
 * @ingroup DataTypes
 */
typedef struct ift_band {
    float *val;
} iftBand;


/**
 * @brief 3D vector/point in the real space.
 * @ingroup DataTypes
 */
typedef struct ift_vector {
    float x, y, z;
} iftVector, iftPoint;


/**
 * @brief 3D point in the integer space.
 * @ingroup DataTypes
 */
typedef struct ift_voxel {
    int x, y, z;
} iftVoxel, iftSize;


/**
 * @brief Image Bounding Box.
 * @author Samuel Martins
 * @date Mar 2, 2016 
 * @ingroup DataTypes
 */
typedef struct ift_bounding_box {
    /** The leftmost and topmost point of the Bounding Box */
    iftVoxel begin;
    /** The rightmost and bottommost point of the Bounding Box */
    iftVoxel end;
} iftBoundingBox;


/**
 * @brief Complex number.
 * @ingroup DataTypes
 */
typedef struct ift_complex {
    double r;
    double i;
} iftComplex;


/**
 * @brief An abstraction of Value datatype. It holds the value and its datatype. Each value field shares the same memory address.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
typedef struct ift_gval {
    /** C Data type of the assigned value. */
    iftCDataType type;
    /**
     * Assigned value. All fields shares the same memory address.
     * Since that <b>long</b> type is has memory size greater or equal to all integer datatypes (and their unsigned versions),
     * we've just used an only field to represent integer number.
     * The same holds to double and float, and to pointers.
     * 
     * @warning THERE IS NO SUPPORT TO LONG LONG AND LONG DOUBLE.
     */
    union {
        bool   bool_val;
        char   char_val;
        uchar  uchar_val;
        char*  str_val;
        long   long_val;
        ulong  ulong_val;
        double dbl_val;
        iftStrArray* str_array_val;
        iftDblArray* dbl_array_val;
        iftDict* dict_val;
        void*  ptr_val;
    };
} iftGVal;



////////////// GENERIC TYPE MANIPULATION //////////////
/**
 * @brief A Generic Function to Create a Generic Value (iftGVal) from a value of any datatype.
 * @author Samuel Martins
 * @date Jan 18, 2016
 * @ingroup DataTypes
 * 
 * Pointer types MUST BE CAST to (void*).
 * Just pass a value VAL that the correct function will be called.
 * A demo can be found in demo/Miscellaneous/iftTestGVal.c
 */
#define iftCreateGVal(VAL) _Generic((VAL), \
    bool:        iftInitBoolGVal, \
    char:        iftInitCharGVal, \
    uchar:       iftInitUCharGVal, \
    char*:       iftInitStrGVal, \
    const char*: iftInitStrGVal, \
    short:       iftInitLongGVal, \
    ushort:      iftInitULongGVal, \
    int:         iftInitLongGVal, \
    uint:        iftInitULongGVal, \
    long:        iftInitLongGVal, \
    ulong:       iftInitULongGVal, \
    float:       iftInitDblGVal, \
    double:      iftInitDblGVal, \
    iftDblArray*: iftInitDblArrayGVal, \
    iftStrArray*: iftInitStrArrayGVal, \
    iftDict*: iftInitDictGVal, \
    void*:       iftInitPtrGVal, \
    iftGVal:     iftCopyGVal, \
    default:     iftInitPtrGVal \
    )(VAL)

/**
 * @brief Initializes a generic value with a bool.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitBoolGVal(bool val);

/**
 * @brief Initializes a generic value with a char. Literal char must use the cast (char) before.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitCharGVal(char val);


/**
 * @brief Initializes a generic value with an unsigned char.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitUCharGVal(uchar val);


/**
 * @brief Initializes a generic value COPYING a string.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitStrGVal(const char *val);


/**
 * @brief Initializes a generic value with a long.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitLongGVal(long val);


/**
 * @brief Initializes a generic value with an unsigned long.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitULongGVal(ulong val);


/**
 * @brief Initializes a generic value with a double.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitDblGVal(double val);

iftGVal iftInitDblArrayGVal(iftDblArray* array);

iftGVal iftInitStrArrayGVal(iftStrArray* array);

iftGVal iftInitDictGVal(iftDict* dict);


/**
 * @brief Initializes a generic value with a void pointer.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
iftGVal iftInitPtrGVal(void *val);


/**
 * @brief Returns the name (literal string) of the CDataType enumerator.
 * @author Samuel Martins
 * @date Feb 4th, 2016
 * @ingroup DataTypes
 */
const char *iftCDataTypeToString(iftCDataType datatype);


/**
 * @brief Converts the GVal to literal string.
 * @author Samuel Martins
 * @date Feb 7th, 2016
 * @ingroup DataTypes
 */
char *iftGValToString(iftGVal gval);


/**
 * @brief Compares two generic values. If they are equals, it returns true.
 * @author Samuel Martins
 * @date Jan 19, 2016
 * @ingroup DataTypes
 */
bool iftCompareGVal(iftGVal val1, iftGVal val2);


/**
 * @brief Copies a Generic Value. It the input is a string, a new string is allocated and copied.
 * @author Samuel Martins
 * @date Feb 18, 2016
 * @ingroup DataTypes
 */
iftGVal iftCopyGVal(iftGVal gval);


/**
 * @brief Gets the boolean value from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a boolean value.
 */
bool iftGetBoolVal(iftGVal gval);


/**
 * @brief Gets the char value from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a char value.
 */
char iftGetCharVal(iftGVal gval);


/**
 * @brief Gets the unsigned char value from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not an unsigned char value.
 */
uchar iftGetUCharVal(iftGVal gval);


/**
 * @brief Gets the string value (copying it) from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a string value.
 */
char *iftGetStrVal(iftGVal gval);


/**
 * @brief Gets the string value (without copying it) from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a string value.
 */
char *iftGetStrPtrVal(iftGVal gval);


/**
 * @brief Gets the constant string value (without copying it) from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a string value.
 */
const char *iftGetConstStrVal(iftGVal gval);


/**
 * @brief Gets the long value from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a long value.
 */
long iftGetLongVal(iftGVal gval);


/**
 * @brief Gets the unsigned long value from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not an unsigned long value.
 */
ulong iftGetULongVal(iftGVal gval);


/**
 * @brief Gets the double value from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a double value.
 */
double iftGetDblVal(iftGVal gval);

/**
 * @brief Gets the double array from gval.
 * @author Peixinho
 * @ingroup DataTypes
 */
iftDblArray * iftGetDblArrayVal(iftGVal gval);

/**
 * @brief Gets the string array from gval.
 * @author Peixinho
 * @ingroup DataTypes
 */
iftStrArray * iftGetStrArrayVal(iftGVal gval);

/**
 * @brief Gets the dictionary from gval.
 * @author Peixinho
 * @ingroup DataTypes
 */
iftDict * iftGetDictVal(iftGVal gval);


/**
 * @brief Gets the pointer value (void*) from gval.
 * @author Samuel Martins
 * @date Feb 12, 2016
 * @ingroup DataTypes
 * @exception gval is not a pointer (void*) value.
 */
void *iftGetPtrVal(iftGVal gval);
///////////////////////////////////////////////////////


/**
 * @brief Copies the src voxel to the destination voxel. Both references are passed.
 * @author Samuka
 * @date Dec 20, 2016
 * @ingroup DataTypes
 */
void iftCopyVoxel(iftVoxel *src, iftVoxel *dst);


#ifdef __cplusplus
}
#endif

#endif //IFT_COMMON_DATATYPES_H
