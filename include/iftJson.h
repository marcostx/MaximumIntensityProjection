/**
 * @file iftJson.h
 * @brief Json parser and manipulation.
 * @author Samuel Martins
 * @date Oct 8, 2015
 * @ingroup Json
 * 
 * @note An example of Json reading and writing and how to get elements from it can be found in
 * ift/demo/Miscellaneous/DataStrutuctures/iftTestJson1.c
 * @note An example of how to create a Json and add elements to it can be found in
 * ift/demo/Miscellaneous/DataStrutuctures/iftTestJson2.c
 * 
 * This is an adaption and extension from the cJSON parser:
 * https://github.com/kbranigan/cJSON
 */
/*
  Copyright (c) 2009 Dave Gamble

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#ifndef _IFT_JSON_H_
#define _IFT_JSON_H_

#ifdef __cplusplus
extern "C" {
#endif


#include "iftCommon.h"
#include "iftMatrix.h"
#include "iftFile.h"

/********************** cJSON Headers *************************/
/**
 * @brief Json Data Types.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 */
typedef enum {
    IFT_JSON_FALSE,
    IFT_JSON_TRUE,
    IFT_JSON_NULL,
    IFT_JSON_NUMBER,
    IFT_JSON_STRING,
    IFT_JSON_ARRAY,
    IFT_JSON_DICT
} iftJsonType;


/**
 * @brief Adapted cJson Structure.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * This Json structure is basically a Linked List. For sake of clarity, we have defined three types to reference it:
 * --> iftJson references a complete Json, that is, the root of the Json.
 * --> iftJDict references a dictionary of the Json. Note that the root of the Json is also a dictionary.
 * --> iftJNode references a node/element of the Json. Each node has a iftJsonType that indicates the type of the element.
 */
typedef struct _ift_json {
    /** Type of Json node, see above - they also store the BOOL values. */
    int type;
    /** Key of the json node. The root of the json has key = "(null)" */
    char *key;
    /** Integer value, if type == IFT_JSON_NUMBER */
    int int_val;
    /** Double value, if type == IFT_JSON_NUMBER */
    double dbl_val;
    /** String value, if type == IFT_JSON_STRING */
    char *str_val;

    /** An Array or Dictionary will have a first child and a last child pointer pointing to a chain of the items in the array/dict. */
    struct _ift_json *first_child, *last_child;
    /* Next/Prev allow to walk array/dict chains. */
    struct _ift_json *next, *prev;
} iftJson, iftJDict, iftJNode;


/**
 * @brief Reads a Json file from the disk.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * Reads a Json file from the disk. If the Json is bad formatted, an Error message is going to be raised.
 *
 * @param json_pathname The pathname from the Json to be read.
 * @return The iftJson.
 */
iftJson *iftReadJson(const char *json_path);


/**
 * @brief Destroys an iftJson.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json Reference to the iftJson to be destroyed.
 */
void iftDestroyJson(iftJson **json);


/**
 * @brief Prints an iftJson following a printing format.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json iftJson to be printed.
 */
void iftPrintJson(iftJson *json);


/**
 * @brief Prints an iftJson following a minified format.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json iftJson to be printed.
 */
void iftPrintJsonMinified(iftJson *json);


/**
 * @brief Writes an iftJson on the disk following a printing format.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json iftJson to be stored.
 */
void iftWriteJson(iftJson *json, const char *path);


/**
 * @brief Writes an iftJson on the disk following a minified format.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param json iftJson to be stored.
 */
void iftWriteJsonMinified(iftJson *json, const char *path);


/**
 * @brief Copies a Json (or a JDict or a JNode).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * @param json The iftJson to be copied.
 * @return The copied iftJson.
 */
iftJson *iftCopyJson(iftJson *json);

/**
 * @brief Copies a Json as a string that can be used for display and/or printing to a file.
 * @author Thiago V. Spina
 * @date Oct 28, 2016
 *
 * @param json The iftJson to be copied.
 * @return A string representing the json file.
 */
char *iftCopyJsonAsString(iftJson *json);

/**
 * @brief Returns an iftJNode REFERENCE from the key <key>. BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 * @sa iftGetDictReference
 *
 * Returns an iftJNode reference from the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jnode> outside, the Json will be inconsistent.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the an intermediate key is not a dict, an Error message is going to be raised.
 *
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The reference of the iftJNode required.
 */
iftJNode *iftGetJNodeReference(iftJson *json, const char *key);


/**
 * @brief Returns an iftJNode (COPYINT IT) from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 * @sa iftGetDict
 *
 * Returns an iftJNode (COPYINT IT) from the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the an intermediate key is not a dict, an Error message is going to be raised.
 *
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The iftJNode required.
 */
iftJNode *iftGetJNode(iftJson *json, const char *key);


/**
 * @brief Returns a boolean value from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a boolean value from the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Boolean, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The boolean value required.
 */
unsigned int iftGetJBool(iftJson *json, const char *key);


/**
 * @brief Returns a NULL pointer from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a NULL pointer from the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not "null", an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return NULL, if the value for the key is "null".
 */
void *iftGetJNull(iftJson *json, const char *key);


/**
 * @brief Returns an Integer value from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns an Integer value with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key do not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an Integer, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Integer value required.
 */
int iftGetJInt(iftJson *json, const char *key);


/**
 * @brief Returns a double value from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a double value with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Double, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The double value required.
 */
double iftGetJDouble(iftJson *json, const char *key);


/**
 * @brief Returns a string value (COPYING IT) from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a string value (COPYING IT) with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a String, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The string value required.
 */
char *iftGetJString(iftJson *json, const char *key);


/**
 * @brief Returns the reference of the iftJDict with the key <key>. BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns the reference of the iftJDict with the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jnode> outside, the Json will be inconsistent.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an iftJDict, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The reference of the required iftJDict.
 */
iftJDict *iftGetJDictReference(iftJson *json, const char *key);


/**
 * @brief Returns the iftJDict (COPYING IT) with the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns the the iftJDict (COPYING IT) with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an iftJDict, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The required iftJDict.
 */
iftJDict *iftGetJDict(iftJson *json, const char *key);


/**
 * @brief Returns an Integer Array from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns an Integer Array with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an Integer Array, an Error message is also going to be raised.
 * If some value of the array is not an Integer, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Integer Array required.
 */
iftIntArray *iftGetJIntArray(iftJson *json, const char *key);


/**
 * @brief Returns a Double Array from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a Double Array with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Double Array, an Error message is also going to be raised.
 * If some value of the array is not a Double, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Double Array required.
 */
iftFloatArray *iftGetJFloatArray(iftJson *json, const char *key);


/**
 * @brief Returns a String Array from the key <key>.
 *
 * @author Samuel Martins
 * @date October 19, 2015
 *
 * Returns a String Array with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a String Array, an Error message is also going to be raised.
 * If some value of the array is not a String, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The String Array required.
 */
iftStrArray *iftGetJStringArray(iftJson *json, const char *key);


/**
 * @brief Returns an Integer Matrix from the key <key>.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Returns an Integer Matrix with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not an Integer Matrix, an Error message is also going to be raised.
 * If some value of the matrix is not an Integer, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Integer Matrix required.
 */
iftIntMatrix *iftGetJIntMatrix(iftJson *json, const char *key);


/**
 * @brief Returns a Double Matrix from the key <key>.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Returns a Double Matrix with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a Double Matrix, an Error message is also going to be raised.
 * If some value of the matrix is not a Double, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The Double Matrix required.
 */
iftMatrix *iftGetJFloatMatrix(iftJson *json, const char *key);


/**
 * @brief Returns a String Matrix (COPYING each string) from the key <key>.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Returns a String Matrix (COPYING each string) with the key <key>.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 * If the value for the key <key> is not a String Matrix, an Error message is also going to be raised.
 * If some value of the matrix is not a String, an Error message is also going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @return The String Matrix required.
 */
iftStringMatrix *iftGetJStringMatrix(iftJson *json, const char *key);


// TODO: iftGetJDictArray


/**
 * @brief Creates a Json root which consists of an only Json node with its elements pointing to NULL.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * Creates a Json root which consists of an only Json node with its elements pointing to NULL.
 * It is necessary to add new elements to Json.
 *
 * @return The created iftJson root.
 */
iftJson *iftCreateJsonRoot(void);


/**
 * @brief Creates an EMPTY Json Node with its elements pointing to NULL (Similar to iftCreateJsonRoot()).
 *
 * @author Samuel Martins
 * @date October 20, 2015
 * @sa iftCreateJsonRoot()
 *
 * @return The created iftJNode.
 */
iftJNode *iftCreateJNode();


/**
 * @brief Creates a Json Node with type IFT_JSON_ARRAY and its elements pointing to NULL.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * @return The created Json Array Node.
 */
iftJNode *iftCreateJArray();


/**
 * @brief Creates a JDict with type IFT_JSON_DICT and its elements pointing to NULL (Similar to iftCreateJsonRoot()).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 * @sa iftCreateJsonRoot()
 *
 * @return The created iftJDict.
 */
iftJDict *iftCreateJDict();


/**
 * @brief Adds a JNode REFERENCE to a Json <json> in the key <key>. BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds a JNode <jnode> REFERENCE to Json in the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jnode> outside, the Json will be inconsistent.
 * If the key already exists already exists, it is going to be destroyed and the new JNode will be assigned instead.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @param jnode The reference of the JNode to added.
 */
void iftAddJNodeReferenceToJson(iftJson *json, const char *key, iftJNode *jnode);


/**
 * @brief Adds a JNode (COPYING IT) to a Json <json> in the key <key>.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds a JNode <jnode> (COPYING IT) to Json in the key <key>.
 * If the key already exists already exists, it is going to be destroyed and the new JNode will be assigned instead.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be used.
 * @param key The key to be used.
 * @param jnode The JNode to be copied and added to Json.
 */
void iftAddJNodeToJson(iftJson *json, const char *key, iftJNode *jnode);


/**
 * @brief Adds or sets a new Boolean value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets a new Boolean value <bool_val> to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param bool_val The Boolean value to be added into.
 * @param key The key to be used.
 */
void iftAddBoolToJson(iftJson *json, const char *key, unsigned int bool_val);


/**
 * @brief Adds or sets a new Null value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets a new Null value to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param key The key to be used.
 */
void iftAddNullToJson(iftJson *json, const char *key);


/**
 * @brief Adds or sets a new integer value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets a new integer value <int_val> to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param int_val The integer to be added into.
 * @param key The key to be used.
 */
void iftAddIntToJson(iftJson *json, const char *key, int int_val);


/**
 * @brief Adds or sets a new double value to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets a new double value <dbl_val> to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param dbl_val The double to be added into.
 * @param key The key to be used.
 */
void iftAddDoubleToJson(iftJson *json, const char *key, double dbl_val);


/**
 * @brief Adds or sets a new string value (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets a new string value <str_val> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param str_val The string to be added into.
 * @param key The key to be used.
 */
void iftAddStringToJson(iftJson *json, const char *key, const char *str_val);


/**
 * @brief Adds or sets the REFERENCE of a JDict to a Json (or a JDict). BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY.
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets the REFERENCE of the JDict <jdict> to the iftJson or iftJDict <json> in the key <key>.
 * BE CAREFUL, BECAUSE IT'S A REFERENCE AND NOT A COPY. Then, if you destroy the <jdict> outside, the Json will be inconsistent.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param jdict The reference of the iftJDict to be added.
 * @param key The key to be used.
 */
void iftAddJDictReferenceToJson(iftJson *json, const char *key, iftJDict *jdict_ref);


/**
 * @brief Adds or sets a JDict (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 22, 2015
 *
 * Adds or sets the JDict <jdict> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param jdict The iftJDict to be added.
 * @param key The key to be used.
 */
void iftAddJDictToJson(iftJson *json, const char *key, iftJDict *jdict);


/**
 * @brief Adds or sets an Integer Array (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Integer Array <iarr> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param iarr The iftIntArray to be added.
 * @param key The key to be used.
 */
void iftAddIntArrayToJson(iftJson *json, const char *key, const iftIntArray *iarr);


/**
 * @brief Adds or sets a Double Array (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Double Array <darr> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param darr The iftDblArray to be added.
 * @param key The key to be used.
 */
void iftAddFloatArrayToJson(iftJson *json, const char *key, iftFloatArray *darr);


/**
 * @brief Adds or sets a String Array (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the String Array <sarr> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param sarr The iftStrArray to be added.
 * @param key The key to be used.
 */
void iftAddStringArrayToJson(iftJson *json, const char *key, iftStrArray *sarr);


/**
 * @brief Adds or sets an Integer Matrix (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Integer Matrix <iM> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param iM The iftIntArray to be added.
 * @param key The key to be used.
 */
void iftAddIntMatrixToJson(iftJson *json, const char *key, iftIntMatrix *iM);


/**
 * @brief Adds or sets a Double Matrix (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the Double Matrix <dM> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param dM The iftDblArray to be added.
 * @param key The key to be used.
 */
void iftAddDoubleMatrixToJson(iftJson *json, const char *key, iftMatrix *dM);


/**
 * @brief Adds or sets a String Matrix (COPYING IT) to a Json (or a JDict).
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds or sets the String Matrix <sM> (COPYING IT) to the iftJson or iftJDict <json> in the key <key>.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the key already exists, the new value is assigned in the key.
 * If the key already exists and the element is a String, it will be deallocated.
 * If the key already exists and the element is a JDict or an Array, its children will be destroyed.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param sM The iftStrArray to be added.
 * @param key The key to be used.
 */
void iftAddStringMatrixToJson(iftJson *json, const char *key, iftStringMatrix *sM);


/**
 * @brief Adds (APPENDING IT) a new Integer into Json Array Node (node with type = IFT_JSON_ARRAY). BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Integer into Json Array Node (node with type = IFT_JSON_ARRAY).
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * No key is necessary.
 * The value is inserted at the end of the array.
 * If the JNode used is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 *
 * @param jarray The JArray node where the value will be added.
 * @param int_val The Integer to be added.
 */
void iftAddIntToJArrayNode(iftJNode *jarray, int int_val);


/**
 * @brief Adds (APPENDING IT) a new Double into Json Array Node (node with type = IFT_JSON_ARRAY). BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Double into Json Array Node (node with type = IFT_JSON_ARRAY).
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * No key is necessary.
 * The value is inserted at the end of the array.
 * If the JNode used is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 *
 * @param jarray The JArray node where the value will be added.
 * @param dbl_val The Double to be added.
 */
void iftAddDoubleToJArrayNode(iftJNode *jarray, double dbl_val);


/**
 * @brief Adds (APPENDING IT) a new String into Json Array Node (node with type = IFT_JSON_ARRAY). BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new String into Json Array Node (node with type = IFT_JSON_ARRAY).
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * No key is necessary.
 * The value is inserted at the end of the array.
 * If the JNode used is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 *
 * @param jarray The JArray node where the value will be added.
 * @param str_val The String to be added.
 */
void iftAddStringToJArrayNode(iftJNode *jarray, const char *str_val);


/**
 * @brief Adds (APPENDING IT) a new Integer into Json Array with key <key>. BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Integer into Json Array with key <key>.
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * The value is inserted at the end of the array.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the JNode with key <key> is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The Json (or JDict) to be used.
 * @param key The key of the JArray to be used.
 * @param int_val The Integer to be added.
 */
void iftAddIntToJArray(iftJson *json, const char *key, int int_val);


/**
 * @brief Adds (APPENDING IT) a new Double into Json Array with key <key>. BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new Double into Json Array with key <key>.
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * The value is inserted at the end of the array.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the JNode with key <key> is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The Json (or JDict) to be used.
 * @param key The key of the JArray to be used.
 * @param dbl_val The Double to be added.
 */
void iftAddDoubleToJArray(iftJson *json, const char *key, double dbl_val);


/**
 * @brief Adds (APPENDING IT) a new String into Json Array with key <key>. BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Adds (APPENDING IT) a new String into Json Array with key <key>.
 * BE CAREFUL BECAUSE THERE IS NO DATATYPE CHECKING IN THE ARRAY. Then, a value can be added in an array with elements of different datatype.
 * The value is inserted at the end of the array.
 * The <json> does not need to be a Json root, a JDict also works.
 * If the JNode with key <key> is not a JArray (type = IFT_JSON_ARRAY), an Error Message is going to be raised.
 * It is possible to pass multiple keys, according to the Json hierarchy. For that, use the following format:
 * "key1:key2:key3:keyn".
 *
 * If an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The Json (or JDict) to be used.
 * @param key The key of the JArray to be used.
 * @param str_val The String to be added.
 */
void iftAddStringToJArray(iftJson *json, const char *key, const char *str_val);


/**
 * @brief Deletes (deallocates) the JNode with the key <key> and its descendant from a Json or a JDict.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 *
 * Deletes (deallocates) the JNode with the key <key> and its descendant from a Json or a JDict.
 * If the key or an intermediate key does not exist, an Error message is going to be raised.
 *
 * @param json The iftJson to be read.
 * @param key The key from the node to be deleted.
 */
void iftDeleteJNode(iftJson *json, const char *key);


bool iftJsonContainKey(const char *key, iftJson *json);
/*********************************************************************************/




#ifdef __cplusplus
}
#endif

#endif //_IFT_JSON_H_
