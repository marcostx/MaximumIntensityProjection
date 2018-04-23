/**
 * @file iftList.h
 * @brief Integer Doubly Linked List definitions and functions.
 *
 * An example can be found in ift/demo/Miscellaneous/iftTestList.c
 *
 * @author Samuel Martins
 * @date September 18, 2015
 */

#ifndef _IFT_LIST_H_
#define _IFT_LIST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"

/**
 * @brief Node of the Integer Doubly Linked List.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 */
typedef struct _ift_node {
    int elem;
    struct _ift_node *previous;
    struct _ift_node *next;
} iftNode;


/**
 * @brief Doubly Linked List of Integers.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 */
typedef struct _ift_list {
    /** Number of Nodes of the Integer Doubly Linked List */
    int n;
    /** Pointer to the Head (first) Node */
    iftNode *head;
    /** Pointer to the Tail (last) Node */
    iftNode *tail;
} iftList;



/**
 * @brief Creates an iftList.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 */
iftList *iftCreateList();


/**
 * @brief Destroys an iftList.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 *
 * @param L Reference to the iftList to be destroyed.
 */
void iftDestroyList(iftList **L);


/**
 * @brief Insert an integer into the head of the iftList.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 *
 * Insert the integer <elem> into the head of the iftList.
 * If the iftList is not allocated, an error will be raised.
 *
 * @param L The iftList.
 * @param elem The integer to be inserted into the head of the iftList.
 */
void iftInsertListIntoHead(iftList *L, int elem);


/**
 * @brief Insert an integer into the tail of the iftList.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 *
 * Insert the integer <elem> into the tail of the iftList.
 * If the iftList is not allocated, an error will be raised.
 *
 * @param L The iftList.
 * @param elem The integer to be inserted into the tail of the iftList.
 */
void iftInsertListIntoTail(iftList *L, int elem);


/**
 * @brief Removes the Head element from the iftList and removes its node.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 *
 * Removes the Head integer node of the iftList and removes its node.
 * If there is no element, the NIL value is returned.
 * If the iftList is not allocated, an error will be raised.
 *
 * @param L The iftList.
 * @return The Head integer from the iftList.
 */
int iftRemoveListHead(iftList *L);


/**
 * @brief Removes the Tail element from the iftList and removes its node.
 *
 * @author Samuel Martins
 * @date September 18, 2015
 *
 * Removes the Tail integer node of the iftList and removes its node.
 * If there is no element, the NIL value is returned.
 * If the iftList is not allocated, an error will be raised.
 *
 * @param L The iftList.
 * @return The Tail integer from the iftList.
 */
int iftRemoveListTail(iftList *L);


/**
 * @brief Prints an iftList (from the Node to the Tail node).
 *
 * @author Samuel Martins
 * @date September 18, 2015
 *
 * Prints an iftList from the Head to Tail node.
 *
 * @param L The iftList.
 */
void iftPrintList(iftList *L);


/**
 * @brief Prints an iftList in Reveresed Order (from the Tail to the Head node).
 *
 * @author Samuel Martins
 * @date September 18, 2015
 *
 * Prints an iftList from the Tail to Head node.
 *
 * @param L The iftList.
 */
void iftPrintReversedList(iftList *L);




#ifdef __cplusplus
}
#endif

#endif //_IFT_LIST_H_

