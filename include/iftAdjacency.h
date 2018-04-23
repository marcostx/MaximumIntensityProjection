#ifndef IFT_ADJACENCY_H_
#define IFT_ADJACENCY_H_

/**
 * @file iftAdjacency.h
 * @brief Structs and function prototypes for adjacency relation.
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * 
 */

/** @addtogroup Image
 * @{ */

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"


/**
 * @brief Defines a neighborhood around a voxel.
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 *
 */
typedef struct ift_adjrel {
    int *dx, *dy, *dz;
    /* displacements to achieve the n adjacent voxels. */
    int n; /* number of adjacent voxels. */
} iftAdjRel;

/**
 * @brief Defines a neighborhood around a voxel, with a faster Iteration, 
 * Disconsidering the image's border
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning Can only be used for images with same size. See also ::iftAdjRel.
 *
 */
typedef struct ift_fastadjrel {
    int n;
    /* number of adjacent voxels */
    int *dq;
    /* displacements to reach adjacent voxels for a given image */
    int bx, by, bz; /* sizes of the image's border to be disconsidered */
} iftFastAdjRel;

/**
 * @brief Allocates memory for a 3D adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 *
 */
iftAdjRel *iftCreateAdjRel(int n); /*  */


/**
 * @brief Deallocates memory for a adjacency relation object.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 *
 */
void iftDestroyAdjRel(iftAdjRel **A);

/**
 * @brief Creates a 3D ball of radius @a r as adjacency relation.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning This function is used for 3D images.
 * @param r Sphere radius.
 */
iftAdjRel *iftSpheric(float r);

/**
 * @brief Creates a 3D half-ball of radius @a r as adjacency relation, in the corresponding axis and direction.
 * This adjacency is useful for segmenting a volume in a single direction, e.g., a video-volume where z axis is time.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Sphere radius.
 * @param axis Axis to be considered.
 * @param direction The direction to follow in @a axis. (-1, 1)
 * @return The hemispheric adjacency.
 * */
iftAdjRel *iftHemispheric(float r, char axis, int direction);


/**
 * @brief Creates a 3D ball surface of radius r as adjacency relation.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning This function is used for 3D images.
 */
iftAdjRel *iftSphericEdges(float r);

/**
 * @brief Creates a 2D ball of radius r on the xy plane as adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning This function is used for 2D images.
 * @param r Circle radius.
 * */
iftAdjRel *iftCircular(float r);


/**
 * @brief Creates a circunference of radius r on the xy plane as adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning This function is used for 2D images.
 * @param r Circle radius.
 * */
iftAdjRel *iftCircularEdges(float r);

/**
 * @brief Creates a 2D ball of radius @a r on the xy plane as adjacency relation for contour pixel labeling.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Circle radius.
 **/
iftAdjRel *iftClockCircular(float r);

/**
 * @brief Creates a 2D adjacency relation with only 
 * the right-hand side of the input adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Circle radius.
 **/
iftAdjRel *iftRightSide(iftAdjRel *A, float r);

/**
 * @brief Creates a 2D adjacency relation with only 
 * the left-hand side of the input adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Circle radius.
 **/
iftAdjRel *iftLeftSide(iftAdjRel *A, float r);

/**
 * @brief Creates a rectangle adjacency with specified dimensions.
 * @warning This function is used for 2D images.
 * @param xsize Rectangle width.
 * @param ysize Rectangle height.
 * @return The rectangular adjacency.
 */
iftAdjRel *iftRectangular(int xsize, int ysize);

/**
 * @brief Creates a Cuboid adjacency with specified dimensions.
 *
 * @warning This function is used for 3D images.
 *
 * @param xsize Rectangle width.
 * @param ysize Rectangle height.
 * @param zsize Rectangle depth.
 * @return The cuboid adjacency.
 */
iftAdjRel *iftCuboid(int xsize, int ysize, int zsize);

/**
 * @brief Creates a new copy of the given adjacency.
 * @param A Adjacency to be copied.
 * @return Copy of adjacency @a A.
 */
iftAdjRel *iftCopyAdjacency(iftAdjRel *A);

iftFastAdjRel *iftCreateFastAdjRel(iftAdjRel *A, int *tby, int *tbz);

/* create an adjacency relation to speed up implementations for a given image by computing the displacements to the adjaceny voxels based on the look-up tables tby and tbz of the image. The fast implementation must disconsider the image's border */

void iftDestroyFastAdjRel(iftFastAdjRel **F);

void iftMaxAdjShifts(const iftAdjRel *A, int *dx, int *dy, int *dz);

void iftWriteAdjRel(iftAdjRel *A, char *filename);

iftAdjRel *iftReadAdjRel(char *filename);

/** Read an Adjacency Relation from a Binary File */
iftAdjRel *iftReadAdjRelBinFile(const char *path);

/** Writes an Adjacency Relation from a Binary File */
void iftWriteAdjRelBinFile(const iftAdjRel *A, const char *path);

iftVoxel iftGetAdjacentVoxel(const iftAdjRel *A, iftVoxel u, int adj);

/**
 * @brief Check if the adjacency is suitable for 3D images.
 * @param A The checked adjacency.
 * @return 1 if it is a 3D adjacency, 0 otherwise.
 */
int iftIsAdjRel3D(iftAdjRel *A);


/**
 * @brief Given an adjacency relation <A>, preferably created with iftCircular or iftSpheric, this function extracts
 * all displacements whose magnitude is within the closed interval [radius_min, radius_max] into a
 * new adjacency relation. Hence, this function may be used to create a shell (ring) adjacency relation. NOTE THAT: this
 * function generalizes the concept of ring adjacency to any kind of adjacency relations other than iftCircular and
 * iftSpheric. We could replace it by a simpler functiond called iftCircularRing (2D) and iftSphericRing (3D), but we
 * implemented this one instead for the sake of generalization. Moreover, by filtering an existing adjacency relation
 * we ensure that <A> may be broken into several non-overlapping rings, if we call iftExtractAdjacencyShell with
 * different non-overlapping intervals.
 *
 * @author Thiago Vallin Spina
 * @date Feb 10, 2016
 *
 * @param A The input adjacency relation (preferably a circular or spherical one).
 * @param radius_min The minimum length that the displacement must have to be part of the new adjacency relation.
 * @param radius_max The maximum length that the displacement must have to be part of the new adjacency relation.
 * @param include_center This parameter forces that displacement (0,0,0) to be added in the new adjacency relation, since
 * it may not be contained within [radius_min, radius_max]
 *
 * @return The filtered adjacency relation with displacement (0,0,0) in the first position, if required and/or part of
 * the interval.
 */
iftAdjRel *iftExtractAdjacencyShell(iftAdjRel *A, double radius_min, double radius_max, uchar include_center);

#ifdef __cplusplus
}
#endif

/** @} */

#endif
