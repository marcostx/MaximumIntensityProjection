/**
 * @file iftRegistration.h
 * @brief Image Registration module.
 * @author Samuel Martins
 * @date Nov 5th, 2015
 * @ingroup Registration
 *
 * @note Examples of Image Registration can be found in ift/demo/Registration/iftRegisterImageByElastix.c and
 * ift/demo/Registration/iftRegisterImageSetByElastix.c.
 * note Examples of Image Registration given deformation fields can be found in ift/demo/Registration/iftTransformImageByTransformix.c and
 * ift/demo/Registration/iftTransformImageSetByTransformix.c
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef IFT_REGISTRATION_H
#define IFT_REGISTRATION_H

#include "iftCommon.h"
#include "iftDict.h"
#include "iftFile.h"
#include "iftImage.h"
#include "iftMatrix.h"
#include "iftObjectModels.h"
#include "iftRadiometric.h"
#include "iftRepresentation.h"
#include "iftSeeds.h"
#include "iftUtil.h"


#define IFT_ELASTIX_AFFINE_DEF_FIELDS_SUFFIX "_DefFields_Affine.txt"
#define IFT_ELASTIX_BSPLINE_DEF_FIELDS_SUFFIX "_DefFields_BSpline.txt"


/**
 * @brief Registration Algorithms which are support on libIFT.
 * @author Samuel Martins
 * @date Nov 14, 2015
 * @ingroup Registration
 */
typedef enum {
    IFT_ELASTIX
} iftRegistrationAlg;


/**
 * @brief Transformation/Mapping Algorithms which are support on libIFT.
 * @author Samuel Martins
 * @date Nov 14, 2015
 * @ingroup Registration
 */
typedef enum {
    IFT_TRANSFORMATION_MATRIX,
    IFT_TRANSFORMIX
} iftTransformationAlg;


/**
 * @brief General Struct for Registration Parameters.
 * @author Samuel Martins
 * @date Nov 14, 2015
 * @ingroup Registration
 */
typedef struct ift_registration_params {
    /** Chosen Registration Algorithm. */
    iftRegistrationAlg alg;
    /** Dictionary with the Registration Method Parameters. */
    iftDict *params;
} iftRegistrationParams;


/**
 * @brief Setups the Parameters of the Elastix Registration.
 * @author Samuel Martins
 * @date November 14, 2015
 * @ingroup Elastix
 * 
 * @param  fixed_mask_path         Pathname from the mask for the fixed image (it can be NULL).
 * @param  affine_params_pathname  Pathname from the Affine Parameter File.
 * @param  bspline_params_pathname Pathname from the BSpline Parameter File.
 * @param  img_depth               Depth of the images to be registered.
 * @return                         The Elastix Registration Parameters.
 */
iftRegistrationParams *iftSetupElastixRegistation(const char *fixed_mask_path, const char *affine_params_pathname,
                                                  const char *bspline_params_pathname, int img_depth,
                                                  const char *out_basename);


/**
 * @brief Destroys Registration Parameters.
 * @author Samuel Martins
 * @date Nov 15, 2015
 * @ingroup Registration
 */
void iftDestroyRegistrationParams(iftRegistrationParams **reg_params);


/**
 * @brief Register a Moving Image onto a Fixed Image.
 * @author Samuel Martins
 * @date Nov 15, 2015
 * @ingroup Registration
 *
 * @param fixed_img Reference Image (space) where the Moving (target) Image will be registered.
 * @param moving_img Image to be registered.
 * @param reg_params_in Input Registration Parameters.
 * @param reg_params_out Reference of the Registration Parameters outputted from the Registration Method.
 * @return The Registered Image.
 *
 * @attention If a NULL pointer is passed to the parameter <code><b>reg_params_out</b></code>,
 * the resulting parameter is not returned, being deallocated instead.
 * 
 * @exception reg_params_in->alg does have a valid Registration Method.
 */
iftImage *iftRegisterImage(const iftImage *moving_img, const iftImage *fixed_img,
                           const iftRegistrationParams *reg_params_in,
                           iftRegistrationParams **reg_params_out);


/**
 * @brief Run the program iftRegisterImageByElastix and returns the resulting registered image.
 * 
 * @param  moving_img          Moving Image pathname.
 * @param  fixed_img           Fixed Image pathname.
 * @param  affine_params_path  Affine Parameter File.
 * @param  bspline_params_path BSpline Parameter File.
 * @param  out_def_fields_dir  If != NULL, it returns a directory with the resulting deformation fields files.
 * @return                     Resulting Registered Image.
 *
 * @author Samuka
 * @date Dec 28, 2016
 * @ingroup Registration
 */
iftImage *iftRunProgRegisterImageByElastix(const iftImage *moving_img, const iftImage *fixed_img,
                                           const char *affine_params_path, const char *bspline_params_path,
                                           char **out_affine_def_field, char **out_bspline_def_field);


/**
 * @brief Run the program iftRegisterSetImageByElastix.
 * @author Samuka
 * @date Jan 6, 2017
 * @ingroup Registration
 */
void iftRunProgRegisterImageSetByElastix(const char *moving_img_entry, const char *fixed_img_path, int img_depth,
                                         const char *affine_params_path, const char *bspline_params_path,
                                         const char *out_dir);


/**
 * @brief Run the program iftTransformImageSetByTransformix.
 * @author Samuka
 * @date Jan 9, 2017
 * @ingroup Registration
 */
void iftRunProgTransformImageSetByTransformix(const char *img_entry, const char *affine_def_fields_entry,
                                              const char *bspline_def_fields_entry, const char *out_dir);


/**
 * @brief Applies a Transformation/Mapping on to an Image.
 * @author Samuel Martins
 * @date Dec 7, 2015
 * @ingroup Registration
 *
 * @param img Image to be transformed/mapped.
 * @param transformation The transformation to be applied on to the image. It can be a Matrix, configuration file, etc...
 * @param alg Transformation/Registration alg.
 * @return The Transformed/Mapped Image.
 *
 * @note The parameter <code><b>transformation</b></code> will be casted according to the transformation alg.
 * @sa iftTransformLabelImage()
 */
iftImage *iftTransformImage(const iftImage *img, const void *transformation, iftTransformationAlg method);


/**
 * @brief Centralizes a set of Label Images in a same domain. The new images are written in an output directory.
 * @author Samuel Martins
 * @date Oct 7, 2015
 * @ingroup Registration
 *
 * Centralizes a set of Label Images into a same domain. The new images are saved in the <output_dir>.\n
 * In order to put all label images into a same domain, the following steps are executed for each image:\n
 * 1. Gets the Minimum Bounding Box (MBB) and its Image Geometric Center;
 * 2. Expands the MBB in order to become the Image Geometric Center the Central Point from the MBB.
 *
 * The new domain, which will fit all objects of the label images, will consist of the largest components
 * from all expanded MBBs of the images.
 *
 * @param label_imgs_paths A File Array with the pathnames from the label images to be centralized.
 * @param output_dir Pathname from the output directory where the new image will be saved/stored.
 */
void iftCentralizeLabelImages(const iftFileSet *label_imgs_paths, const char *output_dir);



void iftDeformDirectoryTransformix(char *masksDirectory, char *outDirectory);
void iftAlignImagesForSimilarity(char *dirIn, char *dirOut);


/**
 * @brief Register a cloud point in order to maximize the given score image.
 *
 * Given an array of points this function tries to apply affine transformations in order to better fit these points maximizing the score image.
 * The score image can be computed as the Inverted Euclidean Distance Transform as in iftEuclideanScoreImage(), or any other maximization score.
 *
 * @author Peixinho
 *
 * @date 12 Aug 2015
 *
 */
iftMatrix* iftShapeRegister(iftPoint* orig, int norig, iftImage* score);


/**
 * @brief Computes the inverse of Euclidean Distance Transform in a border Image.
 *
 * @param img The input binary border image.
 * @param decay The decay factor for the score, larger decay factors creates a score image that penalizes points too far from the border.
 *
 * @author Peixinho
 * @date 12 Aug 2015
 */
iftImage* iftEuclideanScoreImage(iftImage* img, float decay);
iftMatrix *iftShapeTransform(iftPoint *orig, int norig, float rx, float rz, iftImage* score);


#endif

#ifdef __cplusplus
}
#endif

