/**
 * @file iftObjectModels.h
 * @brief Definitions and functions about Object Models.
 * @author Samuel Martins
 * @date Jan 01, 2017
 * @ingroup ObjectModels
 */



#ifndef IFT_OBJECTMODELS_H_
#define IFT_OBJECTMODELS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftDict.h"
#include "iftImage.h"
#include "iftParamOptimizationProblems.h"
#include "iftFImage.h"
#include "iftRadiometric.h"
#include "iftRegistration.h"
#include "iftSegmentation.h"
#include "iftSimilarity.h"


/**
 * Function type to extract the Reference Voxel of the Test Image.
 * @param test_img Testing Image.
 * @param params Dictionary of Parameters. Each algorithm must provide (and deallocate after)
 *               its parameters.
 * @return The Reference Voxel of the Testing Image.
 */
typedef iftVoxel (*iftRefVoxelFinderFunc)(const iftImage *test_img, const iftDict *params);


/**
 * @brief Flag of the supported Algorithms on libIFT to find the Reference Voxel on Test Image.
 * @author Samuel Martins
 * @date Dec 16, 2016
 */
typedef enum {
    IFT_IMAGE_CENTER,
    IFT_GEO_CENTER_AFTER_THRESHOLDING
} iftRefVoxelFinderAlg;


/**
 * @brief Algorithms to find the Reference Voxel on Test Image.
 * @author Samuel Martins
 * @date Dec 16, 2016
 */
typedef struct ift_ref_voxel_finder {
    /** Finder algorithm flag */
    iftRefVoxelFinderAlg alg;
    /** Finder function */
    iftRefVoxelFinderFunc func;
    /** Parameters of the finder function */
    iftDict *params;
} iftRefVoxelFinder;


/**
 * @brief Model of a single Object.
 * @author Samuka
 * @date Dec 13, 2016
 * @ingroup ObjectModel
 */
typedef struct ift_obj_model {
    /** Label of the target Object */
    int label;
    /** Model: Prior Probability Map */
    iftFImage *prob_map;
    /**
     * Displacement vector from the Reference Voxel (on test image) to the coordinate where the center voxel of the
     * object model must be positioned (on test image domain).
     * It used to recognize, at the test image, the object to be segmented, positioning
     * the center of the model over it
     */
    iftVector disp_vec;
    // search region for the object localization from the reference voxel
    iftBoundingBox search_region;
} iftObjModel;


/**
 * @brief Multi-Object Model.
 * @author Samuka
 * @date Dec 13, 2016
 * @ingroup ObjectModel
 */
typedef struct ift_m_obj_model {
    /** Depth (in bits) of the images considered by the model */
    int img_depth;
    /** Array of labels of each Object Model */
    iftIntArray *labels;
    /** Array of Object Models. */
    iftObjModel **models;
    /** Pointer to function that computes the Reference Voxel of Images. */
    iftRefVoxelFinder *finder;
    /** Reference Image where the object models were built (if required) and/or for histogram matching when 
    segmenting a testing image. */
    iftImage *ref_img;
    /** Dictionary with specific parameters for each kind of Model: (Statistical, Fuzzy, ...) */
    iftDict *extra;
} iftMObjModel;


/**
 * @brief Struct for the MAGeT-Brain algorithm.
 * @author Samuka
 * @date Jan 13, 2017
 * @ingroup ObjectModels
 */
typedef struct ift_maget {
    /** Set of Template Images (grayscale) */
    iftFileSet *template_set;
    /** Set of Template Label Images. Each template image has a set of label images obtained from the
    registration of an atlas set. */
    iftFileSet **template_label_sets;
    /** (Temp) Working Directory where the templates are stored. It will be removed when destroying this struct. */
    char *workdir;
    /** Depth (in bits) of the images considered by the model */
    int img_depth;
} iftMAGeT;


///////////////////////////////////// FUNCTIONS //////////////////////////////////////////////////
/**
 * @brief Returns the center of the Test Image as the Reference Voxel. It does not required parameters.
 * @param  test_img Testing Image.
 * @param  params   Parameters of the function. Not required (use NULL).
 * @return          The central voxel of the image.
 *
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftVoxel iftRefVoxelByImageCenter(const iftImage *test_img, const iftDict *params);


/**
 * @brief Returns the geometric center of the segmented testing image by threshold (t1 <= test_img[p] <= t2)
 * as the Reference Voxel.
 *
 * Required keys of <b>params</b>: \n
 * "t1": lower thresholding (integer) \n
 * "t2": upper thresholding (integer) \n
 * 
 * @param  test_img Testing Image.
 * @param  params   Parameters of the function. Not required (use NULL).
 * @return          The central voxel of the image.
 *
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftVoxel iftRefVoxelByGeoCenterAfterThres(const iftImage *test_img, const iftDict *params);


/**
 * @brief Creates a finder to determine a Reference Voxel in images. The dict <params> are just assigned (not copied).
 *
 * A finder function is associated for each algorithm's flag <alg>.
 * 
 * @param  alg    Required Algorithm flag.
 * @param  params Parameters used by the finder function.
 * @return        A reference voxel finder.
 *
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftRefVoxelFinder *iftCreateRefVoxelFinder(iftRefVoxelFinderAlg alg, iftDict *params);


/**
 * @brief Destroys a reference voxel finder.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftDestroyRefVoxelFinder(iftRefVoxelFinder **finder);


/**
 * @brief Copies a reference voxel finder.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftRefVoxelFinder *iftCopyRefVoxelFinder(const iftRefVoxelFinder *finder);


/**
 * @brief Destroys a Multi-Object Shape Model.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftDestroyMObjModel(iftMObjModel **mmodel);


/**
 * @brief Reads a Multi-Object Shape Model. It does not read the specific attributes of specialized
 * shape models (sosm, fuzzy, ...).
 * 
 * @param  path    Pathname from the Shape Model.
 * @param  tmp_dir Returns (if != NULL) the pathname from the temp directory where the files were unzipped.
 *                 If NULL, such dir will be deleted.
 * @return         Multi-Object Shape Model.
 *
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftMObjModel *iftReadMObjModel(const char *path, char **tmp_dir);


/**
 * @brief Destroys a Statistical Object Shape Model.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftDestroySOSM(iftMObjModel **sosm);


/**
 * @brief Stores on disk the Statistical Object Shape Model.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftWriteSOSM(const iftMObjModel *sosm, const char *path);


/**
 * @brief Reads from disk a Statistical Object Shape Model.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftMObjModel *iftReadSOSM(const char *path);


/**
 * @brief Destroys a Movable and Adaptive Probabilistic Atlas.
 * @author samuka
 * @date Jan 5, 2017
 * @ingroup ObjectModels
 */
void iftDestroyMAPA(iftMObjModel **mapa);


/**
 * @brief Stores on disk a Movable and Adaptive Probabilistic Atlas.
 * @author samuka
 * @date Jan 5, 2017
 * @ingroup ObjectModels
 */
void iftWriteMAPA(const iftMObjModel *mapa, const char *path);


/**
 * @brief Reads from disk a Movable and Adaptive Probabilistic Atlas.
 * @author samuka
 * @date Jan 5, 2017
 * @ingroup ObjectModels
 */
iftMObjModel *iftReadMAPA(const char *path);


/**
 * @brief Stores on disk the Fuzzy Shape Model.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftWriteFuzzyModel(const iftMObjModel *fuzzy, const char *path);


/**
 * @brief Reads from disk a Fuzzy Shape Model.
 * @author samuka
 * @date Jan 3, 2017
 * @ingroup ObjectModels
 */
iftMObjModel *iftReadFuzzyModel(const char *path);


/**
 * @brief Destroys a Fuzzy Shape Model.
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftDestroyFuzzyModel(iftMObjModel **fuzzy);


/**
 * Writes a MAGeT struture.
 * @author samuka
 * @date Jan 13, 2017
 * @ingroup ObjectModels
 */
void iftWriteMAGeT(const iftMAGeT *maget, const char *path);


/**
 * @brief Reads a MAGeT struct.
 * @author samuka
 * @date Jan 16, 2017
 * @ingroup ObjectModels
 */
iftMAGeT *iftReadMAGeT(const char *path);


/**
 * Destroys a MAGeT struture.
 * @author samuka
 * @date Jan 13, 2017
 * @ingroup ObjectModels
 */
void iftDestroyMAGeT(iftMAGeT **maget);


/**
 * @brief Builds a Statistical Object Shape Model from a set of registered segmentation masks.
 *
 * @param  label_paths Set of pathnames from registered segmentation mask (one for each image).
 * @param  ref_img     Reference Image (coordinate space).
 * @param  labels      Array with the object labels
 * @return             Statistical Object Shape Model.
 *
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftMObjModel *iftBuildSOSM(const iftFileSet *label_paths, const iftImage *ref_img, iftIntArray *labels);



/**
 * @brief Builds a Movable and Adaptive Probabilistic Atlas (MAPA), proposed in [1] (there named SOSM-A).
 *
 * Trained OPF Texture classifier and smooth_factor are assigned in extra dict, respectively, in keys
 * "texture-clf" and "smooth-factor". \n
 * [1] Martins, 2017 - SPIE - A Multi-Object Statistical Atlas Adaptive for Deformable Registration Errors in Anomalous Medical Image Segmentation
 * 
 * @param label_paths   Set of pathnames from registered segmentation mask (one for each image).
 * @param ref_img       Reference Image (coordinate space).
 * @param labels        Array with the object labels
 * @param marker_path   File with the markers used to train the OPF texture classifier.
 * @param smooth_factor Factor (> 0) used to smooth the image before texture classifier (e.g., 0.3).
 * @return              MAPA.
 *
 * @author samuka
 * @date Jan 5, 2017
 * @ingroup ObjectModels
 */
iftMObjModel *iftBuildMAPA(const iftFileSet *label_paths, const iftImage *ref_img, iftIntArray *labels,
                             const char *marker_path, float smooth_factor);



/**
 * @brief Builds a Fuzzy Shape Model from a set of images and their segmentation masks.
 * 
 * @param  img_paths   Set of images.
 * @param  label_paths Segmentation masks of the input images.
 * @param  ref_img     Reference Image used for histogram matching when segmenting a testing image.
 * @param  labels      Array with the object labels.
 * @param  finder      Algorithm to find the reference voxel in images.
 * @param  img_depth   Depth (in bits) of the
 * @return             Fuzzy Shape Model.
 *
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftMObjModel *iftBuildFuzzyModel(const iftFileSet *img_paths, const iftFileSet *label_paths,
                                 const iftImage *ref_img, iftIntArray *labels, iftRefVoxelFinder *finder);



/**
 * @brief Builds a MAGeT-Brain with elastix for registration for the labels <labels>.
 *
 * [1] Pipitone, 2014 - Neuroimage - Multi-atlas segmentation of the whole hippocampus and subfields
 *     using multiple automatically generated templates.
 * 
 * @param  train_img_entry     Dir with the Train. Images or a CSV file with their pathnames.
 * @param  train_label_entry   Dir with the Train. Label Images or a CSV file with their pathnames.
 * @param  template_img_entry  Dir with the Template Images or a CSV file with their pathnames.
 * @param  img_depth           Depth of the Images: e.g: 12 -> [0,4095]
 * @param  labels              Array with the target object labels. 
 * @param  affine_params_path  Elastix Affine (or Rigid or another) Parameter File
 * @param  bspline_params_path Elastix BSpline Parameter File (required for Non-Rigid Registration)
 * @return                     MAGeT-Brain
 *
 * @author samuka
 * @date Jan 17, 2017
 * @ingroup ObjectModels
 */
iftMAGeT *iftBuildElastixMAGeT(const char *train_img_entry, const char *train_label_entry,
                               const char *template_img_entry, int img_depth, const iftIntArray *labels,
                               const char *affine_params_path, const char *bspline_params_path);



/**
 * @brief Finds the Displacement Vector between the Reference Voxel (from each Image)
 * and the Geometric Center of the target object. Histogram of the images are matched with the reference image
 * before.
 *
 * Resulting displacement vector is obtained by the weighted average displacement vector from all images. \n
 * PS: gcs[i] is the geometric center of the object from the image img_paths[i].
 *
 * 
 * @param  img_paths Set of Images.
 * @param  gcs       Array of Geometric Centers of the target object for each image. 
 * @param  finder    Algorithm to find the reference voxel in images.
 * @param  ref_img   Reference Image whose input histograms of the input images are matched.
 * @return           Resulting displacement vector.
 *
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftVector iftFindDispVector(const iftFileSet *img_paths, iftVoxel *gcs, const iftImage *ref_img,
                            const iftRefVoxelFinder *finder);



/**
 * @brief Segments a Testing Image using the Statistical Multi-Object Shape Model SOSM-S (Phellan, 2016).
 *
 * It applies an Object Location by MSPS translating the seed models over the test image gradient. \n
 * Phellan, 2016 - Medical physics - Medical image segmentation via atlases and fuzzy object models
 * 
 * @param  test_img Test Image to be segmented.
 * @param  sosm_s   Statistical Object Shape Model.
 * @return          Segmented Image.
 * 
 * @author samuka
 * @date Jan 3, 2016
 * @ingroup ObjectModels
 */
iftImage *iftSegmentBySOSM_S(const iftImage *test_img, const iftMObjModel *sosm_s, float e_radius, float d_radius);



/**
 * @brief Segments a Testing Image using Movable and Adaptive Probabilistic Atlas (MAPA), proposed in (Martins, 2017) (there called SOSM-A).
 *
 * It applies an Object Location by MSPS translating the seed models over the test image gradient. \n
 * Final Test Gradient is: \n
 * Gt = grad of Test Image \n
 * Gm = gradient of the object map (resulting from test image classification) \n
 * Gi = (alpha * Gt) + ((1-alpha) * Gm) \n
 * Gp = gradient of Prob. Atlas after optimum search \n
 * Gf = (beta * Gi) + ((1-beta) * Gp) \n\n
 * 
 * Martins, 2017 - SPIE - A Atlas Multi-Object Statistical Adaptive for Deformable Registration Errors
 * in Anomalous Medical Image Segmentation.
 * 
 * @param  test_img Test Image to be segmented.
 * @param  mapa     Movable and Adaptive Probabilistic Atlas.
 * @param  grad_alg Gradient Algorithm.
 * @param  alpha    Factor used to combine linearly the Test Image Gradient and the Classified Image Gradient.
 * @param  beta     Factor used to combine linearly the Test Image Gradient and the MAPA gradient.
 * @param  e_radius Erosion Radius that used to find the Inner Seeds (Certain Region Border of the Target Object).
 * @param  d_radius Radius Dilation that is used to find the Outer Seeds (Certain Region Border of the Background).
 * @return          Segmented Image.
 * 
 * @author samuka
 * @date Jan 4, 2016
 * @ingroup ObjectModels
 */
iftImage *iftSegmentByMAPA(const iftImage *test_img, const iftMObjModel *mapa, iftGradientAlg grad_alg,
                           float alpha, float beta, float e_radius, float d_radius);


/**
 * @brief Segments a Testing Image using the Fuzzy Model proposed in (Phellan, 2016).
 *
 * It applies an Object Location by MSPS translating the seed models over the test image gradient. \n
 * Actually, it is a slightly different from [1], because it applies a histogram matching with a 
 * reference image from training set, finding the reference voxel from it. \n\n
 * 
 * Phellan, 2016 - Medical physics - Medical image segmentation via atlases and fuzzy object models
 * 
 * @param  test_img Test Image to be segmented (or the final test image's gradient used for delineation).
 * @param  fuzzy    Fuzzy Model.
 * @return          Segmented Image.
 * 
 * @author samuka
 * @date Jan 3, 2016
 * @ingroup ObjectModels
 */
iftImage *iftSegmentByFuzzyModel(const iftImage *test_img, const iftMObjModel *fuzzy);


/**
 * @brief Segments an Image by Classical MALF. This functions considers the training label set is already registered on testing image.
 *
 * Classical MALF registers all atlases (img + label img) on test image's space. \n
 * The label of each voxel from segmentated image is the most frequent label (Majority Voting) \n
 * 
 * @param  test_img        Testing Image to be segmented.
 * @param  train_label_set Training label set already registered on testing image's space.
 * @return                 Segmented image.
 */
iftImage *iftSegmentByClassicalMALF(const iftImage *test_img, const iftFileSet *train_label_set);


/**
 * @brief Gets the set of registered label images, used to segment a testing image, by MAGeT-Brain
 * (Multiple Automatically Generated Templates) (Pipitone, 2014).
 *
 * It registers each template image (grayscale) on test image's space, and maps the set of template label images
 * after. \n
 * We use Elastix for registration, instead of ANTs as the paper proposes.
 * 
 * Pipitone, 2014 - Neuroimage - Multi-atlas segmentation of the whole hippocampus and subfields
 * using multiple automatically generated templates.
 * 
 * @param  test_img            Testing image which will be segmented.
 * @param  maget               MAGeT approach.
 * @param  affine_params_path  Affine Elaxtix Registation params.
 * @param  bspline_params_path BSpline Elaxtix Registation params.
 * @param  out_dir             Directory where the registered template label images will be stored.
 * @return                     File Set with the registered label images.
 */
iftFileSet *iftGetLabelSetByMAGeT(const iftImage *test_img, const iftMAGeT *maget, const char *affine_params_path,
                                  const char *bspline_params_path, char **tmp_dir);


/**
 * @brief Puts the (cropped) shape model on test image's domain by aligning its center voxel with
 * the reference testing voxel by the displacement vector <disp_vec>.
 * 
 * @param  model          (Cropped) Shape Model
 * @param  test_img       Testing Image.
 * @param  ref_voxel_test Reference Voxel of testing image.
 * @param  disp_vec       Displacement Vector
 * @return                Model on Test Image's domain.
 *
 * @author samuka
 * @date Jan 5, 2017
 * @ingroup ObjectModels
 */
iftFImage *iftModelOnTestImageDomain(const iftFImage *model, const iftImage *test_img,
                                     iftVoxel ref_voxel_test, iftVector disp_vec);


/**
 * @brief Puts the (cropped) shape model image on test image's domain by aligning its center voxel with
 * the reference testing voxel by the displacement vector <disp_vec>.
 * 
 * @param  model_img      (Cropped) Shape Model Image
 * @param  test_img       Testing Image.
 * @param  ref_voxel_test Reference Voxel of testing image.
 * @param  disp_vec       Displacement Vector
 * @return                Model Image on Test Image's domain.
 *
 * @author samuka
 * @date Jan 5, 2017
 * @ingroup ObjectModels
 */
iftImage *iftModelImageOnTestImageDomain(const iftImage *model_img, const iftImage *test_img,
                                          iftVoxel ref_voxel_test, iftVector disp_vec);


/**
 * Computes the gradient of a Shape Model
 * 
 * @param  model       Shape model.
 * @param  max_img_val Maximum value of the images, based on image depth: (255, 4095, ...)
 * @param  grad_alg    Gradient Algorithm.
 * @return             Gradient of the Shape Model.
 * 
 * @author samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftImage *iftShapeModelGradient(const iftFImage *model, int max_img_val, iftGradientAlg grad_alg);


/**
 * @brief Finds the Object Model's Seeds (on test image space) for delineation.
 *
 * Given a Model (probabilistic map), background's seeds are those with prob 0 and object's seeds
 * with prob 1, after eroding and dilating the target object. \n
 * Finally, resulting seeds are filtered by a membership map (if != NULL), which is the test image previously classified,
 * where 0 are the background and 1 is the target object. \n
 * Filtering follows the rule: Inner seeds classified as backround become outer seeds.
 * 
 * @param  test_img       Testing Image.
 * @param  ref_voxel_test Reference Voxel from the test image.
 * @param  obj_model      Object Model.
 * @param  e_radius       Erosion radius used to get the inner seeds for delineation.
 * @param  d_radius       Dilation radius used to get the external seeds for delineation.
 * @param  membership_map Membership Map (test image previsously classified) used to filtered the seeds. If NULL, no filtering.
 * @return                Labeled seeds for delineation.
 *
 * @author samuka
 * @date Jan 3, 2017
 * @ingroup ObjectModels
 */
iftLabeledSet *iftFindModelSeeds(const iftImage *test_img, iftVoxel ref_voxel_test, const iftObjModel *obj_model,
                                 double e_radius, double d_radius, const iftImage *membership_map);



/**
 * @brief Register a test image on reference image's space of a multi-object shape model.
 *
 * Deformation fields are not saved.
 * 
 * @param test_img            Testing Image.
 * @param mmodel              Multi-Object Statistical Shape Model.
 * @param affine_params_path  Affine Elaxtix Registation params.
 * @param bspline_params_path BSpline Elaxtix Registation params.
 * @return                    Registered Testing Image Space
 *
 * @author samuka
 * @date Jan 6, 2017
 * @ingroup ObjectModels
 */
iftImage *iftRegisterImageOnMObjModelByElastix(const iftImage *test_img, const iftMObjModel *mmodel,
                                               const char *affine_params_path, const char *bspline_params_path);


/**
 * @brief Register (by Elastix) a Multi-Object Statistical Shape Model on a test image's space. Resulting registrations
 * and mapping are assigned in the own input model.
 * 
 * Reference image is registerd with the test one. This will be the new reference image. \n
 * Then, all object models are mapped to new space using the deformation fields. \n
 * Deformation fields are not saved.
 * 
 * @param mmodel              Multi-Object Statistical Shape Model.
 * @param test_img            Testing Image.
 * @param affine_params_path  Affine Elaxtix Registation params.
 * @param bspline_params_path BSpline Elaxtix Registation params.
 *
 * @author samuka
 * @date Jan 6, 2017
 * @ingroup ObjectModels
 */
void iftRegisterMObjModelByElastix(iftMObjModel *mmodel, const iftImage *test_img, const char *affine_params_path,
                                   const char *bspline_params_path);


/**
 * @brief      Selects the best atlases to a test image by Normalized Mutual Information. Train set and test
 *             image must be already registered.
 * 
 * The higher the NMI, more similar (better) the train image is. \n
 * If a train label set is passed, it considers only the bounding box around the objects to compute the NMI. \n
 * Be sure that train_label_set[i] corresponds to label image from image train_img_set[i].
 *
 * @param[in]  test_img        (Registered) Testing Image.
 * @param[in]  train_img_set   (Registered) Training Image Set.
 * @param[in]  train_label_set (Registered) Training Label Set. If != NULL, it considers tmasks the train_img_set before NMI.
 * @param[in]  n_atlases       Number of atlases to be selected.
 * @param      chosen_idxs_out If != NULL, returns the indexes of the selected atlas in the train image set.
 *
 * @return     A file set with the pathnames from the selected images.
 * 
 * @author Samuka
 * @date Jan 27, 2017
 * @ingroup ObjectModels
 */
iftFileSet *iftAtlasSelectionByNMI(const iftImage *test_img, const iftFileSet *train_img_set,
                                   const iftFileSet *train_label_set, int n_atlases, iftIntArray **chosen_idxs);


/**
 * @brief      Runs the program iftSegmentByClassicalMALF and returns the segmented Image.
 *
 * @param[in]  test_img_path        Test Image pathname.
 * @param[in]  img_entry            Train Image Set entry.
 * @param[in]  label_entry          Train Label Set entry.
 * @param[in]  affine_params_path   Affine parameter file (it can be NULL).
 * @param[in]  bspline_params_path  BSpline parameter file (it can be NULL).
 *
 * @return     Segmented Image By MALF.
 * 
 * @author Samuka
 * @date Jan 30, 2017
 * @ingroup ObjectModels
 */
iftImage *iftRunProgSegmentByClassicalMALF(const char *test_img_path, const char *img_entry,
                                           const char *label_entry, const char *affine_params_path,
                                           const char *bspline_params_path);



/**
 * @brief Apply the N4 algorithm for Inhomogeneity Correct in a MRI Image.
 *
 * It is used the default parameters, suggested by 3D Slicer tool [2], for N4 [1].\n
 * [1] Tustison, Nicholas J., et al. \"N4ITK: improved N3 bias correction.\" IEEE transactions on medical imaging 29.6 (2010): 1310-1320.\n
 * [2] https://www.slicer.org/
 * 
 * @param  input_img Image to be corrected.
 * @param  mask      Binary mask that defines the structure of your interest. If NULL, N4 considers the entire image.
 * @param  out_bias  The resulting Bias Field Image by the N4 correction. If NULL, it is not considered. 
 * @return           Correct image by N4.
 */
iftImage *iftN4BiasFieldCorrection(const iftImage *input_img, const iftImage *mask, iftImage **out_bias);



#ifdef __cplusplus
}
#endif

#endif


