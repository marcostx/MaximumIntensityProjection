#ifndef IFT_DATASET_H_
#define IFT_DATASET_H_

/**
 * @file iftDataSet.h
 * @brief Definitions and functions for Machine Learning datasets.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftLabeledSet.h"
#include "iftFImage.h"
#include "iftRadiometric.h"
#include "iftMatrix.h"
#include "iftSort.h"
#include "iftSet.h"
#include "iftKernel.h"
#include "iftParseInput.h"
#include "iftDescriptors.h"

typedef enum ift_sample_status {
    UNKNOWN = 0x01,
    TEST = 0x02,
    TRAIN = 0x04,
    OUTLIER = 0x08,
    ERROR = 0x10,
    ARTIFICIAL = 0x20,
    ALL = TEST | TRAIN | OUTLIER | ERROR | UNKNOWN | ARTIFICIAL
} iftSampleStatus;

#define WEIGHT    0
#define LABEL     1
#define CLASS     2
#define POINT     3
#define STATUS    4

#define K_distance5 20.0

/**
 * @brief Precompute (in a discrete approximated way) the Minkowski distance between samples.
 * @author Peixinho (but it is not mine, I just docummented this, I swear)
 * @warning Requires the data to be norm one. See also iftNormOneDataSet().
 * @date May, 2016
 */
typedef struct ift_minkowski_table {
    float *dist;    /* look-up table for Minkowski distance computation */
    int nelems;  /* number of elements: 2*mult_factor + 1 */
    int mult_factor; /* multiplication factor: indicates the resolution of the precomputed distances approximation*/
    float exp_factor;  /* exponential factor */
} iftMinkowskiTable;

/**
 * @brief Global variable to store the precomputed approximated minkowski function.
 * Used in the iftDistance7()
 * @warning Requires the samples to have norm one. See also iftNormOneDataSet()
 * @author Peixinho
 * @date May, 2016
 */
extern iftMinkowskiTable *iftMinkowski; /* Use iftSetMinkowskiTable to initialize it */

/**
 * @brief Distance function between two arrays prototype.
 * The alpha coefficients alpha can be used to modify the arc weight
computation in many different ways. These coefficients can be found
by optimization.
 * @date may, 2016
 * @author Peixinho
 */
typedef float (*iftArcWeightFun)(float *f1, float *f2, float *alpha, int n);


/**
 * @brief Precomputed distances between dataset samples.
 * @author Peixinho
 * @date May, 2016
 */
typedef struct ift_distance_table {
    float **distance_table; /* n x n table of distances, for n samples */
    int nsamples1;
    int nsamples2;
} iftDistanceTable;

/**
 * @brief Global variable to store the precomputed distances for a dataset.
 * @author Peixinho
 * @date May, 2016
 */
extern iftDistanceTable *iftDist; /* Use iftSetDistanceFunction
					  to change it to one of the
					  iftDistanceX, X=1,2...,7 */

/**
 * @brief Sampling strategy to split a dataset.
 * Splits dataset in TRAIN and TEST for several runs.
 * @author Peixinho
 */
typedef struct ift_sampler {
    char **status;
    int niters;
    int nsamples;
} iftSampler;

/**
 * @brief A sample in a Pattern Recognition problem.
 * Contains the feature descriptor and, for supervised cases, the sample class.
 * @date May, 2016
 * @author Peixinho
 */
typedef struct ift_sample {
    float *feat;    // feature values
    int truelabel;   // 1,2,...,nclasses
    int label;   // 1,2,....,nlabels
    int id;      // identification which may be either a voxel of the
    int cluster_id; // 1,2,...,nClusters (used for OPF_unsup)
    // related image or the position of an image in the
    // directory with a related list of images.
    bool isSupervised;

    float weight;  // probability density, path cost, etc
    uchar status;  // Status of the sample: TEST, TRAIN, ERROR, OUTLIER.

} iftSample;

/**
 * @brief Structure used for feature space transformations based on PCA
 * @author Peixinho
 * @date May, 2016
 * */
typedef struct _featspaceparam {
    float *mean, *stdev; // parameters used for feature space normalization and centralization
    int nfeats; // size of mean, stdev
    char *w; // binary vector that indicates selected components for SupPCA
    int ncomps; // size of w (number of PCA components before
    // feature space reduction by SupPCA)
    iftMatrix *R; // Projection/rotation matrix (the columns are the eigenvectors)
    iftMatrix *W; // Whitening Transformation
} iftFeatSpaceParam;


typedef enum {
  IFT_REF_DATA_NULL, IFT_REF_DATA_IMAGE, IFT_REF_DATA_FIMAGE, IFT_REF_DATA_MIMAGE, IFT_REF_DATA_FILESET, IFT_REF_DATA_MMKERNEL
} iftRefDataType;

/**
 * @brief Dataset structure for pattern recognition.
 * @author Peixinho
 * @date May, 2016
 * Contains a list of samples. Each sample containing the description features. And, for supervised cases, the sample class.
 * See also ::iftSample.
 * @ingroup MachineLearning
 */
typedef struct ift_dataset {
    iftSample *sample;   // sample
    int capacity; //max number of samples (pre allocates this space)
    int nfeats;   // number of features
    int nsamples; // number of samples
    int nlabels;  // number of clusters
    int nclasses; // number of ground truth classes
    int ntrainsamples; // number of training samples
    iftArcWeightFun iftArcWeight; // distance function
    int function_number; // This field is used to store the ID of the distance function set with iftSetDistanceFunction. It allows the
    // proper setting of the distance function when using iftWriteOPFDataset and iftReadOPFDataset
    float *alpha;   // coefficients used for arc weight computation
    iftRefDataType ref_data_type; // type of the Reference Data (mainly used to read and write datasets)
    void *ref_data; // it migh be an image, for voxel datasets, a text file with image filenames, for image datasets, a region graph, for supervoxel datasets, etc.
    iftFeatSpaceParam fsp; // parameters of the feature scape transformation
    iftMatrix* data;
    iftDoubleMatrix* projection;
} iftDataSet;

/**
 * @brief Creates a Sampler Strategy object
 * @param nsamples Number of samples.
 * @param niters Number of iterations.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 */
iftSampler *iftCreateSampler(int nsamples, int niters);

/**
 * @brief Destroys a Sampler Strategy object.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 */
void iftDestroySampler(iftSampler **sampler);

/**
 * @brief Implements a KFold Sampling method.
  * @param nsamples Number of samples.
  * @param k Fold size.
  * @author Peixinho
  * @date May, 2016
  * @ingroup MachineLearning
 */
iftSampler *iftKFold(int nsamples, int k);

/**
 * @brief Repeats N times a KFold Sampling method. The default usage is a 5x2 cross validation.
 * @param nsamples Number of samples.
 * @param n Number of times to run kfold.
 * @param k Fold size.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 */
iftSampler *iftNKFold(int nsamples, int n, int k);

/**
 * @brief Leave One Out Sampling method. Uses all data to train, except one sample used to test. Repeats the process to each sample.
 * @param nsamples Number of samples.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 */
iftSampler *iftLeaveOneOut(int nsamples);

/**
 * @brief Repeats N times a Random Subsampling method.
 * @param nsamples Number of samples.
 * @param n Number of times to run subsampling.
 * @param ntrain Number of training samples.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 */
iftSampler *iftRandomSubsampling(int nsamples, int n, int ntrain);


/**
 * @brief Repeats N times a Random Subsampling method, keeping the classes a priori distribution.
 * @param nsamples Number of samples.
 * @param n Number of times to run subsampling.
 * @param ntrain Number of training samples.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 * @return The sampler.
 */
iftSampler *iftStratifiedRandomSubsampling(const int *labels, int n_samples, int n_iters, int n_train);


/**
 * @brief Repeats N times a Balanced Subsampling method.
 *
 * Balances the classes distribution by setting all classes with the same number of samples X that the smaller one has. \n
 * After that, it applies a random sampling by selecting X samples for each classes, resulting in a sampled \n
 * dataset with C*X samples, where C is the number of classes. \n
 * From this sampled dataset, it applies a new sampling in order to select the number of training samples <b>n_train</b>. \n
 * If <b>n_train</b> is greater than C*X, all sampled dataset (C*X samples) will be used. \n
 * Repeat all steps for each iteration.
 * 
 * @param  labels    Labels from the samples.
 * @param  n_samples Number of Samples.
 * @param  n_iters   Number of iterations.
 * @param  n_train   Number of train. samples.
 * @return           The resulting sampler.
 * 
 * @author Samuel Martins
 * @date Jun 22, 2016
 * @ingroup MachineLearning
 */
iftSampler *iftBalancedRandomSubsampling(const int *labels, int n_samples, int n_iters, int n_train);


/**
 * @brief Selects TRAIN and TEST samples according to the sampler strategy.
 * @param Z Dataset to split.
 * @param sampler Sampler strategy
 * @param iteration Current sampling iteration.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 */
void iftDataSetSampling(iftDataSet *Z, const iftSampler *sampler, int iteration);

/**
 * @brief Get true label array (one true label for each sample) from DataSet samples.
 * @author Peixinho
 * @date May, 2016
 * @ingroup MachineLearning
 *
 * @param dataset Dataset.
 * @return The array of true labels.
 */
iftIntArray *iftGetDataSetTrueLabels(const iftDataSet *dataset);

/**
 * @brief Check if the dataset contains class information and therefore is able to perform supervised classification.
 * @param Z The dataset to be checked
 * @return True if is a supervised dataset, false otherwise.
 * @date May, 2016
 * @author Peixinho
 * @ingroup MachineLearning
 */
bool iftIsSupervisedDataSet(const iftDataSet *Z);

/**
 * @brief Check if the dataset is built of training samples.
 * @param Z The dataset to be checked
 * @return True if is a training dataset, false otherwise.
 * @date May, 2016
 * @author Peixinho
 * @ingroup MachineLearning
 */
bool iftIsTrainingDataSet(const iftDataSet *Z);

/**
 * @brief Check if the dataset is built of testing samples.
 * @param Z The dataset to be checked.
 * @return true if is a testing dataset, false otherwise.
 * @date May, 2016
 * @author Peixinho
 * @ingroup MachineLearning
 */
bool iftIsTestingDataSet(const iftDataSet *Z);

/**
 * @brief Check if the dataset has normalized data. See also iftNormalizeDataSet() and iftNormalizeDataSet2().
 * @param Z The dataset to be checked
 * @date May, 2016
 * @author Peixinho
 * @return true if is a normalized dataset, false otherwise.
 */
bool iftIsNormalizedDataSet(const iftDataSet *Z);

/**
 * @brief Check if the dataset has centralized data. See also iftCentralizeDataSet().
 * @param Z The dataset to be checked
 * @return true if is a centralized dataset, false otherwise.
 * @ingroup Peixinho
 */
bool iftIsCentralizedDataSet(const iftDataSet *Z);

/**
 * @brief Check if the dataset has whitened data. See also iftWhiteningTransform.
 * @param Z The dataset to be checked
 * @return true if is a whitened dataset, false otherwise.
 */
bool iftIsWhitenedDataSet(const iftDataSet *Z);

/**
 * @brief Check if the dataset is preprocessed by PCA. See also iftTransFeatSpaceByPCA(), iftTransFeatSpaceByPCA2() and iftTransFeatSpaceByPCA3().
 * @param Z The dataset to be checked
 * @return true if is a whitened dataset, false otherwise.
 */
char iftIsTransformedDataSetByPCA(iftDataSet *Z);

/**
 * @brief Creates an empty dataset.
 * @param nsamples Number of samples in the dataset.
 * @param nfeats Number of features in dataset samples.
 * @date May, 2016
 * @return The new dataset.
 */
iftDataSet *iftCreateDataSet(int nsamples, int nfeats);
iftDataSet* iftCreateDataSetAux(int nsamples, int nfeats,bool allocData);


/**
 * @brief append samples in the dataset
 * @param dataset
 * @param sample
 * @author Peixinho
 * @date Jul, 2017
 */
void iftAddSampleDataset(iftDataSet* dataset, iftSample sample);

/**
 * @brief Expands dataset capacity, to allow add more samples
 * @param dataset
 * @param capacity
 * @date Jul, 2017
 * @author Peixinho
 */
void iftGrowDataset(iftDataSet* dataset, int capacity);



/**
 * @brief Creates an empty dataset, without feature vectors allocation. See also iftCreateDataSet().
 * Create an iftDataSet with no Feature allocation, used for address copy between datasets.
 * @param nsamples Number of samples in the dataset.
 * @param nfeats Number of features in dataset samples.
 * @date May, 2016
 * @return The new dataset.
 */
iftDataSet *iftCreateDataSetWithoutFeatsVector(int nsamples, int nfeats);


/**
 * @brief Creates/Allocates an array of Data Set Pointers (without allocated them).
 * @author Samuel Martins
 * @date Jun 8, 2016
 * 
 * @param  n_datasets Number of Data Sets
 * @return            Allocated Array of Data Sets.
 */
iftDataSet **iftCreateDataSetArray(size_t n_datasets);


/**
 * @brief Destroys a dataset object (and destroy their Sample Feature Vectors).
 * @author Peixinho
 * @date May, 2016
 * @param Z Dataset to be destroyed.
 */
void iftDestroyDataSet(iftDataSet **Z);


/**
 * @brief Destroys a dataset object without deallocating their Sample Feature Vectors.
 * @author Samuel Martins
 * @date Jun 6, 2016
 * @param Z Dataset to be destroyed.
 */
void iftDestroyDataSetWithoutFeatVectors(iftDataSet **Z);


/**
 * @brief Destroys all datasets from an array (deallocating or not their sample feat. vectors) and
 * destroys the array.
 * @author Samuel Martins
 * @date Jun 8, 2016
 * 
 * @param Z_array       Array of Data Sets
 * @param n_datasets    Number of Datasets
 * @param destroy_feats If true, it destroys the feature vector from the features of all datasets.
 */
void iftDestroyDataSetArray(iftDataSet ***Z_array, size_t n_datasets, bool destroy_feats);


/**
 * @brief Creates a copy of the dataset (copying or just pointing their feature vectors).
 * @date May, 2016
 *
 * @warning If <b>copy_feats</b> is false, the feature vectors are only assigned/pointed to the original ones. \n
 * Then, BE CAREFUL when destroying this dataset. \n
 * Use the function @ref iftDestroyDataSetWithNoFeatVectors()
 * 
 * @param Z          The dataset to be copied.
 * @param copy_feats If true, it copies the sample feature vectors. Otherwise, it only points them.
 * @return The new dataset.
 */
iftDataSet *iftCopyDataSet(const iftDataSet *Z, bool copy_feats);


/**
 * @brief Validates if a group of datasets has the same parameters (Number of samples, Number of features and Number of classes).
* @param Zs Array of datasets to be validated.
* @param num_Z Number of datasets.
* @date May, 2016
 * @author Peixinho
*/
void iftValidateDataSets(iftDataSet **Zs, int num_Z);

/**
 * @brief Computes the Object Map from Z for the predicted label <b>label</b>.
 * @authors Alexandre Falcao, Samuel Martins
 * @date Jul 2, 2016
 * 
 * @param  Z       Input Dataset.
 * @param  max_val Max value for the resulting Image (ex: 255, 4095, ...)
 * @param  label   Target predicted label.
 * @return         Object Map for the target predicted label.
 */
iftImage *iftDataSetObjectMap(const iftDataSet *Z, int max_val, int label);


/**
 * @brief Creates a dataset from an image, where each sample has the voxel brightness value (gray image), or the Y, Cb an Cr (color image).
 * @param img Image to be converted.
 * @date May, 2016
 * @return The image dataset.
 */
iftDataSet *iftImageToDataSet(iftImage *img);


/**
 * @brief Try to breaks a multidimensional image in <b>n_tiles</b> tiles, and creates for each tile an independent dataset
 * @param mimg A dataset which referenced data contains a multidimensional image
 * @param n_tiles Number of tiles
 * @param out_datasets A new array of datasets
 * @param copy_features True to copy the features vectors and false for not copy them
 * @return Number of datasets created
 */
int iftMImageTilesToDataSetArray(iftDataSet *orig_dataset, int n_tiles, iftDataSet *** out_datasets, bool copy_features);

/**
 * @brief Creates a dataset with the elements contained in a bounding box <b>bb</b> of an image <b>img</b>
 * @param img An image
 * @param bb A bounding box
 * @return A new dataset
 */
iftDataSet *iftImageBoundingBoxToDataSet(iftImage *img, iftBoundingBox bb );

/**
 * @brief Creates a dataset with the elements contained in a bounding box <b>bb</b> of a multidimensional image
 * <b>img</b>
 * @param iftDataset A dataset containing a multidimensional image
 * @param bb A bounding box
 * @param copy_feats True to copy the features vectors and false for not copy them
 * @return A new dataset
 */
iftDataSet *iftMImageBoundingBoxToDataSet(iftDataSet *orig_dataset, iftBoundingBox bb, bool copy_feats);

/**
 * @brief Split a big dataset into an array of <b>nb_partitions</b> datasets where in each one the samples are picked randomly without repetition.
 * @param original_dataset A big dataset
 * @param nb_partitions Number of partitions to divide the dataset
 * @param out_datasets Out parameter that will contain the dataset array
 * @param copy_features True to copy the features vectors and false for not copy them
 * @return The number of datasets created, this value is <b>ceil(original_dataset->nelements/n_elem_for_dataset)</br>
 */
int iftSplitDatasetIntoRandomlyDataSetArray(iftDataSet *original_dataset, int nb_partitions,
                                            iftDataSet ***out_datasets, bool copy_features);

/**
 * @brief Copy the ground truth image label to the dataset.
 * @author peixinho
 * @note Adan extended this function to support ground truth images with more than two labels
 * @date May, 2016
 * @param imgGT GT image with labels information.
 * @param Z Output dataset.
 */
void iftImageGTToDataSet(iftImage *imgGT, iftDataSet *Zpixels);

/**
 * @brief Creates a dataset composed by the image, considering an adjacency A.
 * The output dataset contains as many samples as the number of pixels/voxels, and each sample feature is composed by the the pixel/voxel information (Gray/YCbCr) along the specified adjacency.
 * @param img Input image
 * @param A Adjacency to be considered as sample features.
 * @return The image dataset
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftImageToDataSetUsingAdjacency(iftImage *img, iftAdjRel *A);

/**
 * @brief Creates a dataset composed by the pixel/voxel locations.
 * The output dataset contains as many samples as the number of pixels/voxels, and each sample feature is composed by the the pixel location (xy/xyz) in the image domain.
 * @param label
 * @return
 */
iftDataSet *iftObjectToDataSet(iftImage *label);


iftDataSet *iftReadXYDataSet(char *filename);

/**
 * @brief Creates a dataset with seed samples from an image.
 * @param img Input image
 * @param S Seeds to be copied to the dataset.
 * @return The dataset
 * @sa iftMImageSeedsToDataSet()
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftImageSeedsToDataSet(iftImage *img, iftLabeledSet *S);

/**
 * @brief Creates a dataset with seed samples from a multiband image.
 * @param img Input image
 * @param S Seeds to be copied to the dataset.
 * @return The dataset
 * @sa iftImageSeedsToDataSet
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftMImageSeedsToDataSet(iftMImage *img, const iftLabeledSet *S);


/**
 * @brief Creates a dataset from a specified image region, where each sample has the voxel/pixel value (gray/YCbCr), accor
 * @param img Input image
 * @param region Image defining the region of interest. All non zero values are considered.
 * @return The dataset
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftImageRegionToDataSet(iftImage *img, iftImage *region);

/**
 * @brief Creates a dataset from an image. Each sample is the average color of the region.
 * @param img Input image
 * @param supervoxel Supervoxel information
 * @return The dataset
 */
iftDataSet *iftSupervoxelsToDataSet(const iftImage *img, const iftImage *supervoxel,
                                    const iftIntArray *true_labels);


/**
 * @brief Creates a dataset from a GRAY Image where each sample comes from a superpixel/supervoxel and has
 * the features: [average value, stardard deviation, minimum value, maximum value]. 
 * 
 * If the int array <b>true_labels</b> is passed (!= NULL), each supervoxel sample on dataset will have
 * the true label assigned for the array. \n
 * The same indexing scheme holds for it: ex: the true label for the supervoxel 2 will be in true_labels[1]. \n\n
 * If it is not passed (= NULL), the true label of the Supervoxels will be the own supervoxel labels.
 * 
 * @author Samuel Martins
 * @date Jun 1, 2016
 * 
 * @param  img           Original Gray Image where the feature vectors will be extracted.
 * @param  super_img     Supervoxel Image that determines the coordinates from each supervoxel.
 * @param  max_range     Maximum range value of the image <b>img</b> used for feature normalization.
 *                       Ex: For 12 bits: it is 4095 (2^12 - 1).
 *                       If its value is <= 0, this function tries to find it automatically.
 * @param  true_labels   Array of the True Labels of the Supervoxels.
 *                       If != NULL, they are assigned to the dataset samples. 
 * @return               Data Set with the feature vector for each supervoxel.
 */
iftDataSet *iftSupervoxelsToGrayAvgStdevMinMaxDataSet(const iftImage *img, const iftImage *super_img, int max_range,
                                                      bool normalize, const iftIntArray *true_labels);


iftDataSet *iftSupervoxelsToGrayValsVoxelCoordsDataSetInRegion(const iftImage *img, const iftImage *super_img,
                        iftBoundingBox bb, int max_range, bool normalize, iftIntArray *true_labels);


iftDataSet *iftSupervoxelsToGrayValsVoxelCoordsDataSet(const iftImage *img, const iftImage *super_img, int max_range,
                                                      bool normalize, iftIntArray *true_labels);

iftDataSet *iftSupervoxelsToVoxelCoordsDataSetInRegion(const iftImage *super_img, iftBoundingBox bb,
                                bool normalize, iftVoxel *centroids_in, iftIntArray *true_labels);

/** Builds a Dataset where each supervoxel is a sample and its feature vector is the voxel coordinates
from its cluster/supervoxel centroid (root)

centroids = indices (0..n-1) in the super_img from the centroids (roots). If NULL, the centroids
will be computed.
*/
iftDataSet *iftSupervoxelsToVoxelCoordsDataSet(const iftImage *super_img, bool normalize,
                                        iftVoxel *centroids_in, iftIntArray *true_labels);

/**
 * @brief Creates a dataset from an image. Each sample is composed by the region histogram.
 * @param img Input image
 * @param supervoxel Supervoxel information
 * @param nbins Number of bins in the histogram
 * @return The dataset
 */
iftDataSet *iftSupervoxelsToHistogramDataSet(const iftImage *img, const iftImage *supervoxel, int nbins, iftIntArray *true_labels);


/**
 * @brief Creates a dataset from a supervoxel image. Each sample is composed by the Lab histogram of a corresponding superpixel. The histograms are not normalized. The dataset if set with the distance 12 (Bhattacharyya coefficient for discrete probability distributions).
 * @param img Input image
 * @param label A label image
 * @param nbins Number of bins to compose each histogram (preferred a cubic number)
 * @params region_size Out parameter to save the count of pixels that has each superpixel
 * @return The dataset
 */
iftDataSet *iftSupervoxelsToLabOrGrayHistogramDataSet(const iftImage *img, const iftImage *label, int nbins, int *region_size);

/**
 * @brief Creates a dataset from an image. Each sample is composed by the sum of color bands, the proportional area of the image and the number of pixels.
 * @param image Input image
 * @param supervoxel Supervoxel information
 * @param colorspace
 * @return The dataset
 */

iftDataSet *iftSupervoxelsToMeanSizeDataSet(iftImage *image, iftImage *supervoxel, iftColorSpace colorspace);

/**
 * @brief Creates a dataset from a supervoxel image. Features in the detailed description.
 * Features: Band1 sum, Band2 sum, Band3 sum,
 * relative area, npixels,
 * xmax, ymax, xmin, ymin,
 * Band1 histogram, Band2 histogram, Band3 histogram and
 * Local Binary Pattern
 * @param image The input image
 * @param supervoxels Supervoxel information
 * @param bins_per_band Number of bins for histograms
 * @param colorspace The considered colorspace
 * @return the dataset
 *
 */

iftDataSet *iftSupervoxelsToSelectiveSearchDataset(iftImage *image, iftImage *supervoxels, int bins_per_band, iftColorSpace colorspace);

/**
 * @brief Creates a dataset from a supervoxel image. Features in the
 * detailed description.
 * Features in Lab space (9):
 * (normalized) mean L, mean a, mean b,
 * std L, std a, std b,
 * skewness L, skewness a, skewness b
 * @param image The input image
 * @param supervoxels Supervoxel information
 * @param bins_per_band Number of bins for histograms
 * @param colorspace The considered colorspace
 * @return The dataset
 *
 */

iftDataSet *iftSupervoxelsToLabColorMomentsDataset(iftImage *image, iftImage *supervoxels);

iftDataSet *iftSupervoxelsToLabHistogramDataset(iftImage *image, iftImage *label_image, int bins_per_band);

/**
 * @brief Extracts a Feature Vector using BIC for each superpixel from an image and returns them in a Data set.
 *
 * For each supervoxel, a BIC feature vector will be extracted from <b>img</b>, being stored into a data set. \n
 * The Supervoxel Image MUST HAVE ONLY labels from 1 to n. Zero labels will be ignored. \n
 * Thus, for example, the supervoxel <b>i</b> will be the sample <b>i-1</b> in the data set. \n\n
 *
 * If the int array <b>true_labels</b> is passed (!= NULL), each supervoxel sample on dataset will have
 * the true label assigned for the array. \n
 * The same indexing scheme holds for it: ex: the true label for the supervoxel 2 will be in true_labels[1]. \n\n
 * If it is not passed (= NULL), the true label of the Supervoxels will be the own supervoxel labels.
 * 
 * 
 * @authors Somebody, Samuel Martins
 * @date Jun 1, 2016
 * 
 * @param  img           Original Image where the BIC feature vectors will be extracted.
 * @param  super_img     Supervoxel Image that determines the coordinates from each supervoxel.
 * @param  bins_per_band Number of bins for the BIC histograms.
 * @param  max_range     Maximum range value of the image <b>img</b>. Ex: For 12 bits: it is 4095 (2^12 - 1).
 *                       If its value is <= 0, this function tries to find it automatically.
 * @param  true_labels   Array of the True Labels of the Supervoxels.
 *                       If != NULL, they are assigned to the dataset samples. 
 * @return               Data Set with the feature vector for each supervoxel.
 */
iftDataSet *iftSupervoxelsToBICDataset(iftImage *img, iftImage *super_img, int bins_per_band,
                                       int max_range, iftIntArray *true_labels);

iftDataSet *iftSupervoxelsToUniformLBP2D(iftImage *image, iftImage *label_image);

iftDataSet *iftSupervoxelsToSimplifiedHOG2D(iftImage *image, iftImage *label_image, int nbins);

iftDataSet *iftConcatDatasetFeatures(iftDataSet **datasets, int ndatasets);

iftDataSet *iftRegionMergingDataSet(iftMImage *image, iftImage *label_image);

iftDataSet *iftMSupervoxelsToDataSet(iftMImage *mimg, iftImage *label);

iftDataSet *iftImageBorderDataSet(iftDataSet *Z1, int size, int nsamples);

/**
 * Creates a binary matrix (nsamples x nclasses) in a one hot encoding of labels
 * @param Z The dataset
 * @param status Status of the considered samples
 * @return The labels binary matrix
 */
iftMatrix* iftDataSetToLabelsMatrix(iftDataSet *Z, char status);

/**
 * @brief Read a CSV formatted dataset
 * @author Peixinho
 * @date Nov, 2016
 * @param filepath Path to the CSV file
 * @param separator Character separator (usually ',' or ' ').
 * @param isSupervised Indicates if the dataset has a truelabel information
 * @param truelabelColumn The column containing the truelabel information, supports negative index starting from the last column.
 * @return Loaded dataset.
 */
iftDataSet* iftReadCSVDataSet(const char* filepath, char separator, bool isSupervised, int truelabelColumn, bool hasId, int idColumn);

iftDataSet* iftReadCSVImageDataSet(const char* filepath, char separator,
                                   bool isSupervised, int truelabelColumn,
                                   const char* imagesPath);
/**
 * Creates a binary matrix (nsamples x nfeats+1) including a bias column of ones in the feature vector.
 * @param Z The dataset
 * @param status Status of the considered samples
 * @return The features matrix
 */
iftMatrix*iftDataSetToFeatureMatrixHomogCoord(iftDataSet *Z, char status);

/**
 * @brief Linear Discriminant Analisys for supervised dimensionality reduction.
 * Reduces the feature space to a R^(c-1) space, where c is the number of classes.
 * Following https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/discriminant_analysis.py
 * @author Peixinho
 * @date Jul, 2016
 * @return The basis transformation matrix. See also iftDataSetTransformBasis()
 */
iftMatrix* iftLinearDiscriminantAnalysis(const iftDataSet *Z);


void iftWriteCSVDataSet(iftDataSet* dataset, const char* filename);

/**
 * @brief Reads an OPF Dataset.
 * @author Samuel Martins
 * @date Jun 28, 2016
 * 
 * @param  pathname Dataset's pathname.
 * @return          Loaded Dataset.
 */
iftDataSet *iftReadOPFDataSet(const char *pathname);


/**
 * @brief Writes an OPF Dataset.
 * @author Samuel Martins
 * @date Jun 28, 2016
 * 
 * @param Z        Dataset to be stored.
 * @param pathname Dataset's pathname.
 */
void iftWriteOPFDataSet(const iftDataSet *Z, const char *pathname);


/**
 * @brief Reads the feature space parameter information.
 * @author Samuel Martins
 * @date Jun 28, 2016
 *
 * @param path File path to save the feature space param.
 * @return     Feature space param.
 */
iftFeatSpaceParam iftReadFeatSpaceParam(const char *path);

/**
 * @brief Stores the feature space parameter information.
 * @author Samuel Martins
 * @date Jun 28, 2016
 *
 * @param fsp  Reference of the Feature space param will be stored.
 * @param path File path to save the feature space param.
 */
void iftWriteFeatSpaceParam(const iftFeatSpaceParam *fsp, const char *path);


/**
 * @brief Remove information added to the dataset (labels and train/test status).
 * @param Z Dataset to reset.
 */
void iftResetDataSet(iftDataSet *Z);

/**
 * @brief Select samples to train, respecting a priori distribution of classes. See also iftSelectUnsupTrainSamples().
 * The function marks each sample in the dataset with a @a status of TRAIN or TEST
 * @param Z Dataset to be splitted.
 * @param perc Percentage of training data (0,1).
 * @return Number of training samples.
 */
int iftSelectSupTrainSamples(iftDataSet *Z, float perc);

/**
 * @brief Select samples to train. See also iftSelectSupTrainSamples().
 * @note Adan Echemendia add the param marked_times to allow a compromise between the desired normal distribution and the time. The original value set for this function is 100 but this can compromise the runtime of the function working with big datasets
 * The function marks each sample in the dataset with a @a status of TRAIN or TEST
 * @param Z Dataset to be splitted.
 * @param perc Percentage of training data (0,1) or number of training samples (if perc is > 1).
 * @return Number of training samples.
 * @param marked_times Number of marked times needed by a sample to be selected
 */
int iftSelectUnsupTrainSamples(iftDataSet *Z, float perc, int marked_times);

/**
 * @brief Select samples to train. According to a weighted sample selection. See also iftSelectUnsupTrainSamples().
 * The function marks each sample in the dataset with a @a status of TRAIN or TEST.
 * Each sample has a @a weight field to determine the weight to be considered.
 * @param Z Dataset to be splitted.
 * @param perc Percentage of training data (0,1).
 * @return Number of training samples.
 */
int iftSelectUnsupTrainSamplesByWeight(iftDataSet *Z, float train_perc);

void iftCopyClassifResult(iftDataSet *Z1, iftDataSet *Z2);

iftFeatSpaceParam iftComputeOverallMeanAndStdevFromDatasets(iftDataSet **Zs, int num_Z);

/**
 * @brief Normalize the dataset features. If a FeatSpaceParam struct is passed as a parameter it will not compute the
 * mean and the standard deviation of the dataset and get this values from the struct, otherwise pass NULL
 * @param Z Dataset to be normalized.
 * @param fsp Pointer to a feature space parameter struct. Pass NULL if we need to compute the mean and standard
 * deviation from the dataset samples.
 * @modified by Adan Echemendia
 * @return The normalized dataset.
 */
iftDataSet *iftNormalizeDataSet(iftDataSet *Z, iftFeatSpaceParam *fsp);

/**
 * @brief Normalize the dataset samples. Creating unit vector samples.
 * @param Z Dataset to be normalized.
 */
void iftNormalizeSamples(iftDataSet *Z);

/**
 * @brief Transforms the basis of dataset features according to the transformation matrix T
 * @date Jul, 2016
 * @author Peixinho
 * @param Z Input dataset
 * @param T Transformation matrix
 * @return Dataset in the new basis.
 */
iftDataSet* iftDataSetTransformBasis(const iftDataSet *Z, const iftMatrix *T);

iftDataSet *iftNormalizeContrastDataSet(iftDataSet *Z);

iftDataSet *iftCentralizeDataSet(iftDataSet *Z);

/**
 * @brief Normalize the dataset samples. Creating unit vector samples.
 * @param Z Dataset to be normalized.
 * @return The normalized dataset.
 */
iftDataSet *iftUnitNormDataSet(iftDataSet *Z);

/**
 * @brief Removes ambiguous samples in a dataset, by looking for pair of samples with distance 0.
 * @param Z The dataset.
 * @return Returns a dataset free of the ambiguous samples.
 */
iftDataSet *iftEliminateAmbiguousSamples(iftDataSet *Z);

/**
 * @brief Creates a new iftDataSet with samples that belongs to the given class.
 *
 * This function creates a new dataset with all the samples that belongs to the specified class.
 * The samples are copied in a new dataset without any change to the original one.
 *
 * @param Z Original dataset.
 * @param truelabel The class label to be selected.
 * @return A dataset only containing samples from the specific class.
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftExtractClass(iftDataSet *Z, int truelabel);

/**
 * TODO: No idea why this exists.
 */
iftDataSet *iftExtractObjectClasses(iftDataSet *Z);

/**
 * Creates a new dataset with all the samples with an specified status.
 * @param Z The original dataset
 * @param status The specified status
 * @return A dataset with a copy of the specified samples.
 * @author Peixinho
 * @date May, 2016
 */
iftDataSet *iftExtractSamples(const iftDataSet *Z, iftSampleStatus status);

/**
 * @brief Joins an array of Datasets copying or just pointing their Sample Feature Vectors.
 * @author Samuel Martins
 * @date Jun 6, 2016
 *
 * @warning If <b>copy_feats</b> is false, the feature vectors are only assigned/pointed to the original ones. \n
 * Then, BE CAREFUL when destroying this dataset.
 * 
 * @param  Z_array    Array of Datasets to be merged.
 * @param  n_datasets Number of Datasets
 * @param copy_feats  If true, it copies the sample feature vectors. Otherwise, it only points them.
 * @return            The joined dataset.
 */
iftDataSet *iftJoinDataSets(iftDataSet **Z_array, int n_datasets, bool copy_feats);


/**
 * @brief Joins an array of Datasets by sampling/selecting only some samples from them and copying or
 *        just pointing their Sample Feature Vectors.
 * @author Samuel Martins
 * @date Jun 9, 2016
 *
 * Each dataset must be sampled previously. The selected/chosen samples of each dataset have status TRAIN
 * in the corresponding sampler from the sampler array. \n
 * Only the first iteration of each sampler is considered.
 *
 * @warning If <b>copy_feats</b> is false, the feature vectors are only assigned/pointed to the original ones. \n
 * Then, BE CAREFUL when destroying this dataset.
 * 
 * @param  Z_array    Array of Datasets to be merged.
 * @param  n_datasets Number of Datasets.
 * @param  samplers   Array of samplers, one for each dataset, with size <b>n_datasets</b>.
 * @param  copy_feats If true, it copies the sample feature vectors. Otherwise, it only points them.
 * @return            The joined dataset.
 */
iftDataSet *iftJoinDataSetsBySampling(iftDataSet **Z_array, int n_datasets, iftSampler **samplers,
                                      bool copy_feats);


/**
 * @brief Counts the number of different classes present in the dataset. Also updates the ::iftDataSet.nclass variable.
 * @author Peixinho
 * @date Aug, 2015
 * @param Z Input dataset.
 * @return Number of classes in the dataset.
 */
int iftCountNumberOfClassesDataSet(iftDataSet *Z);




/**
 * @brief Gets the Minimum and Maximum weight values from a DataSet.
 * @author Samuel Martins
 * @date Jul 2, 2016
 * 
 * @param Z          Dataset to be analized.
 * @param min_weight Returns by Reference the min weight.
 * @param max_weight Returns by Reference the max weight.
 */
void iftGetMinMaxWeightInDataSet(const iftDataSet *Z, float *min_weight, float *max_weight);


/**
 * @brief Normalize the testing dataset according to the data distribution found in the training dataset.
 * @param Z1 Training dataset.
 * @param Z2 Testing dataset.
 * @return Normalized testing dataset.
 */
iftDataSet *iftNormalizeTestDataSet(iftDataSet *Z1, iftDataSet *Z2);

iftDataSet *iftNormalizeTestDataSet2(iftDataSet *Z, iftFeatSpaceParam fsp);

/**
 * @brief Centralize the testing dataset according to the data distribution found in the training dataset.
 * @param Z1 Training dataset.
 * @param Z2 Testing dataset.
 * @return Centralized testing dataset.
 */
iftDataSet *iftCentralizeTestDataSet(iftDataSet *Z1, iftDataSet *Z2);

/**
 * Multiplies the dataset features by a scalar.
 * @param Z
 * @param scalar
 * @date May, 2016
 * @author Peixinho
 */
void iftMultDataSetByScalar(iftDataSet *Z, float scalar);

/**
 * Creates a dataset form a multiband kernel.
 * @param K Input kernel.
 * @return The dataset.
 * @sa iftDataSetToMMKernel()
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftMMKernelToDataSet(iftMMKernel *K);

/**
 * Creates a multiband kernel from a dataset.
 * @param Z Input dataset.
 * @return The dataset.
 * @sa iftDataSetToMMKernel()
 * @date May, 2016
 * @author Peixinho
 */
iftMMKernel *iftDataSetToMMKernel(iftDataSet *Z);

/**
 * @brief Apply the PCA dimensionality reduction in the dataset features.
 * @param Z The original dataset.
 * @param num_of_comps Number of principal components.
 * @return The reduced dataset.
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftTransFeatSpaceByPCA(iftDataSet *Z, int num_of_comps);

/**
 * @brief Whiten the dataset through PCA.
 * @param Z The original dataset.
 * @return The whitened dataset.
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftWhiteningTransform(iftDataSet *Z);

/**
 * @brief Compute and apply the Supervised PCA dimensionality reduction in the dataset features.
 * @param Z The original dataset.
 * @param num_of_comps Number of principal components.
 * @return The reduced dataset.
 * @date May, 2016
 * @author Peixinho
 * @ingroup Machine Learning
 */
iftDataSet *iftTransFeatSpaceBySupPCA(iftDataSet *Z, int num_of_comps);

/**
 * @brief Given the PCA transform pre computed in Z1, applies the dimensionality reduction in Z2.
 * @param Z1 The PCA trained dataset.
 * @param Z2 Input dataset.
 * @param num_of_comps Number of principal components.
 * @return The reduced dataset.
 * @date May, 2016
 * @author Peixinho
 * @sa iftInverseTransformTestDataSetByPCA().
 */
iftDataSet *iftTransformTestDataSetByPCA(iftDataSet *Z1, iftDataSet *Z2);

/**
 * @brief Given the PCA transform pre computed in Z1, applies the inverse transform in Z2.
 * @param Z1 The PCA trained dataset.
 * @param Z2 Input dataset.
 * @return The restored dataset.
 * @date May, 2016
 * @author Peixinho
 * @sa iftTransformTestDataSetByPCA().
 */
iftDataSet *iftInverseTransformTestDataSetByPCA(iftDataSet *Z1, iftDataSet *Z2);

/**
 * @brief Converts a pixel based dataset into a label image. According to the specified data_type.
 * 
 * Labels are always values in the range [1, +inf) \n
 * If Z represents a Label Image, samples' labels must be decremented, because the background,
 * which will have label 1 in Z, must have value 0 in the resulting image. \n
 * Therefore, use decrement_labels = true 
 * 
 * @param Z Input dataset
 * @param decrement_labels If true, all samples' label values from Z will be decremented in resulting image. 
 * @return The label image.
 */
iftImage *iftDataSetToLabelImage(const iftDataSet *Z, bool decrement_labels);

/**
 * @brief Converts a pixel based dataset into a weight image. According to the specified data_type.
 * @param Z Input dataset
 * @return The weight image
 */
iftImage *iftDataSetWeight(iftDataSet *Z);


/**
 * @brief Computes the classification confusion matrix
 * @param dataset
 * @param normalized false to compute the actual number of correct and incorrect samples, and true to compute as a fraction of the total
 * @return Confusion matrix
 */
iftMatrix* iftComputeConfusionMatrix(iftDataSet *dataset, bool normalized);


/**
 * @brief It computes the confusion matrix of a dataset.
 *
 * Labels are always values in the range [1, +inf)
 *
 * @param Z is the dataset.
 * @param normalized indicates if or not the confusion matrix must be
 * normalized.
 */

//iftMatrix* iftComputeConfusionMatrix(iftDataSet *Z, bool normalized);


/**
 * @brief Set the status of all samples in the dataset.
 * @param Z The dataset.
 * @param status The status applied to the samples. {TRAIN, TEST}
 */
void iftSetStatus(iftDataSet *Z, iftSampleStatus status);

/**
 * @brief Set the status of the specified samples in the dataset.
 * @param Z The dataset.
 * @param set The set of samples to be changed.
 * @param status The status applied to the samples. {TRAIN, TEST}.
 */
void iftSetStatusForSet(iftDataSet *Z, iftSet *set, iftSampleStatus status);

/**
 * @brief Computes the covariance matrix of the dataset features.
 * @param Z Input dataset.
 * @return The Covariance matrix.
 * @date May, 2016
 * @author Peixinho
 */
iftMatrix *iftDatasetCovarianceMatrix(iftDataSet *Z);

/**
 * @brief Computes the PCA transformation matrix for a dataset.
 * @param Z Input dataset.
 * @return The transform matrix.
 * @date May, 2016
 * @author Peixinho
 * @warning All thefeatures are being used to compute the PCA. This dataset MUST contain only training data. See also iftExtractSamples().
 */
iftMatrix *iftRotationMatrixByPCA(iftDataSet *Z);

/**
 * @brief Computes the Rotation and Singular Value Decomposition matrices for the dataset features.
 * In the output rotation matrix, the eigenvectors are the columns.
 * @param Z Input dataset
 * @param S Output SVD matrix
 * @return Rotation matrix.
 * @date May, 2016
 * @author Peixinho
 */
iftMatrix *iftRotationMatrixAndSingularValuesByPCA(iftDataSet *Z, iftMatrix **S);

/**
 * @brief Shows usefull information about a dataset.
 * @param Z Dataset
 * @author Peixinho
 * @date Mar, 2017
 */
void iftPrintDataSetInfo(iftDataSet *Z, bool printStatisticsInfo);

/**
 * @brief Creates a (samples x features) Matrix from the dataset.
 * @param Z Input dataset
 * @return The dataset in matrix.
 * @sa iftFeatureMatrixToDataSet()
 * @author Peixinho
 * @date May, 2016
 */
iftMatrix *iftDataSetToFeatureMatrix(const iftDataSet *Z);

/**
 * @brief Creates a dataset (samples = rows, features = cols) from the input matrix.
 * @param X Input Matrix
 * @return The matrix turned into dataset.
 * @autor Peixinho
 * @date May, 2016
 */
iftDataSet *iftFeatureMatrixToDataSet(iftMatrix *X);

/**
 * @brief Computes the object principal axis, which is the direction of greatest variation.
 * @param label Image containing object.
 * @return The principal axis direction.
 * @date May, 2016
 * @author Peixinho
 */
iftVector iftPrincipalAxis(iftImage *label);

/**
 * Transform the dataset samples in norm one.
 * @param Z Input dataset
 * @return The normalized dataset
 */
iftDataSet *iftNormOneDataSet(iftDataSet *Z);

/**
 * @brief Creates a feature space parameter.
 * @author Peixinho
 * @date May, 2016
 */
iftFeatSpaceParam iftCreateFeatSpaceParam(void);

/**
 * Loads a PCA transform from disk.
 * @param filename The file path
 * @return The PCA feature space parameter
 * @author Peixinho
 * @date May, 2016
 */
iftFeatSpaceParam iftReadPCAMatrix(char *filename);

/**
 * Stores a PCA transform on disk.
 * @param filename The file path
 * @return The PCA feature space parameter
 * @author Peixinho
 * @date May, 2016
 */
void iftWritePCAMatrix(iftFeatSpaceParam fsp, char *filename);

/**
 * @brief Copy the feature space parameter.
 * @param Z
 * @return The space parameter copy.
 * @author Peixinho
 * @date May, 2016
 */
iftFeatSpaceParam iftCopyFeatSpaceParam(const iftDataSet *Z);

//TODO: Shouldnt distance table be able to store the distance between two different datasets?
iftDistanceTable *iftCreateDistanceTable(int nsamples1, int nsamples2);

void iftDestroyDistanceTable(iftDistanceTable **dt);

iftDistanceTable *iftReadDistanceTable(char *filename);

void iftPrintDistanceTable(iftDistanceTable *dt);

void iftWriteDistanceTable(iftDistanceTable *dt, char *filename);

iftDistanceTable *iftCompDistanceTable(const iftDataSet *Z1, const iftDataSet* Z2);

/**
 * @brief Compute the euclidean distance among samples in the dataset.
 * This is equivalent to call iftCompDistanceTable() with euclidean distance and alpha=1, but faster in large datasets.
 * @author Peixinho
 * @date March, 2016
 */
iftDistanceTable *iftCompEuclDistanceTable(const iftDataSet *Z1, const iftDataSet* Z2);

FILE **iftBuildDataSetsFromMImages(int nimages, int nclasses, int nfeatures, int nfilters,
                                   int stride, iftAdjRel *A, char *mimage_directory, char *output_directory,
                                   iftImageNames *img_names);

void iftExtractSamplesOfMImages(iftMImage *img, int truelabel, iftAdjRel *A, int sample_stride, FILE **fp);

void iftSwitchSamples(iftSample *s1, iftSample *s2, int nfeats);

/**
 * @brief Copies a sample (copying or only pointing its feature vector).
 * @author Samuel Martins
 * @date Jun 7, 2016
 * 
 * @param src        Source sample to be copied.
 * @param dst        Destination sample where the source sample will be copied.
 * @param n_feats    Number of Features (feat. vector size). Only necessary if <b>copy_feats</b> is true.
 * @param copy_feats If true, it copies the sample feature vectors. Otherwise, it only points them.
 */
void iftCopySample(iftSample *src, iftSample *dst, int n_feats, bool copy_feats);


void iftSwitchNotSVsPerErrorSamples(iftDataSet *Ztrain, iftDataSet *Ztest, int *not_SVs,
                                    int num_not_SVs, int *error_samples, int num_errors);

iftDataSet *iftSelectSamplesOfTheClass(iftDataSet *Z, int truelabel);

iftDataSet *iftSelectNegativeSamples(iftDataSet *Z, int positive_class);

/**
 * @brief Extracts samples with the specified truelabel
 * @author Peixinho
 * @param Z
 * @param labelClass
 * @return
 */
iftDataSet *iftExtractSamplesFromClass(const iftDataSet *Z, int labelClass);

/* Used by pyift */
void iftSetTrainingSupervoxelsFromSeeds(iftDataSet *dataset, iftImage *label_image, iftLabeledSet *seed);

/**
 * @brief Creates a dataset from a multi-band image. 
 * @param mimg A multi-image
 * @return A dataset
 */
iftDataSet *iftMImageToDataSet(iftMImage *mimg);

/**
 * @brief Computes a dataset from a multi-band image using as features the grayscale values of the bands of adjacent pixels.
 *
 * @author Thiago Vallin Spina
 * @date Feb 20, 2016
 *
 * @param img The input multi-band image.
 * @param A The adjacency relation corresponding to the patch around each voxel.
 * @return The dataset.
 *
 * @sa iftImageToDatasetUsingAdjacency
 */
iftDataSet *iftMImageToDataSetUsingAdjacency(iftMImage *mimg, iftAdjRel *A);

/**
 * Creates a dataset from multiband image voxels/pixels only in the region indicated by the mask.
 * @param mimg Input image
 * @param mask Binary mask
 * @return The dataset
 * @date May, 2016
 * @author Peixinho
 */
iftDataSet *iftMImageToDataSetInRegion(iftMImage *mimg, iftImage *mask);

/**
 * Creates a dataset from the multiband image edges.
 * @param mimg Input image.
 * @param A Adjacency of edges.
 * @return The edges dataset.
 */

iftDataSet *iftMImageToEdgesDataSet(iftMImage *mimg, iftAdjRel *A);

/**
 * TODO: edges here seems to have a diferent meaning than in iftMImageToEdgesDataSet, we have to define name conventions.
 */
iftDataSet *iftMImageToLabeledEdgesDataSet(iftMImage *mimg, iftImage *imgGT, iftImage *imglabels, iftAdjRel *A);

/**
 * @brief Converts edge based dataset into image.
 * @param dataset Input dataset
 * @param mimg input image.
 * @param A adjacency of edges.
 * @return The image.
 * @date May, 2016
 * @author Peixinho
 */
iftFImage *iftEdgesDataSetToFImage(iftDataSet *dataset, iftMImage *mimg, iftAdjRel *A);


/**
 * @brief Converts a pixel based dataset into an image.
 * @param Z Input dataset
 * @param img input image.
 * @return The image.
 * @date May, 2016
 * @author Peixinho
 */
iftImage *iftPixelDataSetToImage(iftDataSet *Z, iftImage *img);

/**
 * @brief Converts a pixel based dataset into a float image.
 * @param Z Input dataset
 * @param img input image.
 * @return The image.
 * @date May, 2016
 * @author Peixinho
 */
iftFImage *iftPixelDataSetToFImage(iftDataSet *Z, iftImage *img);

/**
 * @brief Defines the distance function for the dataset.
 * @param Z dataset
 * @param function_number Function number [1, 11].
 */
void iftSetDistanceFunction(iftDataSet *Z, int function_number);

/**
 * Computes and applies the PCA transform on the dataset.
 * @param Z Dataset
 * @return The aligned dataset.
 * @date May, 2016
 * @author Peixinho
 * @warning All thefeatures are being used to compute the PCA. This dataset MUST contain only training data. See also iftExtractSamples().
 */

iftDataSet *iftAlignDataSetByPCA(iftDataSet *Z);

/*---------------------- Distance functions -----------------------------*/

/** TODO: why not give a proper name for the distance functions???
 * @brief Default distance: Euclidean Distance
 * When alpha is either 0 or 1, it selects features for Euclidean distance. When alpha is from 0 to 1, it becomes a sort of weighted Euclidean distance.
 */
float iftDistance1(float *f1, float *f2, float *alpha, int n);

/**
 * @brief  Log of Euclidean Distance.
 */
float iftDistance2(float *f1, float *f2, float *alpha, int n);

/**
 * @brief In this function, alpha plays the exponents of the absolute differences between each feature. Features must be normalized within [0,1] to use this distance function and alpha >= 0.
 *
 */
float iftDistance3(float *f1, float *f2, float *alpha, int n);

/** @brief Log of iftDistance3(). */
float iftDistance4(float *f1, float *f2, float *alpha, int n);

/**
 * @brief Gaussian like distance.
 * Compute the distance as 1.0 - exp(eucl(f1-f2)/20) */
float iftDistance5(float *f1, float *f2, float *alpha, int n);


/**
 * @brief Inner-product-based distance for centralized datasets. Assumes feature vectors with norm 1.
 */
float iftDistance6(float *f1, float *f2, float *alpha, int n);

/**
 * @brief Minkowsky distance. Assumes feature vectors with norm 1.
 */
float iftDistance7(float *f1, float *f2, float *alpha, int n);

/**
 * @brief Compute the Distance based on the Correlation Coefficient of two vectors.
 * This function computes the Distance based on Correlation Coefficient value 0 <= dist(f1, f2) <= 1.
 * If f1 and f2 are linearly dependent, dist(f1, f2) = 0.
 * Reference: http://dl.acm.org/citation.cfm?id=507477
*/
float iftDistance8(float *f1, float *f2, float *alpha, int n);


/**
 * @brief Compute the Mutual Information Compression Index of two vectors.
 * This function computes the Mutual Information Compression Index value. 0 <= mici(f1,f2) <= 0.5*(var(x) + var(y)).
 * If mici(f1,f2) = 0, f1 and f2 are linearly related.
 * Reference: http://dl.acm.org/citation.cfm?id=507477
*/
float iftDistance9(float *f1, float *f2, float *alpha, int n);

/**
* @brief Chi-Square Distance: computes the distance between two distributions, used commonly to compare two histograms
* The Feature Vectors MUST be between [0,1] for better performance.
 * TODO: if I am not wrong the chi square requires norm one samples (histogram like)
**/
float iftDistance10(float *f1, float *f2, float *alpha, int n);

/**
* @brief Hellinger Distance: computes the distance between two distributions
**/
float iftDistance11(float *f1, float *f2, float *alpha, int n);

/*
 * @brief Computes the Bhattacharyya coefficient for discrete probability distributions (near 1 if the distributions
 * are very similar and near 0 otherwise).
 * */
float iftDistance12(float *f1, float *f2, float *alpha, int n);

/* initialize look-up table for Minkowski
				       distance computation. */
void iftSetMinkowskiTable(int mult_factor, float exp_factor);

void iftDestroyMinkowskiTable(void);

//TODO: Explain better each one of these following functions. I trully have no idea why they do what they do.

/**
 * @brief Pretty specific distance, only used in iftSupervoxelsToMeanDataSet().
 */
float iftDistMeanSizeSupervoxel(float *f1, float *f2, float *alpha, int n);

/**
 * @brief Pretty specific distance, only used in iftSupervoxelsToSelectiveSearchDataSet().
 */
float iftDistSelectiveSearchSupervoxel(float *f1, float *f2, float *alpha, int n);

float iftL1Norm(float *f1, float *f2, float *alpha, int n);

/**
 * @brief Pretty specific distance, only used in iftRegionMergingDataSet().
 */
float iftRegionMergingDist(float *f1, float *f2, float *alpha, int n);


/* --------------------- Merging functions ----------------------------- */

/**
 * TODO:This is just the sum of f1+f2. Move to iftCommon and transform in sum of two arrays.
 **/
float *iftMergeMeanSizeSupervoxel(float *f1, float *f2, float *alpha, int n);

/**
 * @brief Merges the feature vector of two samples extracted by ift iftSupervoxelsToSelectiveSearchDataset().
 * @warning Bad designed function, only works with features extracted by iftSupervoxelsToSelectiveSearchDataset().
 */
float *iftMergeSelectiveSearchSupervoxel(float *f1, float *f2, float *alpha, int n);

/**
 * TODO: computes the mean of two vectors. Creates a function for that in iftCommon
 */
float *iftMergeVector(float *f1, float *f2, float *alpha, int n);

/** TODO: computes the sum of two arrays. Create a function for that in iftCommon */
float *iftRegionMergingFeat(float *f1, float *f2, float *alpha, int n);

/**
 * @brief Counts the number of samples with the specified status.
 * @param Z
 * @param status
 * @return The number of samples with the specified status.
 * @author Peixinho
 * @date May, 2016
 */
int iftCountSamples(const iftDataSet *Z, iftSampleStatus status);

/**
 * @brief Counts the number of samples that satisfies the specified status and true label.
 * @param Z
 * @param truelabel
 * @param status
 * @return The number of samples with the specified status and true label.
 * @author Peixinho
 * @date May, 2016
 */
int iftCountSamplesTrueLabel(const iftDataSet *Z, int truelabel, iftSampleStatus status);

/**
 * @brief Counts the number of samples that satisfies the specified status and label.
 * @param Z
 * @param label
 * @param status
 * @return The number of samples with the specified status and label.
 * @author Peixinho
 * @date May, 2016
 */
int iftCountSamplesLabel(const iftDataSet *Z, int label, iftSampleStatus status);

/**
* @brief Selects TRAIN and TEST samples according to the sampler strategy.
* @param fileSet FileSet to split.
* @param sampler Sampler strategy
* @param iteration Current sampling iteration.
* @author Peixinho
* @note A usage example can be found at @ref iftFileSetSampler.c
* @date May, 2016
*/
void iftSampleFileSet(iftFileSet *fileSet, const iftSampler *sampler, int iteration);

/**
 * @brief Creates a new FileSet with the specified status.
 * @param fileSet The fileSet to extract samples.
 * @param status The sample status to be copied.
 * @author Peixinho
 * @date May, 2016
 */
iftFileSet *iftExtractFileSamples(const iftFileSet *fileSet, iftSampleStatus status);


/**
 * @brief Joins the datasets <b>Z1</b> and <b>Z2</b> and returns it as a new dataset
 * @details The datasets need to have the same number of features, must to be unlabeled, must have the same arc weight function and must have the same referenced data
 * @param Z1 A first dataset
 * @param Z2 A second dataset
 * @return A new dataset
 */
iftDataSet *iftMergeDataSets(iftDataSet *Z1, iftDataSet *Z2);

/**
 * @brief Conquer the outliers(identified with label 0) in an MImage dataset thorough the IFT algorithm
 * @author Adan Echemendia
 * @param Z A dataset
 */
void iftMImageDatasetConqueringOutliersByIFT(iftDataSet *Z);

/**
 * @brief Transforms features by scaling each feature in the range [0...1]
 *
 * @details The transformation is given by
 *
 * X_scaled = (X - Xmin)/(Xmax - Xmin)
 *
 * @author Deangeli
 * @param dataset A dataset
 */
void iftMinMaxFeatureScale(iftDataSet *dataset);

iftDataSet* iftSupervoxelsToLabColorMeanStdAndSizeDataset(iftImage* image, iftImage* label_image);

iftDataSet* iftImageSupervoxelsToYCbCrDataset(iftImage* image, iftImage* label_image, char normalize);

iftDataSet* iftSupervoxelsToYCbCrMeanStdAndSizeDataset(iftImage* image, iftImage* label_image);

iftDataSet* iftImageSupervoxelsToGrayDataset(iftImage* gray_img, iftImage* label_image, char normalize);


void iftWriteDataSet(const iftDataSet *Z, const char *pathname);
void iftWriteField(FILE *fp, const char* field_name, void* field, size_t size, size_t count);

inline unsigned long long iftComputeUniqueIdForSampleInfo(int sampleTrueLabel, int sampleUniqueIdInClass){
    return (sampleTrueLabel*(100000000) + sampleUniqueIdInClass);
}

inline unsigned long long iftComputeUniqueIdForSample(iftSample sample){
    return (sample.truelabel*(100000000) + sample.id);
}


#ifdef __cplusplus
}
#endif


#endif

