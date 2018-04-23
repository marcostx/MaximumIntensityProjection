#ifndef IFT_DESCRIPTORS_H_
#define IFT_DESCRIPTORS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/stat.h>

#include "iftCommon.h"
#include "iftImage.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftRadiometric.h"

/**
 * @brief Feature vector
 * @details [long description]
 * @author Adan Echemendia
 * @date Jun 10, 2016
 */
typedef struct ift_features {

    /** List of feature values */
    float *val;
    /** Size of the feature vector */
    int    n;
} iftFeatures;

typedef struct ift_hog {
    iftVoxel cellSize;
    iftVoxel blockSize;
    iftVoxel blockStride;
    int nbins;
} iftHog;

/**
 * @brief Creates a feature vector
 * 
 * @param n Size of the feature vector
 * @return A feature vector
 */
iftFeatures *iftCreateFeatures(int n);

/**
 * @brief Destroys a feature vector
 * 
 * @param feat Pointer to pointer to a feature vector
 */
void         iftDestroyFeatures(iftFeatures **feat);

/**
 * @brief Reads a set of features from a file and builds a feature vector
 * @details [long description]
 * 
 * @param filename File path
 * @return A feature vector
 */
iftFeatures *iftReadFeatures(const char *filename);

/**
 * @brief Writes a feature vector to a file
 * 
 * @param feats Feature vector
 * @param filename File path
 */
void         iftWriteFeatures(const iftFeatures *feats, const char *filename);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param img [description]
 * @param A [description]
 * 
 * @return [description]
 */
iftImage    *iftLocalBinaryPattern(iftImage *img, iftAdjRel *A);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param img [description]
 * @param A [description]
 * 
 * @return [description]
 */
iftFeatures *iftExtractLBP(iftImage *img, iftAdjRel *A);

/**
 * @brief Builds a feature vector with the normalized brightness values of an image
 * @details [long description]
 * 
 * @param img Target image
 * @return Feature vector with the normalized brightness values of an image
 */
iftFeatures *iftExtractBrightness(iftImage *img);

/**
 * @brief Builds a feature vector with the normalized brightness values of the gradient magnitude image of image <br>img</br>
 * @details [long description]
 * 
 * @param img Target image
 * @param A Adjacency relation
 * 
 * @return Feature vector with the normalized brightness values of the gradient magnitude image of image <br>img</br>
 */
iftFeatures *iftExtractGradient(iftImage *img, iftAdjRel *A);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param img [description]
 * @param bins_per_band [description]
 * 
 * @return [description]
 */
iftFeatures *iftExtractBIC(iftImage *img, int bins_per_band);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param img [description]
 * @param mask [description]
 * @param bins_per_band [description]
 * @return [description]
 */
iftFeatures *iftExtractBICInMask(iftImage *img, iftImage *mask, int bins_per_band);

/**
 * @brief Computes the Manhattan distance between two feature vectors
 * @details [long description]
 * 
 * @param feat1 Feature vector one
 * @param feat2 Feature vector two
 * 
 * @return Manhattan distance value between the two feature vectors
 */
float        iftFeatDistL1(iftFeatures *feat1, iftFeatures *feat2);

/**
 * @brief Builds a feature vector with all values of a multi-band image
 * @details [long description]
 * 
 * @param img Target image
 * @return Feature vector with all values of a multi-band image
 */
iftFeatures *iftMImageToFeatures(iftMImage *img);

/**
 * @brief Writes a feature vector to a file. The first saved value is the true label of the sample.
 * @details [long description]
 * 
 * @param features Feature vector of the sample
 * @param truelabel True label of the sample
 * @param fp File path
 */
void     iftWriteFeaturesInFile(iftFeatures *features, int truelabel, FILE *fp);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param img [description]
 * @param u1 [description]
 * @param u2 [description]
 * @return [description]
 */
iftFeatures *iftIntensityProfile(iftImage *img, iftPoint u1, iftPoint u2);

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param img [description]
 * @param bins_per_band [description]
 * @param normalization_value [description]
 * @return [description]
 */
iftImage    *iftQuantize(iftImage *img, int bins_per_band, int normalization_value);
/**
 * @brief Creates and defines the HOG feature extractor parameters.
 * @param cellSize The cell size considered in the histogram computation.
 * @param blockSize The block size defines the number of cells to be considered in the normalization.
 * @param blockStride The block stride defines the overlap between blocks in the normalization.
 * @author Peixinho
 * @date April, 2016
 */
iftHog* iftCreateHog2D(int cellSize, int blockSize, int blockStride, int nbins);

/**
 * @brief Destroys a HOG
 * 
 * @param Pointer to pointer to a HOG
 */
void iftDestroyHog(iftHog** hog);

/**
 * @brief Extracts the Hog descriptor of an image.
 * @author Peixinho
 * @date April, 2016
 */
iftFeatures* iftExtractHOG2D(const iftHog *hog, iftImage *img);

#ifdef __cplusplus
}
#endif

#endif
