#ifndef IFT_METRICS_H
#define IFT_METRICS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftDataSet.h"


/**
 * @brief Type of the Metric Functions: Distance or Similarity Functions
 * @author Samuel Martins
 * @date Mar 23, 2016
 * @ingroup Metrics
 */
typedef enum {
    IFT_SIMILARITY_FUNCTION,
    IFT_DISTANCE_FUNCTION
} iftMetricType;


//True positives, false positives, false negatives, true negatives.
typedef struct ift_error_classification {
	int tp;
	int fp;
	int fn;
	int tn;
} iftErrorClassification;

float iftPrecisionGivenErrors(iftErrorClassification* error);
float iftRecallGivenErrors(iftErrorClassification* error);

float iftAccuracyGivenErrors(iftErrorClassification* error);
float iftFScoreGivenErrors(iftErrorClassification* error);
float iftProbabilityRatioGivenErrors(iftErrorClassification* error);

/**
 * Computes the median value for a float array.
 * @warning In case of an even size array, this method computes the average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @return median of array
 */
float iftMedianFloatArray(float* array, int length);

/**
 * Computes the weighted median value for a float array.
 * @warning In case of an even size array, this method computes the average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @param weight weight values (must be greater than zero)
 * @return median of array
 */
float iftWeightedMedianFloatArray (float *array, int length, float* weight);

/**
 * Computes the weighted average value for a float array.
 * @warning In case of an even size array, this method computes the average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @param weight weight values (must be greater than zero)
 * @return median of array
 */
float iftWeightedMeanFloatArray (float *array, int length, float* weight);

/**
 * Computes the median value for an int array.
 * @warning In case of an even size array, this method computes the rounded average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @return median of array
 */
int iftMedianIntArray(int* array, int length);

float iftMinFloatArray(float* array, int length);
float iftMaxFloatArray(float* array, int length);
float iftSumFloatArray (float *array, int length);
int iftSumIntArray (int *array, int length);
float iftMeanFloatArray (float *array, int length);
void iftMultScalarFloatArray(float* array, int length, float scalar);
float iftMeanFloatArrayDiffZero(float *array, int length);
float iftVarianceFloatArray (float *array, int length);
float iftStddevFloatArray (float *array, int length);

/**
 * Computes the cohen kappa score on dataset.
 * @param Z
 * @return score
 * @author Peixinho
 * @date Aug, 2017
 */
float iftCohenKappaScore(iftDataSet* Z);

float iftFscore(iftDataSet *Z, int PositiveClass);
float iftTruePositives(iftDataSet *Z);
float iftSkewedTruePositives(iftDataSet *Z);
float iftNormAccuracy(iftDataSet* Z);
float iftClassifAccuracy(iftDataSet *Z);

float iftOPFAccuracy(iftDataSet *Z);

float iftMeasureERR(float *negatives, int total_neg, float *positives, int total_pos, float threshold);
void  iftMeasureFarFrr(iftDataSet *Z, float *FAR, float *FRR);
iftErrorClassification iftMeasureErrorClassification(iftDataSet *Ztest);
float iftWeightedAccuracyBinaryClasses(iftErrorClassification cm);
float iftMeasureMCC(iftErrorClassification cm);

float iftCosineDistance(float *f1, float *f2, int n);
float iftCosineDistance2(float *f1, float *f2, int n);

float iftCorrCoef(float *x, float *y, int n);
float iftCorrCoefDist(float *x, float *y, int n);
float iftMICI(float *x, float *y, int n);


/**
 * @brief It Finds the Most Similar Sample from a Score Matrix.
 * @author Samuel Martins
 * @date Mar 23, 2016
 * @ingroup Metrics
 *
 * Each row i=p and column j=p from the Score Matrix corresponds to the sample.\n
 * The Most Similar Sample is the one whose the sum of all scores is MAXIMUM (for SIMILARITY functions) or MINIMUM (for DISTANCE functions).\n
 *
 * @param score_matrix Score Matrix that will used to figure out the Most Similar Image.
 * @param metric_type Type of the Score Function that generated the considered Score Matrix.
 * @param out_best_score Reference of the variable which will store the score of the Most Similar Sample. If NULL, nothing is returned to it.
 * @return The index of the Most Similar Sample.
 *
 * @note If <b>best_score</b> is NULL, nothing is returned to it.
 */
int iftMostSimilarSample(const iftMatrix *score_matrix, iftMetricType metric_type, double *out_best_score);


#ifdef __cplusplus
}
#endif

#endif
