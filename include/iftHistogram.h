#ifndef IFT_HISTOGRAM_H_
#define IFT_HISTOGRAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"

typedef struct ift_histogram {
  double *val;
  int nbins;
} iftHistogram;

iftHistogram *iftCreateHistogram(int nbins);
void          iftDestroyHistogram(iftHistogram **hist);
iftHistogram *iftReadHistogram(char *filename);
void          iftWriteHistogram(iftHistogram *hist, char *filename);


/**
 * @brief Compute a Histogram for the Integer Values from an integer array. Return a Normalized Histogram if
 * <b>normalize</b> is true.
 * @author Thiago Vallin Spina
 * @date May 12, 2017
 * @ingroup Histogram
 *
 * @param  val The input integer valued array.
 * @param  n   The size of the input integer valued array.
 * @param  nbins     Number of bins from the Output Histogram.
 * @param  normalize If true, the histogram is normalized.
 * @return           The computed histogram from the gray values.
 */
iftHistogram *iftHistogramOfIntArray(const int *val, int n, int nbins, bool normalize);

/**
 * @brief Computes the histogram of integer values for the elements marked by a given binary array.
 *
 * @author Thiago Vallin Spina
 * @date May 12, 2017
 *
 * @param val The input integer valued array.
 * @param n   The size of the input integer valued AND binary arrays.
 * @param bin The binary array.
 * @param nbins The number of bins into which the histogram should be divided. Note: the maximum intensity of the array
 *              is considered for defining the bins.
 * @param normalize If true, normalizes the histogram by dividing it by the number of pixels in the binary mask.
 *
 * @return The computed histogram.
 *""""
 */
iftHistogram *iftHistogramOfIntArrayInRegion(const int *val, int n, const int *bin, int nbins, bool normalize);

/**
 * @brief Computes the histogram of integer values for the elements with a given label.
 * @author Thiago Vallin Spina
 * @date May 12, 2017
 *
 * @param val The input integer valued array.
 * @param n   The size of the input integer valued AND label arrays.
 * @param label The label array.
 * @parm obj_label  The label from the required object whose element will be the only ones considered by the histogram.
 * @param nbins     The number of bins into which the histogram should be divided. Note: the maximum intensity of the array
 *                  is considered for defining the bins.
 * @param normalize If true, normalizes the histogram by dividing it by the number of elements in the label array.
 * @return          The computed histogram.
 */
iftHistogram *iftHistogramOfIntArrayInLabelRegion(const int *val, int n, const int *label, int obj_label,
                                                  int nbins, bool normalize);


/**
 * @brief Compute a Histogram for the Float Values from a float array. Return a Normalized Histogram if
 * <b>normalize</b> is true.
 * @author Thiago Vallin Spina
 * @date May 12, 2017
 * @ingroup Histogram
 *
 * @param  val The input float valued array.
 * @param  n   The size of the input float valued array.
 * @param  nbins     Number of bins from the Output Histogram.
 * @param  normalize If true, the histogram is normalized.
 * @return           The computed histogram from the gray values.
 */
iftHistogram *iftHistogramOfFloatArray(const float *val, int n, int nbins, bool normalize);

/**
 * @brief Compute a Histogram for the Gray Values from an Image. Return a Normalized Histogram if
 * <b>normalize</b> is true.
 * @author Samuel Martins
 * @date May 12, 2016
 * @ingroup Histogram
 *  
 * @param  img       Input Image.
 * @param  nbins     Number of bins from the Output Histogram.
 * @param  normalize If true, the histogram is normalized.
 * @return           The computed histogram from the gray values.
 */
iftHistogram *iftGrayHistogram(const iftImage *img, int nbins, bool normalize);


/**
 * @brief Computes the histogram of pixel intensities for the pixels inside a given binary mask.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param img The input grayscale image.
 * @param bin The binary mask.
 * @param nbins The number of bins into which the histogram should be divided. Note: the maximum intensity of the image
 *              is considered for defining the bins.
 * @param normalize If true, normalizes the histogram by dividing it by the number of pixels in the binary mask.
 *
 * @return The computed histogram.
 *""""
 */
iftHistogram *iftGrayHistogramInRegion(iftImage *img, iftImage *bin, int nbins, bool normalize);

/**
 * @brief Computes the histogram of pixel intensities for the pixels with a given label.
 * @author Samuel Martins
 * @date Apr 14, 2016
 *
 * @param img       The input grayscale image.
 * @param label_img The Label Image.
 * @parm obj_label  The label from the required object whose pixel will be the only ones considered by the histogram.
 * @param nbins     The number of bins into which the histogram should be divided. Note: the maximum intensity of the image
 *                  is considered for defining the bins.
 * @param normalize If true, normalizes the histogram by dividing it by the number of pixels in the binary mask.
 * @return          The computed histogram.
 */
iftHistogram *iftGrayHistogramInLabelRegion(const iftImage *img, const iftImage *label_img,
                                            int obj_label, int nbins, bool normalize);


/**
 * @brief Computes a histogram for each Label from <b>label_img</b> from 1 to n.
 * @author Samuel Martins
 * @date Jun 18, 2016
 *
 * @note The histogram i is from the label i+1.
 * 
 * @param img       The input grayscale image.
 * @param label_img The Label Image.
 * @param nbins     The number of bins into which the histogram should be divided. Note: the maximum intensity of the image
 *                  is considered for defining the bins.
 * @param normalize If true, normalizes the histogram by dividing it by the number of pixels in the binary mask.
 * @return          An array of histograms, one for each label.
 */
iftHistogram **iftGrayHistogramsForLabels(const iftImage *img, const iftImage *label_img, int nbins, bool normalize);


/**
 * @brief Computes the histogram of pixel intensities for the pixels in the given set.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param img The input grayscale image.
 * @param pixels The pixel set that constrains the computation of the histogram.
 * @param nbins The number of bins into which the histogram should be divided. Note: the maximum intensity of the image
 *              is considered for defining the bins.
 * @param normalize If true, normalizes the histogram by dividing it by the number of pixels in the set.
 *
 * @return The computed histogram.
 *
 */
iftHistogram *iftGrayHistogramForSet(iftImage *img, iftSet *pixels, int nbins, char normalize);


/**
 * @brief Compute a Histogram for the RGB Values from a Color Image. Return a Normalized Histogram if
 * <b>normalize</b> is true.
 * @author Samuel Martins
 * @date May 12, 2016
 * @ingroup Histogram
 *  
 * @param  img       Input Image.
 * @param  nbins     Number of bins from the Output Histogram.
 * @param  normalize If true, the histogram is normalized.
 * @return           The computed histogram from the gray values.
 */
iftHistogram *iftColorHistogram(const iftImage *img, int nbins, bool normalize);


iftHistogram *iftAccHistogram(iftHistogram *hist);
int           iftHistogramMode(iftHistogram *hist, char exclude_zero);


/**
 * @brief Computes the mean value of a histogram.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The mean value
 *
 */
double iftHistogramMean(iftHistogram *hist, bool exclude_zero);

/**
 * @ brief Computes the standard deviation around a given value.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist The input histogram
 * @param value The value (e.g., the mean).
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The standard deviation around the value
 */
double iftHistogramStdevAroundValue(iftHistogram *hist, double value, char exclude_zero);

/**
 * @brief Computes the median value of a histogram.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The median value
 *
 */
int iftHistogramMedian(iftHistogram *hist, char exclude_zero);

/**
 * @brief Computes the median absolute deviation of a histogram given a value, which we suppose is the median.
 *
 * The median absolute deviation is the equivalent of the standard deviation for the median. It is a robust statistic
 * as the median. For more information, see https://en.wikipedia.org/wiki/Median_absolute_deviation.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist The input histogram
 * @param value The value (e.g., the median of the histogram).
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The median absolute deviation value
 *
 */

int iftHistogramMedianAbsoluteDeviationAroundValue(iftHistogram *hist, int value, char exclude_zero);


/**
 * @brief Gets the Arg (Bin) Min of a histogram. For ties, the first arg min value found is selected. 
 * @author Samuel Martins
 * @date Jun 2, 2016
 *
 * Ex: For a histogram with 256 values/bins of a Gray Image, the value/bin (in [0, 255]) with the lowest incidence is selected.
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The arg (bin) min value.
 *
 */
int iftHistogramArgMin(iftHistogram *hist, bool exclude_zero);


/**
 * @brief Gets the Arg (Bin) Max of a histogram. For ties, the first arg max value found is selected. 
 * @author Samuel Martins
 * @date Jun 2, 2016
 *
 * Ex: For a histogram with 256 values/bins of a Gray Image, the value/bin (in [0, 255]) with the highest incidence is selected.
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The arg (bin) max value.
 *
 */
int iftHistogramArgMax(iftHistogram *hist, bool exclude_zero);



/**
 * @brief Gets the Bin of minimum bin value that has at least one incidence in the histogram.
 * @author Samuel Martins
 * @date Jun 2, 2016
 *
 * Ex: [0, 0, 40, 20] - the chosen bin is the bin 2 (bin[2] = 40)
 *
 * @param hist         The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return             The chosen bin.
 *
 */
int iftHistogramMinBinWithIncidences(iftHistogram *hist, bool exclude_zero);


/**
 * @brief Gets the Bin of maximum bin value that has at least one incidence in the histogram.
 * @author Samuel Martins
 * @date Jun 2, 2016
 *
 * Ex: [0, 0, 40, 20] - the chosen bin is the bin 3 (bin[3] = 20)
 *
 * @param hist         The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return             The chosen bin.
 *
 */
int iftHistogramMaxBinWithIncidences(iftHistogram *hist, bool exclude_zero);



/**
 * @brief Adds the bins of two histograms.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist1 The first histogram.
 * @param hist2 The second histogram.
 * @return The added histogram.
 */
iftHistogram *iftAddHistograms(iftHistogram *hist1, iftHistogram *hist2);


/**
 * @brief Adds the histogram src to dst.
 * @author Samuel Martins
 * @author Thiago Vallin Spina
 * @date APr 14, 2016
 *
 * @param src The source histogram.
 * @param dst The destination histogram.
 */
void iftAddHistogramsInPlace(iftHistogram *src, iftHistogram *dst);


/**
 * @brief Prints a histogram on console.
 * @author Samuka
 * @date Jan 3, 2017
 */
void iftPrintHistogram(const iftHistogram *hist);


#ifdef __cplusplus
}
#endif

#endif
