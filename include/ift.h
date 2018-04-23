#ifndef IFT_H_
#define IFT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftActiveLearning.h"
#include "iftAdjacency.h"
#include "iftAdjSet.h"
#include "iftBasicDataTypes.h"
#include "iftBMap.h"
#include "iftBagOfFeatures.h"
#include "iftClassification.h"
#include "iftActiveLearning.h"
#include "iftClustering.h"
#include "iftColor.h"
#include "iftCommon.h"
#include "iftGenericLinkedList.h"
#include "iftCompression.h"
#include "iftCompTree.h"
#include "iftCSV.h"
#include "iftCurve.h"
#include "iftDataSet.h"
#include "iftDeepLearning.h"
#include "iftDescriptors.h"
#include "iftDHeap.h"
#include "iftDialog.h"
#include "iftDicom.h"
#include "iftDict.h"
#include "iftDisjointSet.h"
#include "iftDList.h"
#include "iftFHeap.h"
#include "iftFunctions.h"
#include "iftFIFO.h"
#include "iftFile.h"
#include "iftFiltering.h"
#include "iftFImage.h"
#include "iftFImageForest.h"
#include "iftFSet.h"
#include "iftGeometric.h"
#include "iftGQueue.h"
#include "iftGraphics.h"
#include "iftHierarchicCluster.h"
#include "iftHistogram.h"
#include "iftIGraph.h"
#include "iftImage.h"
#include "iftImageForest.h"
#include "iftImageMath.h"
#include "iftInpainting.h"
#include "iftInterpolation.h"
#include "iftIntQueue.h"
#include "iftJson.h"
#include "iftKernel.h"
#include "iftKmeans.h"
#include "iftLabeledSet.h"
#include "iftLIFO.h"
#include "iftList.h"
#include "iftMathMorph.h"
#include "iftMatrix.h"
#include "iftMemory.h"
#include "iftMetrics.h"
#include "iftMImage.h"
#include "iftMSPS.h"
#include "iftMST.h"
#include "iftObjectModels.h" 
#include "iftParamOptimizer.h"
#include "iftParamOptimizationProblems.h"
#include "iftParseInput.h"
#include "iftPlane.h"
#include "iftRadiometric.h"
#include "iftReconstruction.h"
#include "iftRegion.h"
#include "iftRegistration.h"
#include "iftRepresentation.h"
#include "iftSeeds.h"
#include "iftSegmentation.h"
#include "iftSegmentationResuming.h"
#include "iftSet.h"
#include "iftSimilarity.h"
#include "iftSlic.h"
#include "iftSList.h"
#include "iftSort.h"
#include "iftString.h"
#include "iftSVM.h"
#include "iftUtil.h"
#include "iftVideo.h"
#include "iftManifold.h"

//generic data structures. Private use only in order to keep Type safety.
//#include "iftGenericVector.h"
//#include "iftGenericVectorUtil.h"
//#include "iftGenericMatrix.h"
//#include "iftGenericMatrixUtil.h"
//#include "iftArgumentList.h"
//#include "iftGenericLinkedList.h"
//#include "iftDataTransference.h"

#undef I


/**
 * @defgroup Basic
 * @brief Structures and functions used for basic operations.
 * @{
 *
 * @defgroup DataTypes
 * @brief Definition of Basic Data Types.
 * @note Programs:
 * * @note iftTestGVal.c = It shows how to use the Generic Value (iftGVal)
 * @{ @}
 * 
 * @defgroup Dialog
 * @brief Functions used to show Dialog Messages.
 * @{ @}
 *
 * @defgroup Geometry
 * @brief Analytic Geometry operations.
 * @{ @}
 *
 * @defgroup Memory
 * @brief Functions and definitions for memory management.
 * @{ @}
 *
 * @defgroup Numeric
 * @brief Functions about Numeric operations.
 * @{ @}
 *
 * @defgroup Regex
 * @brief Regular Expression functions.
 * @note Programs:
 * * @ref iftRegexExamples.c = Examples and a quick explanation about Regex
 * * @ref iftTestRegex.c = Checks if an expression/string matches with a Regex
 * @{ @}
 *
 * @defgroup String
 * @brief String manipulation functions.
 * @note Programs:
 * * @ref iftTestSplitString.c = This program splits a string at a given position
 * * @ref iftTestReplaceString.c = It returns a copy of the string in which the occurrences of old substring have been replaced
 * with new substring.
 * @{ @}
 * 
 * @defgroup Time
 * @brief A set of functions to compute and show the runtime of programs and/or code blocks.
 * @{ @}
 * 
 * @} 
 */
/**
 * @defgroup DataStructures
 * @brief Module with Data Structure definitions and functions.
 * @{
 *
 * @defgroup BitMap
 * @brief Bit Map structure.
 * @{ @}
 * 
 * @defgroup Dictionary
 * @brief Dictionary structure, based on Hash Tables, that stores pairs of generic keys and values.
 *
 * @note Examples can be found in demo/Miscellaneous/DataStructures/iftTestDict*.c
 *
 * The allowed Datatypes for the keys are: <b>IFT_CHAR_TYPE, IFT_UCHAR_TYPE, IFT_STR_TYPE, SHORT_TYPE, INT_TYPE, IFT_LONG_TYPE, USHORT_TYPE, INT_TYPE, IFT_LONG_TYPE,
 * FLOAT_TYPE, IFT_DBL_TYPE, IFT_PTR_TYPE.</b>\n
 * The allowed Datatypes for the value are the same of the key plus <b>IFT_PTR_TYPE</b>\n.
 * Initially, empty buckets of the dictionary has key type and val type equals to IFT_UNTYPED.
 * @{ @}
 *
 * @defgroup Image
 * @brief Image definition and basic processing functions.
 * @note <b>Programs:</b>
 * * @ref iftExtractROI.c = Extract ROIs of an Image
 * * @ref iftInsertROI.c = Inserts/Copies an Image (Region Of Interest (ROI)) inside a Target Image.
 * * @ref iftExtractObject.c = Extract an Object from an Image
 * * @ref iftExtractAllObjects.c = Extract All Objects (Labels) from an Image
 * * @ref iftExtractImageROI.c = Extract the Min. Bounding Box (ROI) from an Image with all objects/voxels
 *                               (non-zero) inside it. \n
 * *
 * * @ref iftLabelImageAreaVolume.c = Computes the area/volume from each Object from an Label Image. 
 * * @ref iftSetVoxelSize.c = Overwrites the Pixel/Voxel Sizes of an Image or a Set of Images to some given Values.
 * @{ @}
 *
 * @defgroup StringDoublyLinkedList
 * @brief String Doubly Linked List definitions and functions.
 * @note Programs:
 * * @ref iftTestSList.c = It shows how to use the String Doubly Linked List.
 * @{ @}
 * 
 * @}
 */

/**
 * @defgroup MachineLearning
 * @brief Algorithms and data structures for machine learning and pattern recognition.
 * @{
 * @defgroup Supervised
 * @brief Labeled based algorithms.
 * @{ @}
 * @defgroup Algorithms not relying on labeled information.
 * @brief .
 * @{ @}
 * @}
 */

/**
 * @defgroup Dataset
 * @brief All Management function and definitions of Dataset.
 * @{
 *
 * @defgroup Metrics
 * @brief Several functions used to compute Metrics between Samples, Feature Vectors, etc, from a Dataset.
 * @{ @}
 *
 * @}
 */
/**
 * @defgroup ImageProcessing
 * @brief Image processing functions.
 * @{ 
 * 
 * @defgroup Filtering
 * @brief Several Image Filters.
 * @{ @}
 * 
 * @defgroup ImageConverters
 * @brief Several Image Converter functions.
 * @{ @}
 *
 * @defgroup ImageMetrics
 * @brief Functions that computes some metric (distance or similarity) between images.
 * @{ @}
 *
 * @defgroup Interpolation
 * @brief Image reslicing and other interpolation functions.
 * @{ @}
 *
 * 
 * @defgroup Registration
 * @brief Image Registration functions and approaches.
 * @{
 *
 * @defgroup Elastix
 * @brief Functions and definitions of the Elastix Registration approach.
 * @{ @}
 * 
 * @}
 *
 * @}
 */
/**
 * @defgroup ObjectModels
 * @brief Functions and definitions of any type of Object Model, such as Prob. Atlas, Fuzzy models, etc.
 * @{
 *
 * @defgroup Atlas
 * @brief Functions and Definitions about Probabilistic Atlases.
 * @note Programs
 * * @ref iftBuildProbAtlas.c = Builds a Probabilistic Atlas from a set of Images and their Labels.
 * @{ @}
 *
 * @}
 */
/**
 * @defgroup Utilities
 * @{
 *
 * @defgroup Allocation
 * @brief Common functions to allocate memory.
 * @{ @}
 *
 * 
 * @defgroup CommandLineParser
 * @brief Common functions to parser and manipulate Command Line Arguments from the Input.
 * @note Programs
 * * @ref iftTestCmdLineParser1.c, @ref iftTestCmdLineParser2.c = Show how to parse input command line
 *   arguments from the keyboard
 * @{ @}
 *
 * 
 * @defgroup CSV
 * @brief Generic CSV manipulations.
 * @note Programs:
 * * @ref iftTestCSV.c = It shows how to open a CSV file.
 * @{ @}
 *
 * 
 * @defgroup File
 * @brief File and path management functions.
 * @note Programs:
 * * @ref iftTestLoadDir.c = It loads a directory.
 * * @ref iftMakeDirsExamples.c = It shows examples of how to create Pathnames, Files, and Directories.
 * @{ @}
 *
 * 
 * @defgroup Messages
 * @brief Console user messages.
 * @{ @}
 *
 * @defgroup Compression
 * @brief Allows the reading and writing of compressed data into/from files.
 * @{ @}
 *
 * @}
 */

#ifdef __cplusplus
}
#endif

#endif
