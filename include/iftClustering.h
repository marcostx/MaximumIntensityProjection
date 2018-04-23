/**
 * @file iftClustering.h
 * @brief Definitions and functions about clustering methods.
 * @author Adan Echemendia
 * @date May 17, 2016
 *
  */
#ifndef IFT_CLUSTERING_H_
#define IFT_CLUSTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftDataSet.h"
#include "iftAdjSet.h"
#include "iftSet.h"
#include "iftFHeap.h"
#include "iftSet.h"
#include "iftBMap.h"
#include "iftMSPS.h"
#include "iftMST.h"
#include "iftSort.h"
#include "iftMetrics.h"

/**
 * @brief Node of a KNN Graph
 * @author Adan Echemendia
 * @date May 17, 2016
  */
typedef struct ift_knnnode {
    /**  Maximum arc weight from the node to its neighbors */
    float maxarcw;
    /** Corresponding root node in the graph */
    int   root;
    /** Corresponding training sample in the original dataset */
    int   sample;
    /** predecessor node*/
    int pred;
    /** List of adjacent nodes */
    iftAdjSet *adj;
    /** List of adjacent nodes on plateaus of density */
    iftSet    *adjplat;
} iftKnnNode;

/**
 * @brief KNN Graph composed by the training samples of a dataset
 * @author Adan Echemendia
 * @date May 17, 2016
 * 
 */
typedef struct ift_knngraph {
    /** List of nodes in the graph */
    iftKnnNode *node;
    /** List of path value of the nodes */
    float      *pathval;
    /** List of nodes ordered by its path value */
    int        *ordered_nodes;
    /** Number of nodes of the graph */
    int         nnodes;
    /** Maximum arc weight for each value of k */
    float      *maxarcw;
    /** Maximum number of neighbors possible for each node */
    int         kmax;
    /** Best number of neighbors for each node */
    int         k;
    /** Priority queue */
    iftFHeap   *Q;
    /** Corresponding dataset */
    iftDataSet *Z;
} iftKnnGraph;

/**
 * @brief Generic type to define a cut function
 * @details [long description]
 * 
 * @param graph Input clustered graph
 * @return  A graph cut function value
 */
typedef float (*iftKnnGraphCutFun)(iftKnnGraph *graph);

/**
 * @brief Creates a KNN Graph from the training samples of a dataset
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param Z Input dataset
 * @param kmax Maximum number of neighbors possible for each node
 * @return A KNN Graph 
 * 
 */
iftKnnGraph *iftCreateKnnGraph(iftDataSet *Z, int kmax);

/**
 * @brief Destroys a KNN Graph
 * @param graph KNN Graph to destroy
 */
void iftDestroyKnnGraph(iftKnnGraph **graph);

/**
 * @brief It selects a maximum given number of labeled training
 * samples and removes potential outliers from that training
 * set. These outliers are identified as samples that fall into
 * clusters, as computed for a pdf with many domes, whose roots are
 * from different classes. If all samples are removed from any given
 * class, then the set is considered invalid for training. Otherwise,
 * the remaining training samples are considered valid.
 *
 * @author Alexandre Falcao
 * @date Jul 2, 2016
 * @param Input data set for training sample selection.
 * @param Maximum number of training samples.
 * @return Boolean decision about the vality of the data set for training.
 *
 * @note This function needs to be improved in that the Kmax should be
 * selected according to the problem at hand. If it is too high (e.g.,
 * 5 for some problems) it might eliminate too many informative
 * samples.
 */
bool iftSelectSupTrainSamplesWithNoOutliers(iftDataSet *Z, int max_ntrainsamples);

/**
 * @brief Clusters the samples of a KNN Graph by OPF
 * @details Firstly, it finds the best k value (adjacency parameter) in the KNN Graph minimizing the graph cut function, given as a
 * parameter, in the interval [1, graph->kmax]. After that, computes the clusters with OPF.  
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param graph Input KNN Graph
 * @param iftGraphCutFun Graph cut function
 * @return Resulting number of clusters
 */
int  iftUnsupTrain(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Clusters the samples of a KNN Graph by OPF. The difference with the function iftUnsupTrain is that this uses the new and not-tested method iftFastBestKByKnnGraphCut
 * @details Firstly, it finds the best k value (adjacency parameter) in the KNN Graph minimizing the graph cut function, using the function iftFastBestKByKnnGraphCut. After that, computes the clusters with OPF.
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param graph Input KNN Graph
 * @param iftGraphCutFun Graph cut function
 * @return Resulting number of clusters
 */
int iftFastUnsupTrain(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Clusters the samples of a KNN Graph in a specified number of groups by OPF. 
 * @details The returned value is the number of clusters given as result by OPF. This number is then adjusted to 
 * remain with the specified number of groups.  
 * @author Adan Echemendia
 * @date May 17, 2016
 * @param graph Input KNN Graph
 * @param c Desired number of clusters
 * @return Original number of groups given by OPF. 
 */
int  iftUnsupTrainWithCClusters(iftKnnGraph *graph, int c);

/**
 * @brief Propagates group labels in the samples of a dataset given a clustered graph
 * @details [long description]
 * 
 * @param graph Input clustered graph
 * @param Z Input ataset
 * 
 * @return Number of outliers in the dataset
 */
int iftUnsupClassify(iftKnnGraph *graph, iftDataSet *Z);

/**
 * @brief Propagates group labels in the samples of a dataset given a clustered graph. If a sample of <b>Z</b> is recognized as outlier it stays with the label 0 in the dataset. This outlier treatment is the difference of this function with the function <b>iftUnsupClassify</b>
 * @details [long description]
 * @author Adan Echemendia
 * @param graph Input clustered graph
 * @param Z Input ataset
 *
 * @return Number of outliers in the dataset
 */
int iftUnsupClassifyDataSetLessOutliers(iftKnnGraph *graph, iftDataSet *Z);


/**
 * @brief Propagates group labels of a clustered graph, composed by a "representative sample set", in the samples of a large dataset
 * @details At first, the function propagates the group labels using a clustered graph that contains a little initial set. If this result results in a high number of outliers, that means that the "representative sample set" is not good enough and in the next iteration, some of these outliers will be included in the "representative sample set". The process continues until reach a maximum number of iterations or detect a low number of outliers. Note that at the end, the method is trying to learn the PDF of a large dataset from a little sample set.
 *  
 * @param Z Large input dataset
 * @param kmax_perc The maximum number of neighbors possible for each node will be this percent value of the number of training samples in the dataset 
 * @param iftGraphCutFun Graph cut function
 * @param niterations Maximum number of iterations
 * @return A clustered graph, composed by a representative sample set of a large dataset
 */
iftKnnGraph *iftUnsupLearn(iftDataSet *Z, float kmax_perc, iftKnnGraphCutFun iftGraphCutFun, int niterations);

iftKnnGraph *iftUnsupLearn2(iftDataSet *Z, float kmax_perc, iftKnnGraphCutFun iftGraphCutFun, int niterations);

/**
 * @brief Clusters a big dataset using a divide and conquer strategy and the unsupervised OPF method.
 * @note The samples of each partition are selected randomly without repetition.
 * @param dataset A big dataset to be clustered
 * @param nb_partitions Number of partitions to divide the big dataset
 * @param iftGraphCutFun A graph cut function
 * @param kmax_perc1 k max percent in the first levels
 * @param kmax_perc2 k max percent in the last level
  * @returns The number of clusters found
 */
void iftDivideAndConquerRandomlyUnsupOPF2Levels(iftDataSet *dataset, int nb_partitions,iftKnnGraphCutFun iftGraphCutFun, float kmax_perc1, float kmax_perc2);


/**
 * @brief Clusters a big image using a divide and conquer strategy, tile division and the unsupervised OPF method. Here the hierarchy will have two levels, the first with the image block partition and the second with the roots of the first
 * @note The difference with the function iftDivideAndConquerRandomlyUnsupOPF is that in this one is designed to have less partitions because we are not clustering with all samples in each partition, otherwise in each partition we are training with a set of samples and classifying the rest.
 * @param orig_dataset A big dataset (that references an image) to be clustered
 * @param num_blocks Number of tiles to partition the image
 * @param train_percent Percent of train samples for each partition
 * @param iftGraphCutFun A graph cut function
 * @param k_max_percent The maximum number of neighbors possible for each node will be this percent value of the number of training samples in the each partitioned dataset
 * @returns The number of clusters found
 */
iftKnnGraph *iftClusterImageByDivideAndConquerBlocksUsingOPF(iftDataSet *orig_dataset, int num_blocks,
                                                             float train_percent,
                                                             iftKnnGraphCutFun iftGraphCutFun, float kmax_perc1,
                                                             float kmax_perc2);

/**
 * @brief Clusters a big image using a divide and conquer strategy, tile division and the unsupervised OPF method. Here we are only working in the first level, I means I'm clustering each partition independently and I'm not joining the cluster roots in the second level. In each partition some samples are used as train samples and the others are classified by the method OPFUnsupClassify
 * @param orig_dataset A big dataset (that references an image) to be clustered
 * @param num_blocks Number of tiles to partition the image
 * @param train_percent Percent of train samples for each partition
 * @param iftGraphCutFun A graph cut function
 * @param k_max_percent The maximum number of neighbors possible for each node will be this percent value of the number of training samples in the each partitioned dataset
 * @returns The number of clusters found
 */
int iftClusterImageByDividingInBlocksAndUsingOPF(iftDataSet *orig_dataset, int num_blocks,
                                                 float train_percent, iftKnnGraphCutFun iftGraphCutFun,
                                                 float k_max_percent);

/**
 * @brief Propagates group labels in the samples of a dataset given a supervised clustered graph
 * @details Firstly, it replaces the root labels by the root classes and propagate them in the clustered graph. Afterwards,
 * it propagates these root labels in the remaining samples of the dataset.
 * @param graph Input clustered graph composed by a supervised training set
 * @param Z Input dataset
 * @return Number of outliers in the dataset Z
 */
int  iftSemiSupClassify(iftKnnGraph *graph, iftDataSet *Z);

/**
 * @brief Selects the best features in a dataset by multi-scale parameter search 
 * @details [long description]
 * 
 * @param Z Input dataset
 * @param kmax Maximum number of neighbors possible for each node
 */
void  iftUnsupFeatSelecByMSPS(iftDataSet *Z, int kmax);

/**
 * @brief Computes a normalized cut in a clustered graph
 * @details [long description]
 * 
 * @param graph Input clustered graph
 * @return The normalized cut value
 */
float  iftNormalizedCut(iftKnnGraph *graph);

/**
 * @brief Finds the best k value (adjacency parameter) for a clustered graph given a cut function  
 * @details For each k in the interval [1, graph->kmax], the function computes the PDF of all nodes, groups them by OPF, and 
 * determines the cut function value. It returns the k value that minimizes the referred cut function value.
 * @param graph Graph containing the samples to cluster
 * @param iftGraphCutFun Cut Function
 * @return Optimum k value for the clustered graph
 */
int  iftBestkByKnnGraphCut(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Finds the best k value (adjacency parameter) for a clustered graph given a cut function. In this function we are not searching the best k in all the interval from [0,kmax] but we begin with k=kmax and we are decrementing k til we stop to improve the cut function. See also function iftBestkByKnnGraphCut
 * @param graph Graph containing the samples to cluster
 * @param iftGraphCutFun Cut Function
 * @return Optimum k value for the clustered graph
 */
int iftFastBestKByKnnGraphCut(iftKnnGraph *graph, iftKnnGraphCutFun iftGraphCutFun);

/**
 * @brief Builds a new dataset with the boundary samples of an input clustered dataset
 * @details [long description]
 * 
 * @param graph Input clustered graph
 * @param Z Input clustered dataset
 * 
 * @return A new dataset with the boundary samples of the clustered dataset Z
 */
iftDataSet  *iftGraphBoundaryReduction(iftKnnGraph *graph, iftDataSet *Z);

/**
 * @brief Gets the boundary samples from the dataset contained in a clustered graph
 * @details [long description]
 * 
 * @param graph Input clustered graph
 * @return Set of boundary samples from the dataset contained in a clustered graph
 */
iftSet  *iftGetKnnBoundarySamples(iftKnnGraph *graph);

/**
 * @brief Gets the samples corresponding to the root nodes of a clustered graph
 * @details [long description]
 * 
 * @param graph Input clustered graph
 * @return The sample set corresponding to the root nodes of a clustered graph
 */
iftSet  *iftGetKnnRootSamples(iftKnnGraph *graph);

/**
 * @brief Computes an array containing, in each position, a list of sample indexes of each cluster sorted by their distance to the root node. 
 * @details The array has size Z->nlabels+1 containing, in each position, a list of sample indexes of each cluster sorted by their distance 
 * to the root node (increasing distance). The root node is the first element of each list. This function assumes that each cluster 
 * contains only one root node and labels start at 1.
 * 
 * @param graph Input clustered graph
 * @param nNodes Output array of size Z->nlabels+1 with the number of nodes in each cluster
 * 
 * @return Array that contains in each position a list of sample indexes of each cluster, sorted by their distance to the root node
 */
int  **iftGetNodesPerClusterOrderedByDistanceToRoot(iftKnnGraph *graph, int **nNodes);

/**
 * @brief Computes an array containing, in each position, a list of sample indexes of each cluster sorted by their decreasing weight
 * @details The array has size Z->nlabels+1 containing, in each position, a list of sample indexes of each cluster sorted by their 
 * decreasing weight. The root node is the first element of each list. This function assumes that each cluster contains only one root 
 * node and labels start at 1.
 * 
 * @param graph Input clustered graph
 * @param nNodes Output array of size Z->nlabels+1 with the number of nodes in each cluster
 * 
 * @return Array that contains in each position a list of sample indexes of each cluster, sorted by their decreasing weight
 */
int  **iftGetNodesPerClusterOrderedByWeight(iftKnnGraph *graph, int **nNodes);

/**
 * @brief Propagates, for each cluster of the graph, the true label of the corresponding root  
 * @details [long description]
 * 
 * @param graph Input clustered graph
 */
void  iftPropagateClusterTrueLabels(iftKnnGraph *graph);

void  iftPDFByRange(iftKnnGraph *graph);

int   iftUnsupOPF(iftKnnGraph *graph);

#ifdef __cplusplus
}
#endif

#endif
