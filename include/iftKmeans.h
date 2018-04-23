#ifndef IFT_KMEANS_H_
#define IFT_KMEANS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftDataSet.h"

void iftSimpleKmeansRun(iftDataSet* Z, iftDataSet** pZk, int maxIterations,float minImprovement);

/**
 * @brief Cluster a dataset by KMeans method
 * @author Adan Echemendia
 * @param Z A dataset to be clustered
 * @param k Desired number of clusters
 * @param maxIterations Maximun Number of iterations
 * @param minImprovement Threshold of error difference to consider between the iterations
 * @param kmedoids 1 if we want to do kmedoids, 0 otherwise
 */
void iftClusterDataSetByKMeans(iftDataSet *Z, int nb_cluster, int maxIterations, float minImprovement,char kmedoids, char grid_sampling);

void  iftKmeansRun(int bKmedoids,iftDataSet* Z,iftDataSet** pZk,int maxIterations,float minImprovement);
void  iftSphericalKmeansRun(iftDataSet* Z,iftDataSet** pZk,int maxIterations);

/**
 * @brief Initializes k centroid getting random samples from dataset Z
 * @param Z A dataset
 * @param k Number of centroids
 * @return A new dataset with the centroids
 */
iftDataSet* iftKmeansInitCentroidsFromSamples(iftDataSet* Z, int k);

iftDataSet* iftKmeansInitCentroidsRandomNormal(iftDataSet* Z, int k);

iftDataSet *iftNearestSamplesFromCentroids(iftDataSet* Zorig, iftDataSet *Zc, iftDataSet* Zk, int m, int include_centroids);
iftDataSet *iftBuildDataSetAsNearestSamplesFromCentroids(iftDataSet* Zorig, iftDataSet *Zc,
			iftDataSet* Zk, int include_centroids, char *method);

#ifdef __cplusplus
}
#endif

#endif
