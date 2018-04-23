/**
 * @file iftDeepLearning.h
 * @brief Deep Learning structures and algorithms
 * @author Peixinho
 * @date Jun, 2016
 */


#ifndef IFT_DEEPLEARNING_H_
#define IFT_DEEPLEARNING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftGraphics.h"
#include "iftColor.h"
#include "iftFiltering.h"
#include "iftInterpolation.h"
#include "iftClustering.h"
#include "iftSVM.h"
#include "iftKmeans.h"
#ifdef IFT_GPU
#include "iftDeepLearning.cuh"
#endif

#define ACTIV_THRES_LIMIT  1000000
#define ALPHA_LIMIT        100.

typedef struct ift_neural_network iftNeuralNetwork;
typedef struct ift_neural_layer iftNeuralLayer;

typedef void (*iftNeuralForward) (iftNeuralNetwork*, int);
typedef void (*iftNeuralBackward) (iftNeuralNetwork*, int);
typedef void (*iftNeuralUpdate) (iftNeuralNetwork*, int);

typedef float (*iftActivationFunction) (float);

typedef enum ift_neural_layer_type { FC_LAYER, CONV_LAYER } iftNeuralLayerType;

struct ift_neural_layer {

	iftDict* imgs;

    iftTensor* weight;
    iftTensor* delta;
    iftTensor* deltaW;
    iftTensor* bias;

    iftTensor* data;
    iftTensor* error;

    int idxCounter;
    iftIntArray* idx;
	iftFileSet* fileset;
    iftDataSet* dataset;
	iftSampleStatus status;
	int nsamples;

	iftTensor* mean;
	iftTensor* std;

    iftActivationFunction act;
    iftActivationFunction actDeriv;

    iftNeuralForward forward;
    iftNeuralBackward backward;
	iftNeuralUpdate update;

	iftNeuralLayerType type;

	int nclasses;

    int ninput, noutput;//nbands for convolution
    int kernelx, kernely, kernelz;
    int xsize, ysize, zsize;//convolutional only
};

struct ift_neural_network {
    int nlayers;
    iftNeuralLayer** layers;

    float learnRate;
	float weightDecay;

    int miniBatch;

    float error;

    int epoch;

    bool verbose;

    iftTensor* label;
};

void iftPrintTensor(iftTensor* t);

void iftNetForward(iftNeuralNetwork* net);

void iftNetBackward(iftNeuralNetwork* net);

void iftNetUpdateWeight(iftNeuralNetwork *net);

void iftNeuralNetGradientDescent(iftNeuralNetwork* net, int epochs);

iftNeuralLayer* iftCreateNeuralLayer();

  void iftDestroyNeuralLayer(iftNeuralLayer** l);

iftNeuralNetwork* iftCreateNeuralNetwork(int nlayers);
  
void iftDestroyNeuralNetwork(iftNeuralNetwork** net);

iftNeuralLayer* iftDiscriminantLossLayer(int minibatch, int ninput);

void iftDiscriminantLossBackward(iftNeuralNetwork* net, int l);

void iftImageForward(iftNeuralNetwork* net, int l);

void iftDataSetForward(iftNeuralNetwork* net, int l);

void iftNoForward(iftNeuralNetwork* net, int l);

void iftNoBackward(iftNeuralNetwork* net, int l);

void iftNoUpdate(iftNeuralNetwork* net, int l);

void iftLossBackward(iftNeuralNetwork* net, int l);

void iftFullyConnectedUpdate(iftNeuralNetwork* net, int l);

void iftConvolutionUpdate(iftNeuralNetwork* net, int l);

void iftFullyConnectedForward(iftNeuralNetwork* net, int l);

void iftFullyConnectedBackward(iftNeuralNetwork* net, int l);

iftNeuralLayer* iftLossLayer(int minibatch, int ninput);

iftNeuralLayer* iftConvolutionLayer2D(int minibatch, int xsize, int ysize, int nbands, int kernelXSize, int kernelYSize,
                                      int nkernels, iftActivationFunction act, iftActivationFunction deriv);

iftNeuralLayer* iftImageLayer(iftNeuralNetwork* net, const iftFileSet* fileset);

iftNeuralLayer* iftDataSetLayer(iftNeuralNetwork*, const iftDataSet* Z);

iftNeuralLayer* iftFullyConnectedLayer(int minibatch, int ninput, int noutput, iftActivationFunction act, iftActivationFunction deriv);

iftTensor* iftNeuralNetClassifyImageSet(iftNeuralNetwork* net, iftFileSet* fileSet, uchar status);

void iftNeuralNetClassifyDataSet(iftNeuralNetwork* net, iftDataSet* Z, uchar status);

iftDataSet* iftNeuralNetDataSetFeatures(iftNeuralNetwork* net, int l);

void iftNeuralCorrelation(iftTensor* data, iftTensor* kernel, iftTensor* out);
void iftNeuralConvolution(iftTensor* data, iftTensor* kernel, iftTensor* out);
void iftNeuralConvolutionUpdate(iftTensor* error, iftTensor* d_out, iftTensor* out);

float iftRelu(float f);

float iftReluDeriv(float f);

float iftSigm(float p);

float iftSigmDeriv(float p);

float iftTanh(float f);

float iftTanhDeriv(float f);

/**
 * @brief Convolution Network Descriptor
 * Contains the ConvNet architecture information along with the kernel weights.
 * @author Peixinho
 * @date Jun, 2016
 **/
typedef struct ift_convnetwork {
    int nlayers; // number of layers
    int nstages; // number of stages is 3*layers + 2
    iftAdjRel *input_norm_adj; // normalization adjacency of the input image
    int input_norm_adj_param; // parameter to create the normalization adjacency of the input image
    int input_xsize, input_ysize, input_zsize, input_nbands; // input image dimensions
    iftMMKernel **k_bank; // one kernel bank per layer
    int *k_bank_adj_param; // parameter to create adjacency of each kernel bank
    int *nkernels; // number of kernels per layer
    iftAdjRel **pooling_adj; // one pooling adjacency per layer
    int *pooling_adj_param; // parameters to create one pooling adjacency per layer
    int *stride; // one pooling stride per layer
    float *alpha; // one parameter of the pooling metric per layer
    iftAdjRel **norm_adj; // one normalization adjacency for the end of each layer
    int *norm_adj_param; // parameter to create one normalization adjacency for the end of each layer
    iftMatrix **img_ind;  // one image index matrix per layer for fast filtering
    int rescale;  // Output the image of the last layer with its 0- actual resolution or 1- with the resolution of the input image
    int with_weights; // Write 1-with or 0-without kernel weights */
    float *activ_thres; // one activation threshold per layer
    int *activ_option; // one option per layer
} iftConvNetwork;

/**
 * @brief Multi Scale Convolution Network.
 * Contains many ConvNets, analysing the image in different input scales.
 * @date Jun, 2016
 * @author Peixinho
 */
typedef struct ift_msconvnetwork {
    int nscales; // number of scales
    float *reduction_factor; // one reduction factor per scale
    iftAdjRel *output_norm_adj;  // normalization adjacency of the output image
    float output_norm_adj_param; // parameter to create the normalization adjacency of the output image
    iftConvNetwork **convnet; // one convolution network per scale
} iftMSConvNetwork;

/**
 * @brief Creates a Convolution Network with nlayers.
 * @param nlayers Number of layers in the ConvNet.
 * @return The Convolution Network structure.
 * @author Peixinho
 */
iftConvNetwork *iftCreateConvNetwork(int nlayers);

/**
 * @brief Destroys the Convolution Network.
 * @param convnet
 * @author Peixinho
 * @date Jun, 2016
 */
void iftDestroyConvNetwork(iftConvNetwork **convnet);

/**
 * @brief Create Multi Scale Convolution Network
 * @param nscales Number of input scales
 * @return Multi Scale Convolution Network
 * @author Peixinho
 * @date Junt, 2016
 */
iftMSConvNetwork *iftCreateMSConvNetwork(int nscales);

/**
 * @brief Destroy Multi Scale Convolution Network
 * @param msconvnet
 * @author Peixinho
 * @date Jun, 2016
 */
void iftDestroyMSConvNetwork(iftMSConvNetwork **msconvnet);

/**
 * @brief Create Convolution Network from a dictionary of parameters
 * @param params Dictionary of parameters
 * @return Convolution Network
 * @author Peixinho
 * @date Nov, 2016
 */
iftConvNetwork* iftCreateConvNetworkFromDict(const iftDict* params);

/**
 * @brief Read Convolution Network from disk.
 * The Convolution Network can be defined as a .convnet or .json (Describes the )
 * @param filename Convolution Network definition.
 * @return
 */
iftConvNetwork *iftReadConvNetwork(const char *filename);

/**
 * @brief Print the Convolution Network architecture.
 * @param convNet Convolution Network
 * @author Peixinho
 * @date Jun, 2016
 */
void iftPrintConvNetArch(iftConvNetwork *convnet);

/**
 * @brief Check if the architecture is feasible.
 * The pooling and convolution layers may cause the image to have the size reduced.
 * This funcion checks if all the hidden layers will have a valid size (greater than zero), and also checks if the ammount of memory needed to represent this network can be allocated.
 * @param convnet
 * @return True if the architecture is feasible, False otherwise
 */
bool iftValidArchitecture(iftConvNetwork *convnet);

/**
 * @brief Saves the Convolution Network in disk.
 * The ConvNet can be saved in two formats (.convnet/.json)
 * @warning .json for now is only able to save the ConvNet architecture, to save the ConvNet weights you must use .convnet format
 * @param convnet Convolution Network
 * @param filename Filepath
 */
void iftWriteConvNetwork(iftConvNetwork *convnet, const char *filename);

/**
 * @brief Appends two Convolution Networks, as one single serialized Convolution Network.
 * The second ConvNet is placed at the end of the first, taking as input the first ConvNet output.
 * @param convnet1 First ConvNet
 * @param convnet2 Second ConvNet
 * @return Joined ConvNet.
 * @author Samuel
 * @date Jun, 2016
 */
iftConvNetwork *iftMergeConvNetwork(iftConvNetwork *convnet1, iftConvNetwork *convnet2);

/**
 * @brief Read a Multi Scale Convolution Network from disk.
 * @param filename Filepath
 * @return The ConvNet.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMSConvNetwork *iftReadMSConvNetwork(char *filename);

/**
 * @brief Stores a Multi Scale Convolution Network on disk.
 * @param msconvnet Multi Scale ConvNet
 * @param filename Filepath
 * @author Peixinho
 * @date Jun, 2016
 */
void iftWriteMSConvNetwork(iftMSConvNetwork *msconvnet, const char *filename);

/**
 * @brief Apply ConvNet in image.
 * The ConvNet output is a multiband image
 * @warning The high dimensional output is better suited for linear classifier. See also ::iftSVM and ::iftLogReg.
 * @param img Input image.
 * @param convnet Convolution Network
 * @return The multibandimage.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMImage *iftApplyConvNetwork(iftMImage *img, iftConvNetwork *convnet);


/**
 * @brief Apply Multi Scale ConvNet in image.
 * The ConvNet output is a multiband image
 * @warning The high dimensional output is better suited for linear classifier. See also ::iftSVM and ::iftLogReg.
 * @param img Input image.
 * @param convnet Multi Scale Convolution Network
 * @return The multibandimage.
 * @date Jun, 2016
 * @author Peixinho
 */
iftMImage *iftApplyMSConvNetwork(iftMImage *img, iftMSConvNetwork *convnet);

/**
 * @brief Computes the multiband image dimensions along the hidden layers.
 * @param convnet Convolution Network
 * @param xsize
 * @param ysize
 * @param zsize
 * @param nbands
 * @return True if the architecture is valid, False otherwise. See also iftValidArchitecture()
 */
bool iftImageDimensionsAlongNetwork(iftConvNetwork *convnet, int *xsize, int *ysize, int *zsize, int *nbands);

/**
 * @brief Learns the kernel weights of a one layer Convolution Network.
 * @warning Only works with 1 layer ConvNets.
 * @param img Input image
 * @param convnet Convolution Network
 * @param nsamples Number of sampled patches for clustering
 * @param kmax_perc OPF unsup maximum k
 * @param whitening Apply, or not, whitening in the dataset
 * @author Peixinho
 * @date Jun, 2016
 */
void iftUnsupLearnKernels(iftMImage *img, iftConvNetwork *convnet, int nsamples, float kmax_perc, bool whitening);

/**
 * Learns the kernel weights of a one layer Convolution Network with KMeans clustering.
 * @warning Only works with 1 layer ConvNets.
 * @param img Input image
 * @param convnet Convolution Network
 * @param nsamples Number of sampled patches for clustering
 * @param k Number of clusters in KMeans
 * @param whitening Apply, or not, whitening in the dataset
 * @date Jun, 2016
 * @author Peixinho
 */
void iftUnsupLearnKernelsByKmeans(iftMImage *img, iftConvNetwork *convnet, int nsamples, int k, char whitening);

/**
 * Learns the kernel weights of a one layer Convolution Network with KMedoids clustering.
 * @warning Only works with 1 layer ConvNets.
 * @param img Input image
 * @param convnet Convolution Network
 * @param nsamples Number of sampled patches for clustering
 * @param k Number of clusters in KMedoids
 * @param whitening Apply, or not, whitening in the dataset
 * @author Peixinho
 * @date Jun, 2016
 */
void iftUnsupLearnKernelsByKmedoids(iftMImage *img, iftConvNetwork *convnet, int nsamples, int k, char whitening);

/**
 * @brief Learns the kernel weights of a target layer Convolution Network with OPF Clustering over patches of training images.
 * @param directory_path Directory containing the training images
 * @param convnet Convolution Network
 * @param target_layer Target layer to learn kernel weights
 * @param nsamples Number of patch samples used in clustering
 * @param kmax_perc Maximum k percentage for OPF clustering
 * @param whitening Apply, or not, whitening in the dataset
 * @date Jun, 2016
 * @author Peixinho
 */
void iftUnsupLearnKernelsFromImages(char *directory_path, iftConvNetwork *convnet, int target_layer, int nsamples,
                                    float kmax_perc, char whitening);
/**
 * @brief Learns the kernel weights of a target layer Convolution Network with KMeans Clustering over patches of training images.
 * @param directory_path Directory containing the training images
 * @param convnet Convolution Network
 * @param target_layer Target layer to learn kernel weights
 * @param nsamples Number of patch samples used in clustering
 * @param k Number of clusters in KMeans
 * @param whitening Apply, or not, whitening in the dataset
 * @date Jun, 2016
 * @author Peixinho
 */
void iftUnsupLearnKernelsByKmeansFromImages(char *directory_path, iftConvNetwork *convnet, int target_layer,
                                            int nsamples, int k, char whitening);

/**
 * @brief Learns the kernel weights of a target layer Convolution Network with KMedoids Clustering over patches of training images.
 * @param directory_path Directory containing the training images
 * @param convnet Convolution Network
 * @param target_layer Target layer to learn kernel weights
 * @param nsamples Number of patch samples used in clustering
 * @param k Number of clusters in KMedoids
 * @param whitening Apply, or not, whitening in the dataset
 * @date Jun, 2016
 * @author Peixinho
 */
void iftUnsupLearnKernelsByKmedoidsFromImages(char *directory_path, iftConvNetwork *convnet, int target_layer,
                                              int nsamples, int k, char whitening);

/**
 * @brief Learns the kernel weights of a target layer Convolution Network with Spherical KMeans Clustering over patches of training images.
 * @param directory_path Directory containing the training images
 * @param convnet Convolution Network
 * @param target_layer Target layer to learn kernel weights
 * @param nsamples Number of patch samples used in clustering
 * @param k Number of clusters in Spherical KMeans
 * @param whitening Apply, or not, whitening in the dataset
 * @date Jun, 2016
 * @author Peixinho
 */
void iftUnsupLearnKernelsBySphericalKmeansFromImages(char *directory_path, iftConvNetwork *convnet, int target_layer,
                                              int nsamples, int k, char whitening);

/**
 * Learns a multiband kernel from a dataset of image patches with unsupervised OPF clustering.
 * @param Z Patches dataset
 * @param A Patches adjacency
 * @param kmax_perc Maximum k percentage for OPF clustering
 * @param whitening Apply, or not, whitening in the dataset
 * @return The multiband kernel.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMMKernel *iftUnsupLearnKernelsFromDataSet(iftDataSet *Z, iftAdjRel *A, float kmax_perc, bool whitening);

/**
 * Learns a multiband kernel from dataset of image patches with unsupervised KMedoids clustering.
 * @param Z Patches dataset
 * @param A Patches adjacency
 * @param k Number of clusters for KMedoids clustering
 * @param whitening Apply, or not, whitening in the dataset
 * @return The multiband kernel.
 * @author Peixinho
 * @date Jun, 2016
 */
iftMMKernel *iftUnsupLearnKernelsByKmedoidsFromDataSet(iftDataSet *Z, iftAdjRel *A, int k, char whitening);

/**
 * @brief Learns the multiband kernel from a dataset of image patches with Suport Vector Machine.
 * @param dataset_names List of image paths to be used in the learning stage.
 * @param nclasses Number of classes
 * @param nkernels Number of kernels
 * @param A Kernel adjacency
 * @param nbands Input image number of bands
 * @param nsamples Number of patch samples used in kernel learning.
 * @param[out] predict_mean TODO: discover what this variable returns
 * @return The multiband kernels
 */
iftMMKernel **iftSupLearnKernelBanksFromDataSets(iftImageNames *dataset_names, int nclasses,
                                                 int nkernels, iftAdjRel *A, int nbands, int nsamples,
                                                 float *predict_mean);
/**
 * @brief Supervised Learn dataset with SVM
 * @param Z Dataset
 * @param[out] predict_sum TODO: discover what this variable returns
 * @return The hyperplanes learned
 * @author Peixinho
 * @date Jun, 2016
 */
iftSVMHyperplane **iftSupLearnKernelsFromDataSet(iftDataSet *Z, float *predict_sum);

/**
 * @brief Randomly selects patches from image as kernel.
 * @param directory_path Directory path containing images
 * @param convnet Convolution Network
 * @param target_layer Target layer to be learned
 * @param whitening To apply, or not, whitening in dataset
 * @param nsamples Number of samples extracted from image.
 * @date Jun, 2016
 * @author Peixinho
 */
void iftSelectRandomKernelsFromImages(char *directory_path, iftConvNetwork *convnet, int target_layer, char whitening,
                                      int npatches);

/**
 * @brief Randonly selects samples from a directory of images. The whitened selected samples are used as kernel weights for a target layer.
 * @author Peixinho
 * @date Jun, 2016
 * @param directory_path Directory path containing images
 * @param convnet Convolution Network
 * @param target_layer Target layer to learn weights
 * @param nsamples Number of sample patches
 */
void iftWhiteningFromImages(char *directory_path, iftConvNetwork *convnet, int target_layer, int nsamples);

/**
 * @brief Learn weights for a one layered Convolution Network, with randomly selected patchs from image
 * @param img Input image
 * @param convnet Convolution Network
 * @param whitening Apply, or not, whitening in dataset
 */
void iftSelectRandomKernels(iftMImage *img, iftConvNetwork *convnet, char whitening);

/**
 * @brief Creates the adjacency objects, according to kernel size in convolution layers
 * @param convnet Convolution Network
 */
void iftCreateAdjRelAlongNetwork(iftConvNetwork *convnet);

/**
 * @brief Creates the fast convolution matrices for each layer
 * @param convnet Convolution Network
 */
void iftMImageIndexMatrices(iftConvNetwork *convnet);

/**
 * @brief Creates randomly initialized kernels
 * @param convnet Convolution Network
 * @author Peixinho
 * @date Jun, 2016
 */
void iftCreateRandomKernelBanks(iftConvNetwork *convnet);

/**
 * @brief Creates a ConvNet from the first layers of a ConvNet
 * @param convnet Convolution Network
 * @param nlayers Number of layers to be copied
 * @author Peixinho
 * @date Jun, 2016
 * @return
 */
iftConvNetwork *iftCopyFirstLayersConvNetwork(iftConvNetwork *convnet, int nlayers);

/**
 * @brief Loads kernel weights from file
 * @param convnet Convolution Network
 * @param fp Pointer for file with kernel weights
 * TODO: I believe this shouldbe a private function
 * @author Peixinho
 * @date Jun, 2016
 */
void iftLoadKernelBanks(iftConvNetwork *convnet, FILE *fp);

/**
 * @brief Converts SVM hyperplanes to kernel
 * @param h SVM hyperplanes
 * @param kernel Output kernel
 * @param band_size Input image bands
 * @author Peixinho
 * @date Jun, 2016
 */
void iftConvertHyperplaneInKernelBand(iftSVMHyperplane *h, iftBand *kernel, int band_size);

/**
 * @brief Prints Convolution Network architecture
 * @param convNet Convolution Network
 * @author Peixinho
 * @date Jun, 2016
 */
void iftPrintConvNetArch(iftConvNetwork *convNet);

#ifdef __cplusplus
}
#endif

#endif
