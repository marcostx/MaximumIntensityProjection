/**
 * @file
 * @brief Parameter Optimization utilities
 * @author Peixinho
 */

#ifndef IFT_PARAMOPTIMIZER_H
#define IFT_PARAMOPTIMIZER_H

#include "iftDataSet.h"
#include "iftDict.h"
#include "iftFitnessFunc.h"
#include "iftMetrics.h"

typedef double (*iftDictFitnessFunc) (iftDict*, iftDict *params);
typedef double (*iftClassifierFitnessFunc)   (iftDict*, iftDataSet*, iftDict *params);
typedef double (*iftDescriptorFitnessFunc)   (iftDict*, iftFileSet*, iftDict *params);

//TODO: change all the void* problem in a params dictionary

/*
 * @brief Searchs for the best score parameters among all the possible values.
 * @param opt Space search params
 * @param f Fitness function, evaluates the classifier given the specified parameters.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * */
iftDict* iftGridSearch(iftDict *opt, iftDictFitnessFunc f, iftDict *problem);

/*
 * @brief Searchs for the best classifier parameters in a given dataset using the iftGridSearch() method.
 * @param opt Space search params
 * @param f Fitness function, evaluates the classifier given the specified parameters.
 * @param dataset The train and test dataset, with thre proper TRAIN/TEST split.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * */
iftDict * iftGridSearchClassifier(iftDict *opt, iftClassifierFitnessFunc f, iftDataSet *dataset, iftSampler *sampler,
                                  iftDict *problem);

/*
 * @brief Searchs for the best descriptor parameters in a given fileset using the iftGridSearch() method.
 * @param opt Space search params
 * @param f Fitness function, evaluates the descriptor given the specified parameters.
 * @param dataset The fileset to be used in descriptor evaluation.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * */
iftDict * iftGridSearchDescriptor(iftDict *, iftDescriptorFitnessFunc f, iftFileSet *dataset, iftDict *problem);

/*
 * @brief Searchs for the best score parameters using random parameters trials.
 * @param opt Space search params
 * @param f Fitness function, evaluates the classifier given the specified parameters.
 * @ntrials The number of random trials.
 * @param problem Extra parameters to an specific problem
 * @author Peixinho
 * @return The best params
 * */
iftDict* iftRandomSearch(iftDict *opt, iftDictFitnessFunc f, int ntrials, iftDict *problem);

iftDict * iftRandomSearchClassifier(iftDict *, iftClassifierFitnessFunc, iftDataSet *dataset, iftSampler *sampler,
                                    int ntrials, iftDict *);
iftDict * iftRandomSearchDescriptor(iftDict *, iftDescriptorFitnessFunc fitnessFunc, iftFileSet *, int ntrials,
                                    iftDict *);


#endif //IFT_PARAMOPTIMIZER_H
