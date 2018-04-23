#ifndef _IFT_PARSE_INPUT_H_
#define _IFT_PARSE_INPUT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftFile.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/*
    This file has the goal to implement basic parsing methods, such as
    split, endswith, iftCreateAndLoadImageNames.
    
    It is expected that the image names are: XXXX_YYYY.Extension. 
    Where XXXX is and integer of the class and YYYY the imagem unique identifier.

    The iftImageNames structure is inteded to be created and loaded in the iftCreateAndLoadImageNames method.

    The following code segment, open a directory, count the number of images with a given extension,
    load all the names in the structure, show the usage of the scruture.


    Usage suggestion:
    char * directory_path = argv[1];
    int i=0, number_of_images = iftCountImageNames(directory_path, "ppm");

    iftImageNames * IN;
    iftImageNames * GT;

    IN = iftCreateAndLoadImageNames(number_of_images, directory_path, "ppm");
    GT = iftCreateAndLoadImageNames(number_of_images, directory_path, "pgm");

    for (i = 0; i < number_of_images ; i++)
        iftReadImageP6(IN[i].image_name);
        ....

    iftDestroyImageNames(IN);
    iftDestroyImageNames(GT);
*/


typedef struct ift_iftImageNames {
  char image_name[256];
  int attribute;
} iftImageNames;


iftImageNames * iftCreateImageNames (int nimages);
void iftDestroyImageNames (iftImageNames * IN_structure);

int  iftCountImageNames (char * directory_path, char * file_extension);
/* Allocates the iftImageNames structure. Use iftCountImageNames to count */
iftImageNames *iftCreateAndLoadImageNames (int nimages, char * directory_path, char *file_extension);

/* This function expects that the image names are padronized as follows: XXXX_YYYY.extension. Where XXXX is the class number and YYYY is the name. */
/* If there are less then the specified number of samples, only the number of images available will be loaded */
iftImageNames * iftCreateAndLoadImageNamesSample (char * directory_path, char *file_extension, int number_of_classes, int number_of_samples_per_class);
int iftCountNumberOfClasses(iftImageNames *IN_structure, int nimages);

char * iftSplitStringOld (char *string, char * separator, int position);

/* Directory operations */
void iftCreateDirectory(char *path);
void iftCreateOrCleanDirectory(char *path);
void iftCreateOrCleanFilesInDirectory(char *path, char *token1, char *token2);


#ifdef __cplusplus
}
#endif

#endif
