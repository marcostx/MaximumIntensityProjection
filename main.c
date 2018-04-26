#include <stdio.h>
#include "ift.h"


int main(int argc, char *argv[])
{
    if (argc != 6)
        Error("Run: ./main <filename> <output> <xtheta> <ytheta> <ztheta>", "main");

    char buffer[512];

    float tx, ty, tz;
    tx = atof(argv[3]);
    ty = atof(argv[4]);
    //tz = atof(argv[5]);
    iftImage *img = iftReadExtImage(argv[1]);
    iftImage *output = NULL;

    output = MaximumIntensityProjection(img, tx, ty);
    //sprintf(buffer, "data/%.1f%.1f%.1f%s", tx, ty, tz, argv[2]);

    //WriteImageP2(output, buffer);
    iftDestroyImage(img);
    iftDestroyImage(output);
    return 0;
}
