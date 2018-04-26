#include <stdio.h>
#include "ift.h"



Voxel GetVoxelCoord(iftImage *img, int p)
{
    iftVoxel u;

    u.x = ittGetXCoord(img, p);
    u.y = iftGetYCoord(img, p);
    u.z = iftGetZCoord(img, p);

    return (u);
}

/* this function creas the rotation/translation matrix for the given theta */
iftMatrix *createTransformationMatrix(iftImage *img, int xtheta, int ytheta)
{
    iftMatrix *resMatrix = NULL;

    iftVector v1 = {.x = (float)img->xsize / 2.0, .y = (float)img->ysize / 2.0, .z = (float)img->zsize / 2.0};
    iftMatrix *transMatrix1 = iftTranslationMatrix(v1);

    iftMatrix *xRotMatrix = iftRotationMatrix(IFT_AXIS_X, -xtheta);
    iftMatrix *yRotMatrix = iftRotationMatrix(IFT_AXIS_Y, -ytheta);

    float D = sqrt(img->xsize*img->xsize + img->ysize*img->ysize);
    iftVector v2 = {.x = -(D / 2.0), .y = -(D / 2.0), .z = -(D / 2.0)};
    iftMatrix *transMatrix2 = iftTranslationMatrix(v2);


    resMatrix = iftMultMatricesChain(4, transMatrix1, xRotMatrix,yRotMatrix, transMatrix2);

    return resMatrix;
}

int DDA(iftImage *img, iftVoxel p1, iftVoxel pn)
{
    int n, k;
    //iftVoxel p;
    float px, py;
    float J=0;
    int Dx,Dy;
    float dx=0,dy=0;
    //iftCreateMatrix* J;


    if (p1.x == pn.x && p1.y == pn.y)
        n=1;
    else{
        Dx=pn.x - p1.x;
        Dy=pn.y - p1.y;

        if( abs(Dx) >= abs(Dy) ){
            n = abs(Dx)+1;
            dx = sign(Dx);
            dy = (dx * Dy)/Dx;
        }
        else{
            n = abs(Dy)+1;
            dy = sign(Dy);
            dx = (dy * Dx)/Dy;
        }
    }

    px = p1.x;
    py = p1.y;


    // TODO: calcular I como interpolacao

    for (k = 1; k < n; k++)
    {
        //J+=  (float)LinearInterpolationValue(img, px, py);
        J+=  iftImgVal2D(img, (int)px, (int)py);
        //J+=  (float)LinearInterpolationValue(img, px, py);

        px = px + dx;
        py = py + dy;
    }

    return (int)J;
}


int ComputeIntersection(iftMatrix *Tpo, iftImage *img, iftMatrix *Tn, iftVolumeFaces *vf, int *p1, int *pn)
{
    return;
}


int isValidPoint(iftImage *img, iftVoxel u)
{
    //printf("u = %d, %d\n",u.x,u.y );
    //printf("img size = %d, %d\n",img->xsize,img->ysize );
    if ((u.x >= 0) && (u.x < img->xsize) &&
        (u.y >= 0) && (u.y < img->ysize)){
        return 1;
    }
    else{
        return 0;
    }
}


int LinearInterpolationValue(iftImage *img, iftVoxel v)
{
    int p[8], i, value;

    return value;
}


Image *MaximumIntensityProjection(iftImage *img, float xtheta, float ytheta, float ztheta)
{
    int diagonal = 0;
    double maxIntensity;
    int p=0;
    int Nu, Nv;
    int p1, pn;
    float* lineVals;
    iftVoxel p0, v1, vn;
    //
    iftMatrix *Mt1, *Mt2, *Mrx, *Mry, *Mrz, *Mtemp, *T;
    iftMatrix *Norigin, *Nend, *Tnorigin, *Tnend, *Tn;
    iftMatrix *Tpo;

    diagonal = int(sqrt((double) (img->xsize * img->xsize) + (img->ysize * img->ysize) + (img->zsize * img->zsize)));
    Nu = Nv = diagonal;
    iftImage *output = iftCreateImage(Nu, Nv, 1);
    T  = createTransformationMatrix(img, xtheta, ytheta);

    Norigin =  iftCreateMatrix(4, 1);
    iftMatrixElem(Norigin, 0, 0) = 0;
    iftMatrixElem(Norigin, 0, 1) = 0;
    iftMatrixElem(Norigin, 0, 2) = 0;
    iftMatrixElem(Norigin, 0, 3) = 1;
    Tnorigin = MatrixMultiply(T, Norigin);


    for (p = 0; p < output->n; p++)
    {
        p1 = pn = -1;
        p0 = GetVoxelCoord(output, p);
        Mtemp = iftVoxelToMatrix(p0);
        Tpo = iftMultMatrices(T, Mtemp);


        if (ComputeIntersection(Tpo, img, Tnorigin, vf, &p1, &pn))
        {
            v1 = GetVoxelCoord(img, p1);
            vn = GetVoxelCoord(img, pn);

            maxIntensity = IntensityProfile(img, v1, vn);

            output->val[p] = maxIntensity;
        }
        else
            output->val[p] = 0;

        DestroyMatrix(Mtemp);
        DestroyMatrix(Tpo);
    }

    DestroyMatrix(Mt1);
    DestroyMatrix(Mt2);
    DestroyMatrix(Mrx);
    DestroyMatrix(Mry);
    DestroyMatrix(Mrz);
    DestroyMatrix(T);
    DestroyMatrix(Norigin);
    DestroyMatrix(Nend);
    DestroyMatrix(Tnorigin);
    DestroyMatrix(Tnend);
    DestroyMatrix(Tn);
    DestroyVolumeFaces(vf);

    return output;
}
