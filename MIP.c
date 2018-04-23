#include <stdio.h>
#include "ift.h"

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


int ComputeIntersection(Matrix *Tpo, Image *img, Matrix *Tn, VolumeFaces *vf, int *p1, int *pn)
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


int LinearInterpolationValue(Image *img, FVoxel v)
{
    int p[8], i, value;

    return value;
}


Image *MaximumIntensityProjection(Image *img, float xtheta, float ytheta, float ztheta)
{
    iftImage *output = CreateImage(Nu, Nv, 1);

    return output;
}
