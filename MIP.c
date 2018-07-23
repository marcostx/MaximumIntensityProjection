#include <stdio.h>
#include "ift.h"

#define GetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define GetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define GetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))

#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
#define GetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])

int VoxelValue(iftImage *img, iftVoxel v)
{
    return img->val[GetVoxelIndex(img, v)];
}

int sign( int x ){
    if(x >= 0)
        return 1;
    return -1;
}

typedef struct volume
{
    iftMatrix* orthogonal;
    iftMatrix* center;
}iftVolumeFaces;

int isValidPoint(iftImage *img, iftVoxel u)
{
    if ((u.x >= 0) && (u.x < img->xsize) &&
        (u.y >= 0) && (u.y < img->ysize) &&
        (u.z >= 0) && (u.z < img->zsize)){
        return 1;
    }
    else{
        return 0;
    }
}

float MatrixInnerProduct(iftMatrix *A, iftMatrix *B)
{
  float result = 0;

  for (int i = 0; i < A->ncols; i++)
    result += (A->val[i] * B->val[i]);

  return result;
}

iftVoxel GetVoxelCoord(iftImage *img, int p)
{
    iftVoxel u;

    u.x = GetXCoord(img, p);
    u.y = GetYCoord(img, p);
    u.z = GetZCoord(img, p);

    return u;
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
    iftVoxel p;
    float J=0, max=0;
    int Dx,Dy,Dz;
    float dx=0,dy=0,dz=0;
    //iftCreateMatrix* J;


    if (p1.x == pn.x && p1.y == pn.y && p1.z == pn.z)
        n=1;
    else{
        Dx=pn.x - p1.x;
        Dy=pn.y - p1.y;
        Dz=pn.z - p1.z;

        if( abs(Dx) >= abs(Dy) && abs(Dx) >= abs(Dz) ){
            n = abs(Dx)+1;
            dx = sign(Dx);
            dy = (dx * Dy)/Dx;
            dz = (dx * Dz)/Dx;
        }
        else if ( abs(Dy) >= abs(Dx) && abs(Dy) >= abs(Dz)){
            n = abs(Dy)+1;
            dy = sign(Dy);
            dx = (dy * Dx)/Dy;
            dz = (dy * Dz)/Dy;
        }
        else{
            n = abs(Dz)+1;
            dz = sign(Dz);
            dx = (dz * Dx)/Dz;
            dy = (dz * Dy)/Dz;
        }
    }

    p.x = p1.x;
    p.y = p1.y;
    p.z = p1.z;


    for (k = 1; k < n; k++)
    {
        //printf("debug\n");
        //J+=  (float)LinearInterpolationValue(img, px, py);
        //J+=  iftImgVal2D(img, (int)px, (int)py);
        if (isValidPoint(img,p)){
          //J= iftImgVal(img,p.x,p.y,p.z);

          iftPoint aux;
          aux.x = p.x;
          aux.y = p.y;
          aux.z = p.z;

          // pegando o ponto com interpolacao
          J= iftImageValueAtPoint(img,aux);

          if (J>max)
            max=J;
        }

        p.x = p.x + dx;
        p.y = p.y + dy;
        p.z = p.z + dz;
    }

    return (int)max;
}




iftVector columnVectorMatrixToVector(iftMatrix *m)
{
    iftVector v = {.x = iftMatrixElem(m, 0, 0), .y = iftMatrixElem(m, 0, 1), .z = iftMatrixElem(m, 0, 2)};
    v.x = iftAlmostZero(v.x) ? 0.0 : v.x;
    v.y = iftAlmostZero(v.y) ? 0.0 : v.y;
    v.z = iftAlmostZero(v.z) ? 0.0 : v.z;

    return v;
}


int ComputeIntersection(iftMatrix *Tpo, iftImage *img, iftMatrix *Tn, iftVolumeFaces *vf, iftVoxel *p1, iftVoxel *pn)
{

    float lambda[6] = { -1};
    float max=-9999999.9, min=9999999.9;
    int i;
    p1->x=pn->x=p1->y=pn->y,p1->z=pn->z=-1;
    iftMatrix *Nj = iftCreateMatrix(1, 3);
    iftMatrix *DiffCandP0 = iftCreateMatrix(1, 3);
    float NdotNj = 0, DiffShiftDotNj = 0;
    iftVector v1, v2, v3;
    iftVoxel v;


    for (i = 0; i < 6; i++) {
      iftMatrixElem(Nj, 0, 0) = vf[i].orthogonal->val[0];
      iftMatrixElem(Nj, 0, 1) = vf[i].orthogonal->val[1];
      iftMatrixElem(Nj, 0, 2) = vf[i].orthogonal->val[2];
      v1 = columnVectorMatrixToVector(Nj);
      v2 = columnVectorMatrixToVector(Tn);

      NdotNj = iftVectorInnerProduct(v1,v2);

      //NdotNj = MatrixInnerProduct(Tn,Nj);

      if (NdotNj != 0){

        iftMatrixElem(DiffCandP0, 0, 0) = vf[i].center->val[0] - Tpo->val[0];
        iftMatrixElem(DiffCandP0, 0, 1) = vf[i].center->val[1] - Tpo->val[1];
        iftMatrixElem(DiffCandP0, 0, 2) = vf[i].center->val[2] - Tpo->val[2];
        //iftPrintMatrix(DiffCandP0);


        v3 = columnVectorMatrixToVector(DiffCandP0);

        DiffShiftDotNj = iftVectorInnerProduct(v1,v3);
        //DiffShiftDotNj = MatrixInnerProduct(Nj,DiffCandP0);

        //if(NdotNj == 0.000000){
        //  continue;
        //}

        lambda[i]=(float) DiffShiftDotNj / NdotNj;


        v.x = Tpo->val[0] + lambda[i] * Tn->val[0];
        v.y = Tpo->val[1] + lambda[i] * Tn->val[1];
        v.z = Tpo->val[2] + lambda[i] * Tn->val[2];

        if (isValidPoint(img, v))
        {
          if (lambda[i] < min){
              p1->x = v.x;
              p1->y = v.y;
              p1->z = v.z;
              min = lambda[i];
          }
          if (lambda[i] > max) {
              pn->x = v.x;
              pn->y = v.y;
              pn->z = v.z;
              max = lambda[i];
          }
        }
      }
    }

    iftDestroyMatrix(&Nj);
    iftDestroyMatrix(&DiffCandP0);


    if ((p1->x != -1) && (pn->x != -1)){
      return 1;
    }
    else
      return 0;
}



iftVolumeFaces* createVF(iftImage* I)
{
  int Nx = I->xsize;
  int Ny = I->ysize;
  int Nz = I->zsize;
  int i;

  iftVolumeFaces *vf = (iftVolumeFaces *) malloc(sizeof(iftVolumeFaces) * 6);

  for (i = 0; i < 6; i++)
  {
      vf[i].orthogonal = iftCreateMatrix(1, 4);
      vf[i].center = iftCreateMatrix(1, 4);
  }

  // Face of Plane XY
  vf[0].orthogonal->val[0] = 0;
  vf[0].orthogonal->val[1] = 0;
  vf[0].orthogonal->val[2] = -1;
  vf[0].orthogonal->val[3] = 1;

  vf[0].center->val[0] = Nx / 2;
  vf[0].center->val[1] = Ny / 2;
  vf[0].center->val[2] = 0;
  vf[0].center->val[3] = 1;

  // Face of Plane XZ
  vf[1].orthogonal->val[0] = 0;
  vf[1].orthogonal->val[1] = -1;
  vf[1].orthogonal->val[2] = 0;
  vf[1].orthogonal->val[3] = 1;

  vf[1].center->val[0] = Nx / 2;
  vf[1].center->val[1] = 0;
  vf[1].center->val[2] = Nz / 2;
  vf[1].center->val[3] = 1;

  // Face of Plane YZ
  vf[2].orthogonal->val[0] = -1;
  vf[2].orthogonal->val[1] = 0;
  vf[2].orthogonal->val[2] = 0;
  vf[2].orthogonal->val[3] = 1;

  vf[2].center->val[0] = 0;
  vf[2].center->val[1] = Ny / 2;
  vf[2].center->val[2] = Nz / 2;
  vf[2].center->val[3] = 1;

  // Face of Opposite Plane XY
  vf[3].orthogonal->val[0] = 0;
  vf[3].orthogonal->val[1] = 0;
  vf[3].orthogonal->val[2] = 1;
  vf[3].orthogonal->val[3] = 1;

  vf[3].center->val[0] = Nx / 2;
  vf[3].center->val[1] = Ny / 2;
  vf[3].center->val[2] = Nz - 1;
  vf[3].center->val[3] = 1;

  // Face of Opposite Plane XZ
  vf[4].orthogonal->val[0] = 0;
  vf[4].orthogonal->val[1] = 1;
  vf[4].orthogonal->val[2] = 0;
  vf[4].orthogonal->val[3] = 1;

  vf[4].center->val[0] = Nx / 2;
  vf[4].center->val[1] = Ny - 1;
  vf[4].center->val[2] = Nz / 2;
  vf[4].center->val[3] = 1;

  // Face of Opposite Plane YZ
  vf[5].orthogonal->val[0] = 1;
  vf[5].orthogonal->val[1] = 0;
  vf[5].orthogonal->val[2] = 0;
  vf[5].orthogonal->val[3] = 1;

  vf[5].center->val[0] = Nx-1;
  vf[5].center->val[1] = Ny / 2;
  vf[5].center->val[2] = Nz / 2;
  vf[5].center->val[3] = 1;

  return vf;
}

void DestroyVF(iftVolumeFaces *vf)
{
    int i;

    for (i = 0; i < 6; i++)
    {
        iftDestroyMatrix(&vf[i].orthogonal);
        iftDestroyMatrix(&vf[i].center);
    }
    free(vf);
}


iftMatrix *voxelToMatrix(iftVoxel v)
{
    iftMatrix *voxMat = iftCreateMatrix(1, 4);
    iftMatrixElem(voxMat, 0, 0) = v.x;
    iftMatrixElem(voxMat, 0, 1) = v.y;
    iftMatrixElem(voxMat, 0, 2) = v.z;
    iftMatrixElem(voxMat, 0, 3) = 1;

    return voxMat;
}

iftMatrix *imagePixelToMatrix(iftImage *img, int p)
{
    iftMatrix *pixMat = iftCreateMatrix(1, 4);
    iftMatrixElem(pixMat, 0, 0) = p % img->xsize;
    iftMatrixElem(pixMat, 0, 1) = p / img->xsize;
    iftMatrixElem(pixMat, 0, 2) = 0;
    iftMatrixElem(pixMat, 0, 3) = 1;

    return pixMat;
}



iftImage *MaximumIntensityProjection(iftImage *img, float xtheta, float ytheta)
{
    float diagonal = 0;
    double maxIntensity;
    int p=0;
    int Nu, Nv;
    iftVoxel p1, pn;

    iftVolumeFaces* volumeFaces;
    iftMatrix *Mtemp, *T;
    iftMatrix *Norigin, *Tnorigin;
    iftMatrix *Tpo;

    diagonal = sqrt((img->xsize * img->xsize) + (img->ysize * img->ysize) + (img->zsize * img->zsize));
    Nu = Nv = diagonal;
    iftImage *output = iftCreateImage(Nu, Nv, 1);
    T  = createTransformationMatrix(img, xtheta, ytheta);
    //iftPrintMatrix(T);

    Norigin =  iftCreateMatrix(1, 4);
    iftMatrixElem(Norigin, 0, 0) = 0;
    iftMatrixElem(Norigin, 0, 1) = 0;
    iftMatrixElem(Norigin, 0, 2) = -1;
    iftMatrixElem(Norigin, 0, 3) = 0;
    //iftPrintMatrix(Norigin);

    Tnorigin = iftMultMatrices(T, Norigin);

    volumeFaces = createVF(img);


    for (p = 0; p < output->n; p++)
    {
        printf("Step : %d\n", p);
        Mtemp = imagePixelToMatrix(output,p);
        iftMatrixElem(Mtemp, 0, 2) = diagonal/2;

        Tpo = iftMultMatrices(T, Mtemp);

        if (ComputeIntersection(Tpo, img, Tnorigin, volumeFaces, &p1, &pn))
        {
            maxIntensity  = DDA(img,p1,pn);

            output->val[p] = maxIntensity;
        }


        iftDestroyMatrix(&Mtemp);
        iftDestroyMatrix(&Tpo);
    }

    iftDestroyMatrix(&T);
    iftDestroyMatrix(&Norigin);
    iftDestroyMatrix(&Tnorigin);
    DestroyVF(volumeFaces);

    return output;
}


int main(int argc, char *argv[])
{
    //if (argc != 6)
    //    Error("Run: ./main <filename> <output> <xtheta> <ytheta> , "main");

    char buffer[512];

    float tx, ty;
    tx = atof(argv[3]);
    ty = atof(argv[4]);
    //tz = atof(argv[5]);
    char *imgFileName = iftCopyString(argv[1]);
    iftImage *img = iftReadImageByExt(imgFileName);

    iftImage *output = NULL;

    output = MaximumIntensityProjection(img, tx, ty);
    sprintf(buffer, "data/%.1f%.1f%s", tx, ty, argv[2]);
    iftImage *normalizedImage= iftNormalize(output,0,255);

    iftWriteImageByExt(normalizedImage, buffer);
    iftDestroyImage(&img);
    iftDestroyImage(&output);
    return 0;
}
