#ifndef IFT_GRAPHICS_H_
#define IFT_GRAPHICS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftColor.h"
#include "iftSet.h"
#include "iftImage.h"
#include "iftPlane.h"
#include "iftImageMath.h"
#include "iftFImage.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftSeeds.h"
#include "iftDataSet.h"
#include "iftRepresentation.h"
#include "iftMathMorph.h"
#include "iftInterpolation.h"
#include "iftCurve.h"

#define  NUM_OF_NORMALS    65161
#define  SCENE_NORMAL      0
#define  OBJECT_NORMAL     1
#define  RAYCASTING        0
#define  SPLATTING         1

/* Graphical context with the visualization attributes for ray casting
   and voxel splatting, including current viewing direction
   parameters.

       Octants of the Scene
 z
  /
 /
 ---- x
 |
y|
       4 ------- 5
      /        /|
     /        / |
    /        /  |
  0 -------- 1  |
   |  6     |  / 7
   |        | /
   |        |/
  2 -------- 3


The origin of the image to be rendered on the viewing plane is
initially at (0,0,-diag/2). The observer is initially at (diag/2,
diag/2, -infinity) with a single white light source at the same
position. The transformation T is applied to the scene for voxel
splatting or the inverse of T is applied to the vieweing plane for ray
casting. The front-to-back (FTB) voxel access order for splatting
depends on the pair (closest octant, principal axis). The FTB
possibilities (dx,dy,dz,xo,yo,zo,xf,yf,zf) are stored in a 2D array
whose column is the closest octant and line is the principal axis. Any
change in viewing direction requires to update all viewing parameters below.

*/

/* Extra surface rendering buffers, especially required for voxel
   splatting. They all have the same dimensions of the resulting
   rendered image (diag x diag pixels) to avoid clipping. */

typedef struct ift_surface_rendering_buffers {
  float  depth;         /* distance to the viewing plane, whose center
			   is initially positioned at
			   (diag/2,diag/2,-diag/2) */
  float  opacity;       /* accumulated opacity of the voxels projected
			   onto the viewing plane */
  int    voxel;         /* voxels in the scene that have been projected
			   onto the viewing plane */
  int    object;        /* object in the scene whose voxel has been
			   projected onto the viewing plane */
} iftSRBuffers;

/* Object opacity, visibility, and color attributes in the scene */

typedef struct ift_object_attributs {
  float opacity;          /* opacity value in [0,1] */
  float red, green, blue; /* proportions in [0,1] of red, green and
			     blue components in its color */
  char visibility;        /* visibility flag (0/1 for invisible/visible) */
} iftObjectAttributes;

/* Parameters of the Phong's illumination model */

typedef struct ift_phong_model {
  float      ka;     /* ambient reflection constant  (0.1 by default) */
  float      kd;     /* diffuse reflection constant  (0.7 by default) */
  float      ks;     /* specular refrection constant (0.2 by default) */
  float      ns;     /* exponent of the specular reflection component (5 by default) */
  iftVector *normal; /* normal table used to speed up shading */
  float     *Idist;  /* distance table used to convert voxel's distance to the viewing plane into intensity in [0,1] (also known as depth shading) */
  int        ndists; /* number of distance values, which is the size of the scene's diagonal. */
} iftPhongModel;

/* Viewing direction and FTB voxel access order */

typedef struct ift_ftb_access {
  int dx, dy, dz; /* -1 or 1. Default is 1,1,1 */
  int xo, yo, zo; /* the x,y,z original coordinates of the closest  octant: the initial FTB voxel access coordinates */
  int xf, yf, zf; /* the x,y,z original coordinates of the farthest octant: the final FTB voxel access coordinates */
} iftFTBaccess;

typedef struct ift_viewing_direction {
  iftMatrix           *T;            /* scene's transformation matrix: default is spin=tilt=0.0 */
  iftMatrix           *Tinv;         /* inverse of the scene's transformation matrix */
  iftMatrix           *R;            /* rotation matrix, only for vectors */
  iftMatrix           *Rinv;         /* inverse of the rotation matrix, only for vectors */
  char                 paxis;        /* principal axis, 0 for 'x', 1 for 'y', or 2 for 'z', which is the last one to be visited in the voxel splatting loop. This is the most orthogonal axis to the viewing plane. */
  int                  octant;       /* closest octant, 0, 1, ..., 7, to the vieweing plane, which together with the principal axis indicate the front-to-back voxel access order for splatting. */
  iftFTBaccess        *ftb;          /* gives the front-to-back access order from the closest octant to the farthest octant. Used for voxel splatting only. */
  iftVector            V;            /* visualization direction --- vector that points to the observer. */
} iftViewDir;


/* Graphical context of the scene for 3D visualization */

typedef struct ift_graphic_context {
  iftObjectAttributes *object;       /* list of attributes per object
					in the scene */
  int                  nobjs;        /* number of objects in the scene */
  float                overall_opac; /* it indicates the overall opacity of the objects in the scene */
  char                 proj_mode;    /* voxel projection mode: ray casting is 0 and voxel splatting is 1 */
  iftPhongModel       *phong;        /* rendering parameters of the Phong's illumination model */
  iftViewDir          *viewdir;      /* viewing direction parameters */
  iftSRBuffers        *surf_render;  /* extra surface rendering buffers, especially required for voxel splatting */
  iftFImage           *scene;        /* the 3D image property for visualization (intensity, BMD, etc.) */
  iftImage            *label;        /* the 3D label image with objects for visualization, used only for surface rendering */
  iftImage            *border;        /* the 3D labeled borders of the objects in label, used only for surface rendering */
  iftPlane             face[6];      /* planes of the 6 faces of the scene  */
  iftImage            *normal;       /* normal indices of the scene's voxels */
  iftFImage           *opacity;      /* opacity scene used for volume rendering only */
} iftGraphicalContext;

/*
   Graphical Context for 3D visualization of the image properties in
   scene. The object labels i > 0 (0 is background) are used only for
   surface rendering (i.e., label==NULL for volume rendering). The
   normal vectors may be estimated as the gradient of the scene
   (normal_type=SCENE_NORMAL), as the gradient of the EDT of the
   objects (normal_type=OBJECT_NORMAL) or on-the-fly from the index
   buffer (normal_type=NIL).

*/

void                 iftSetSceneNormal(iftGraphicalContext *gc);
void                 iftSetProjectionMode(iftGraphicalContext *gc, char proj_mode);
void                 iftSetObjectNormal(iftGraphicalContext *gc, char normal_type);
void                 iftSetObjectColor(iftGraphicalContext *gc, int object, float red, float green, float blue);
void                 iftSetSceneOpacity(iftGraphicalContext *gc, float min_val, float max_val, iftImage *grad, int grad_thres, float max_opac);
void                 iftSetObjectOpacity(iftGraphicalContext *gc, int object, float opacity);
void                 iftSetObjectVisibility(iftGraphicalContext *gc, int object, char visibility);
void                 iftSetViewDir(iftGraphicalContext *gc, float tilt, float spin);
char                 iftIntersectionPoints(iftGraphicalContext *gc, iftPoint P0, iftVector n, iftPoint *P1, iftPoint *Pn);
iftImage            *iftSurfaceRender(iftGraphicalContext *gc);
iftImage            *iftVolumeRender(iftGraphicalContext *gc);
iftGraphicalContext *iftCreateGraphicalContext(iftFImage *scene, iftImage *label);
void                 iftDestroyGraphicalContext(iftGraphicalContext *gc);

iftImage  *iftGraphicImage(iftImage *img); /* creates image for overlaying drawings */
iftImage  *iftDrawPath(iftImage *img, iftImage *pred, int last, iftFColor normRGB, iftAdjRel *B); /* Draw a path in the predecessor map with terminus at point last */
void iftDrawPoint(iftImage *img, iftVoxel u, iftColor YCbCr, iftAdjRel *B, int rangeValue);
void       iftDrawPoints(iftImage *img, iftSet *S, iftColor YCbCr, iftAdjRel *B);
void       iftDrawCurve(iftImage *img, iftCurve *curve, iftColor YCbCr, iftAdjRel *B);
void       iftDrawLine(iftImage *img, iftVoxel u1, iftVoxel u2, iftColor YCbCr, iftAdjRel *B);
void       iftDrawBorders(iftImage *img, iftImage *label, iftAdjRel *A, iftColor YCbCr, iftAdjRel *B);
void 	   iftDrawLabels(iftImage *img, iftImage *label, iftColorTable *cmap,	iftAdjRel *A, uchar fill, float alpha);
void       iftDrawObject(iftImage *img, iftImage *bin, iftColor YCbCr, iftAdjRel *B);
void       iftDrawRoots(iftImage *img, iftImage *root, iftColor YCbCr, iftAdjRel *B);
iftImage *iftDraw2DFeatureSpace(iftDataSet *Z, uchar opt, uchar status);
iftImage *iftDrawVoxelSamples(iftDataSet *Z, uchar opt, char *ref_data_type);
iftImage *iftColorizeComp(iftImage *label);
iftImage *iftColorizeCompOverImage(iftImage *orig, iftImage *label);
iftImage  *iftProjectMaxValue(iftImage *img,  iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMinValue(iftImage *img,  iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMeanValue(iftImage *img, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMaxObjectValue(iftImage *img,  iftImage *obj, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMinObjectValue(iftImage *img,  iftImage *obj, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMeanObjectValue(iftImage *img, iftImage *obj, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftCurvilinearSplatting(iftImage *img, iftFImage *dist, iftPlane *cutplane, int viewsize, float depth);
void       iftTiltSpinFromPrincipalAxis(iftImage *bin, float *tilt, float *spin);

/**
 * @brief Draws a line from px1 to pxn recursively using color YCbCr
 * for a given adjacency relation. The main advantage of this function
 * over iftDrawLine is that the line will be "perfectly" rasterized,
 * guaranteeing that endpoints px1 and pxn are part of it and that
 * there will be no "holes" on the line.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img colored Image.
 * @param px1 start voxel.
 * @param pxn end voxel.
 * @param YCbCr YCbCr color.
 * @param B adjacency used to draw the line.
*/
void iftDrawLineRecursive(iftImage *img, iftVoxel px1, iftVoxel pxn, iftColor YCbCr, iftAdjRel *B);

/**
 * @brief Draws the lines of a 2D polygon given by a set of points and then fills it to form a binary mask.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param points The set of vertices.
 * @param n      The number of points.
 * @param xsize  The X image dimensions.
 * @param ysize  The Y image dimensions.
 * @return The X tile coordinate.
 */
iftImage *iftMaskImageFromPolygon2D(iftVoxel *points, int n, int xsize, int ysize);


/**
 * @brief Draws the borders of a given Bounding Box in a pre-allocated Image.
 * The border will have value/color <b>YCbCr</b>.
 * @author Samuel Martins
 * @date May 30, 16
 * 
 * @param img   Pre-allocated Image where the borders will be drawn.
 * @param bb    Bounding box.
 * @param YCbCr Value/Color used to draw the borders.
 */
void iftDrawBoundingBoxBordersInPlace(iftImage *img, iftBoundingBox bb, iftColor YCbCr);


#ifdef __cplusplus
}
#endif

#endif

