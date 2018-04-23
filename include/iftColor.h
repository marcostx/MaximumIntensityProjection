/**
 * @file iftColor.h
 * @brief Structs and function prototypes for color.
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 * @note Programs:
 * @ref iftGenerateColorTable.c (demo/Miscellaneous/iftGenerateColorTable.c) = Generates a Color Table. 
 */

#ifndef IFT_COLOR_H_
#define IFT_COLOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"

/**
 * @brief Color space enumeration
 * @note added by Adan Echemendia. The three color components from LABNorm2_CSPACE are normalized from 0 to 1. Not so
 * the components from LABNorm_CSPACE.
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
typedef enum ift_color_space {
    YCbCr_CSPACE,
    YCbCrNorm_CSPACE,
    RGB_CSPACE,
    RGBNorm_CSPACE,
    GRAY_CSPACE,
    GRAYNorm_CSPACE,
    WEIGHTED_YCbCr_CSPACE,
    LAB_CSPACE,
    LABNorm_CSPACE,
		LABNorm2_CSPACE
} iftColorSpace;

#define WHITEPOINT_X	0.950456
#define WHITEPOINT_Y	1.0
#define WHITEPOINT_Z	1.088754

#define LABF(t)	\
	((t >= 8.85645167903563082e-3) ? \
	pow(t,0.333333333333333) : (841.0/108.0)*(t) + (4.0/29.0))

#define LABINVF(t)	\
	((t >= 0.206896551724137931) ? \
	((t)*(t)*(t)) : (108.0/841.0)*((t) - (4.0/29.0)))


/**
 * @brief Float color space structure
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
typedef struct ift_fcolor {
	float val[3];
} iftFColor;


/**
 * @brief Integer color space structure
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 */
typedef struct ift_color {
  int val[3];
	float alpha;
} iftColor;


/**
 * @brief Array of colors
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 */
typedef struct ift_colortable {
  iftColor *color;
  int ncolors;
} iftColorTable;


/**
 * @brief Creates a color array by generating colors in a smart way.
 * @author Samuel Martins
 * @date Jun 15, 2016
 * @ingroup Color
 */
iftColorTable *iftCreateColorTable(int ncolors, iftColorSpace cspace);
iftColorTable *iftCreateHeatMapColorTable(float delta, int scalingFactor);
void iftConvertRGBColorTableToYCbCrColorTable(iftColorTable *colorTable_rgb,int normalization_value);

/**
 * @brief Destroy color array from memory
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
void iftDestroyColorTable(iftColorTable **ctb);

/**
 * @brief Create a color array from blue to red
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColorTable *iftBlueToRedColorTable(int ncolors);

/**
 * @brief Converts an RGB color t o the iftColor data structure (does not convert it to YCbCr).
 *
 * @author Thiago Vallin Spina
 * @date Feb 14, 2016
 *
 * @param R The Red component.
 * @param G The Green component.
 * @param B The Blue component.
 *
 * @return The iftColor structure with the RGB color.
 */
iftColor  iftRGBColor(int R, int G, int B);

/**
 * @brief Convert color from RGB to YCbCr
 * http://www.equasys.de/colorconversion.html#YCbCr-RGBColorFormatConversion (for 24-bit color images)
 * This code is actually converting to YCgCo 
 *
 * @author Alexandre Falcao
 * @date Jul 30, 2016
 * @ingroup Color
 *
 */
iftColor  iftRGBtoYCbCr(iftColor cin, int normalization_value);

/**
 * @brief Convert color from YCbCr to RGB \n
 * http://www.equasys.de/colorconversion.html#YCbCr-RGBColorFormatConversion \n
 * http://www.bluebit.gr/matrix-calculator/calculate.aspx 
 * This code is actually converting from YCgCo 
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftYCbCrtoRGB(iftColor cin, int normalization_value);

/**
 * @brief Convert color from RGB to HSV
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftRGBtoHSV(iftColor cin, int normalization_value);

/**
 * @brief Convert color from HSV to RGB
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftHSVtoRGB(iftColor cin, int normalization_value);

/**
 * @brief Convert color from RGB to Lab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftFColor iftRGBtoLab(iftColor rgb, int normalization_value);

/**
 * @brief Convert color from RGB to Normalized Lab. The function calls iftRGBtoLab function and normalizes each value to get in the interval [0..1]
 * @author Adan Echemendia
 * @date Jul 26, 2017
 * @ingroup Color
 *
 */
iftFColor iftRGBtoLabNorm2(iftColor rgb, int normalization_value);

/**
 * @brief Convert color from YcbCr to Lab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */  
iftFColor iftYCbCrtoLab(iftColor ycbcr, int normalization_value);

/**
 * @brief Convert color from YCbCr to LabNorm
 *
 * @author Cesar Castelo
 * @date Jan 10, 2018
 * @ingroup Color
 *
 */
iftFColor iftYCbCrtoLabNorm(iftColor ycbcr, int normalization_value);

/**
 * @brief Convert color from YcbCr to LabNorm (all components are in the interval [0..1])
 * @author Adan Echemendia
 * @date Jul, 2017
 * @ingroup Color
 *
 */
iftFColor iftYCbCrtoLabNorm2(iftColor ycbcr, int normalization_value);
/**
 * @brief Convert color from RGB to LabNorm
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftFColor iftRGBtoLabNorm(iftColor rgb, int normalization_value);

/**
 * @brief Convert color from Lab to RGB
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftLabtoRGB(iftFColor lab, int normalization_value);

/**
 * @brief Convert color from RGB to YCbCr BT.2020
 *
 * http://http://en.wikipedia.org/wiki/Rec._2020
 *
 * Book Digital Video Concepts, Methods, and Metrics: Quality, Compression, Performance,
 * and Power Trade-off Analysis, ISBN 9781430267133, page 27
 * https://books.google.com.br/books?id=fr08BQAAQBAJ
 *
 * @param cin Input RGB color
 * @param rgbBitDepth Input bit depth
 * @param yCbCrBitDepth Output bit depth (10 or 12)
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftRGBtoYCbCrBT2020(iftColor cin, const int rgbBitDepth, const int yCbCrBitDepth);

/**
 * @brief Convert color from YCbCr BT.2020 to RGB
 *
 * http://http://en.wikipedia.org/wiki/Rec._2020
 *
 * Book Digital Video Concepts, Methods, and Metrics: Quality, Compression, Performance,
 * and Power Trade-off Analysis, ISBN 9781430267133, page 27
 * https://books.google.com.br/books?id=fr08BQAAQBAJ
 *
 * @param cin Input RGB color
 * @param rgbBitDepth Input bit depth
 * @param yCbCrBitDepth Output bit depth (10 or 12)
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftYCbCrBT2020toRGB(iftColor cin, const int yCbCrBitDepth, const int rgbBitDepth);

/**
 * @brief Convert color from Lab to QLab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftColor  iftLabtoQLab(iftFColor lab,int normalization_value);

/**
 * @brief Convert color from QLab to Lab
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Color
 *
 */
iftFColor iftQLabToLab(iftColor qlab, int normalization_value);

/**
 * @brief Converts values from [0,1] to R, G, B colors according to the heat color map (from blue to red). 
 * @param value: input value 
 * @param R: output red value
 * @param G: output green value 
 * @param B: output blue value 
 * @author Alexandre Falcao
 * @date Aug, 5th, 2016
 * @ingroup Color
 */
void iftHeatColorMapping(float value, float *R, float *G, float *B);


/**
 * @brief Create the Hot Iron Color Table (based on http://dicom.nema.org/dicom/2013/output/chtml/part06/chapter_B.html)
 * with RGB colors.
 * The maximum number of colors is 256.
 * 
 * @param n_colors Number of colors. Maximum of 256
 * @return Hot Iron Color Table.
 * 
 * @author Samuka Martins
 * @date Mar 12, 2018.
 */
iftColorTable *iftCreateHotIronColorTable(int n_colors);



#ifdef __cplusplus
}
#endif

#endif
