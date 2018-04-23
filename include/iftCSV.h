/**
 * @file iftCSV.h
 * @brief Generic CSV manipulations.
 * @author Samuel Martins
 * @date Sep 15, 2015
 *
 * @note Programs:
 * * @ref iftTestCSV.c = It shows how to open a CSV file.
 */

#ifndef IFT_CSV_H
#define IFT_CSV_H

#ifdef __cplusplus
extern "C" {
#endif


#include "iftCommon.h"
#include "iftSList.h"


/**
 * @brief Comma-Separated Values struct.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 */
typedef struct ift_csv {
    /** Number of rows of the CSV matrix. */
    size_t nrows;
    /** Number of columns of the CSV matrix. */
    size_t ncols;
    /** CSV matrix of strings. Each string has 512 characters. */
    char ***data;
} iftCSV;


/**
 * @brief Creates an iftCSV (and allocates their string matrix).
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * @param nrows Number of rows of the CSV matrix.
 * @param ncols Number of columns of the CSV matrix.
 * @return An iftCSV.
 */
iftCSV *iftCreateCSV(size_t nrows, size_t ncols);


/**
 * @brief Destroys an iftCSV.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * @param csv Reference to the iftCSV to be destroyed.
 */
void iftDestroyCSV(iftCSV **csv);

/**
 * @brief Prints an iftCSV.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * @param csv iftCSV to be printed.
 */
void iftPrintCSV(iftCSV *csv);


/**
 * @brief Reads an iftCSV from file.
 * @author Samuel Martins
 * @date Sep 15, 2015
 * @ingroup CSV
 *
 * Reads an iftCSV from file.
 * The CSV must have exactly the following separating pattern:
 *
 * Brazil,100,20.50
 * France,256,35.00
 *
 * @param csv_pathname The pathname from the CSV to be read.
 * @return The iftCSV.
 */
iftCSV *iftReadCSV(const char *csv_pathname, const char separator);


/**
 * @brief Writes a CSV structure into a file.
 * @author Thiago Vallin Spina
 * @date Feb 23, 2016
 * @ingroup CSV
 *
 * @param csv The original CSV structure.
 * @param filename The output filename.
 * @param separator The selected separating character.
 *
 * @warning This function issues an error if separator != ',' and separator != ';'
 */
void iftWriteCSV(iftCSV *csv, const char *filename, const char separator);


#ifdef __cplusplus
}
#endif



#endif //IFT_CSV_H
