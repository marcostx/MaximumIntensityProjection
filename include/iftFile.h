/**
 * @file iftFile.h
 * @brief Structs and function prototypes for file management.
 * 
 * @author Samuel Martins
 * @date Aug 10, 2015
 * @ingroup File
 * 
 * @note Programs:
 * * @ref iftTestLoadDir.c = It loads a directory.
 * * @ref iftMakeDirsExamples.c = It shows examples of how to create Pathnames, Files, and Directories.
 */

#ifndef IFT_FILE_H
#define IFT_FILE_H

#ifdef __cplusplus
extern "C" {
#endif
    
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined(__WIN32) || defined(__WIN64)
#include <Windows.h>
#endif

#include "iftCommon.h"
#include "iftCSV.h"
#include "iftDict.h"
#include "iftSList.h"
#include "iftString.h"
#include "iftUtil.h"


// DIRECTORY SEPARATOR STRING
/**
 * @brief Directory Separator String.
 * @author Samuel Martins
 * @ingroup File
 */
#if defined(__WIN32) || defined(__WIN64)
    #define IFT_SEP_C "\\"
#else
    #define IFT_SEP_C "/"
#endif

/**
 * @brief An abstract representation of file.
 *
 * \ref iftFile contains the file info for creating, using and managing disk files.
 *
 * @author Samuel Martins
 * @ingroup File
 */
typedef struct ift_file {
    /** File path */
    char *path;
    /**Label for machine learning dataset files*/
    int label;
    /** Sample id following the new convention for IFT datasets (Thiago Vallin Spina, Mar 09, 2016) */
    int sample;
    /** Suffix following the new convention for IFT datasets (Thiago Vallin Spina, Mar 09, 2016) */
    char *suffix;

    /**Status for machine learning dataset files
     * iftFile can be marked as TRAIN, TEST, OUTLIER as the iftSample*/
    char status;
} iftFile;


/**
 * @brief A struct that contains an array to iftFile*.
 * @author Samuel Martins
 * @date Oct 6th, 2015
 * @ingroup File
 */
typedef struct ift_fileset {
    /** Number of Files */
    size_t n;
    /** Array of Files */
    iftFile **files;
} iftFileSet;


/**
 * @brief An abstract representation of directories. Contains a list of files and subdirectories under that directory.
 * @author Samuel Martins
 * @ingroup File
 */
typedef struct ift_dir {
    char *path;
    /** Number of files from the directory */
    size_t nfiles;
    /** Vector of iftFiles with the children */
    iftFile **files;
    /** Number of files from the directory */
    size_t nsubdirs;
    /** Vector of iftDir with the subdirs */
    struct ift_dir **subdirs;
} iftDir;

/********************** Pathname Management ***********************/
/**
 * @brief Creates an iftFile.
 * @author Samuel Martins
 * @date Nov 10, 2015
 * @ingroup File
 *
 * Creates an iftFile from a <code><b>pathname</b></code> on memory, BUT IT DOES NOT CREATE IT ON DISK.
 * If the <code><b>pathname</b><code> exists, the field <exists> is set to 1.
 * Otherwise, it is set to 0.
 * If the file exists and is a directory, an error will be raised.
 *
 * 
 * @param pathname The file to be read.
 * @return An iftFile with the <code><b>pathname</b><code>.
 * @sa iftLoadDir()
 */
iftFile *iftCreateFile(const char *pathname);


/**
 * @brief Copies an iftFile.
 * @author Peixinho
 * @ingroup File
 *
 * @param file The iftFile to be copied.
 * @return A copy of the input iftFile.
 */
iftFile *iftCopyFile(const iftFile* file);


/**
 * @brief Copies any <b>File</b> (text, binary, pdf, doc, png, ...) or <b>Directory</b> from DISK.
 * @author Samuel Martins
 * @date Dec 1, 2015
 * @ingroup File
 *
 * @param src Filename from the source file/directory to be copied.
 * @param dst Destination Filename from the copy. It can be a directory filename where the src will be copied.
 *
 * @attention If the <code><b>dst</b></code> filename is an existing directory, the <b>src</b> (file or directory) will be copied into the <b>dst</b> directory with the same src basename.
 *
 * @exception src filename does not exist.
 * @exception src filename is a <i>Directory</i> and the <code><b>dst</b></code> is a <i>File</i>.
 */
void iftCopyFromDisk(const char *src, const char *dst);


/**
 * Count number of unique label on a file set.
 * @author Peixinho
 * @date Jun, 2017
 * @param fileset
 * @return Number of unique labels
 */
int iftFileSetLabelsNumber(const iftFileSet* fileset);


/**
 * @brief Loads a Directory <dir_pathname>.
 * @author Samuel Martins
 * @date Aug 14, 2015
 * @ingroup File
 *
 * Loads a Directory from a <dir_pathname> on memory.
 * THE DIRECTORY MUST EXIST.
 * The parameter @p hier_levels indicates until which hierarchical level will be considered in the
 * loading.
 *
 * If <code>hier_levels=2</code>, all files and subdirs until the 2º level from the root directory are loaded, but
 * the content of the subdirs from the 2º level are not loaded.
 * If  <code>hier_levels=0</code>, all files and subdirs from all levels are loaded.
 * If <code>hier_levels is</code> greater than the real number of levels from the root directory, all files and subdirs
 * will be loaded.
 *
 * A separator char '/' (linux) or '\' (windows) is put at the end of the pathname.
 * If the file exists and is not a directory (it is a file), an error will be raised.
 *
 * @param dir_pathname The directory to be read.
 * @param hier_levels Number of levels that will be loaded. The vale 0 loads all files from the whole hierarchy.
 * @return An iftDir with the directory.
 */
iftDir *iftLoadDir(const char *dir_pathname, size_t hier_levels);


/**
 * @brief Loads all files (not subdirs) in the directory <dir_pathname> from a given <extension>.
 * @author Samuel Martins
 * @ingroup File
 *
 * Loads all files (not subdirs) in the directory <b><code>dir_pathname</code></b> from a given <b><code>extension</code></b>.
 * The filenames are sorted in ascending order.
 * It only gets the files from the 1st file hierarchical level.
 * If <b><code>extension</code></b> = "", it gets all files (not subdirs).
 *
 * @param dir_pathname The directory to be read.
 * @param suffix Suffix from the files.
 * @return An iftDir with all files with <extension> inside <dir_pathname>.
 */
iftDir *iftLoadFilesFromDirBySuffix(const char *dir_pathname, const char *suffix);


/**
 * @brief Loads all files (not subdirs) in the directory <dir_pathname> from a given <regex>.
 * @author Samuel Martins
 * @author Jul 7, 2016
 * @ingroup File
 *
 * Loads all files (not subdirs) in the directory <b><code>dir_pathname</code></b> from a given <b><code>regex</code></b>. \n
 * The filenames are sorted in ascending order. \n
 * It only gets the files from the 1st file hierarchical level. \n
 * Example of Regex: \n
 * regex = ".*Pos.*\\.scn", equivalent to ls *Pos*.scn in bash-linux
 *
 * @param dir_pathname The directory to be read.
 * @param regex        Target Regex.
 * @return             An iftDir with all files with <extension> inside <dir_pathname> filtered by <regex>.
 */
iftDir *iftLoadFilesFromDirByRegex(const char *dir_pathname, const char *regex);


/**
 * @brief Loads an Array of iftFiles from a Directory (in the entire directory hierarchy).
 *
 * The parameter @p hier_levels indicates until which hierarchical level will be considered in the
 * loading.
 *
 * If <code>hier_levels=2</code>, all files and subdirs until the 2º level from the root directory are loaded, but
 * the content of the subdirs from the 2º level are not loaded.
 * If  <code>hier_levels=0</code>, all files and subdirs from all levels are loaded.
 * If <code>hier_levels is</code> greater than the real number of levels from the root directory, all files and subdirs
 * will be loaded.
 *
 * @author Samuel Martins
 * @date Oct 6th, 2015
 * @ingroup File
 *
 * @param dir_pathname The pathname from the directory to be read.
 * @param hier_levels Number of folder levels that will be loaded.
 * @return The array of iftFiles.
 *
 * @exception dir_pathname does not exist or it is a File.
 */
iftFileSet *iftLoadFileSetFromDir(const char *dir_pathname, size_t hier_level);


/**
 * @brief Get the file labels, according to the IFT file format.
 * @author Peixinho
 * @ingroup File
 *
 * The image label is stored as part of the image file name.
 * The IFT default nomenclature is <label>_<id>.<ext>
 *
 * Ex.: 000002_000015.pgm
 *
 * @param dir FileSet containing the files.
 * @return The file labels.
 */
int* iftFileSetLabels(const iftFileSet* fileSet);


/**
 * @brief Loads an Array of iftFiles with extension <b><code>ext</code></b> from a Directory.
 * @author Samuel Martins
 * @date Oct 7th, 2015
 * @ingroup File
 *
 * Loads an Array of iftFiles with extension <ext> from a Directory.
 * The function gets the files in the entire directory hierarchy.
 *
 * @param dir_pathname The pathname from the directory to be read.
 * @return The array of iftFiles.
 */
iftFileSet *iftLoadFileSetFromDirBySuffix(const char *dir_pathname, const char *ext);


/**
 * @brief Loads an Array of iftFiles with Regex <b><code>regex</code></b> from a Directory.
 * @author Samuel Martins
 * @date Jul 7, 2015
 * @ingroup File
 *
 * @note The function gets the files in the entire directory hierarchy.
 *
 * @param dir_pathname The pathname from the directory to be read.
 * @param regex        Target Regex.
 * @return The array of iftFiles filtered by the regex.
 */
iftFileSet *iftLoadFileSetFromDirByRegex(const char *dir_pathname, const char *regex);


/**
 * @brief Reads an Array of iftFiles from a CSV.
 * @author Samuel Martins
 * @date Sep 25th, 2015
 * @ingroup File
 *
 * Reads an iftCSV from file.
 * The CSV must have exactly the following separating pattern:
 *
 * Brazil,100,20.50
 * France,256,35.00
 *
 * @param csv_pathname The pathname from the CSV to be read.
 * @param sort_pathnames It indicates if the pathnames should be sorted.
 * @return The array of iftFiles.
 */
iftFileSet *iftLoadFileSetFromCSV(const char *csv_pathname, bool sort_pathnames);


/**
 * @brief Reads an Array of iftFiles from a Directory or a CSV file.
 *
 * The parameter @p hier_levels indicates until which hierarchical level will be considered in the
 * loading.
 *
 * If <code>hier_levels=2</code>, all files and subdirs until the 2º level from the root directory are loaded, but
 * the content of the subdirs from the 2º level are not loaded.
 * If  <code>hier_levels=0</code>, all files and subdirs from all levels are loaded.
 * If <code>hier_levels is</code> greater than the real number of levels from the root directory, all files and subdirs
 * will be loaded.
 *
 * @author Samuel Martins
 * @date Sep 25th, 2015
 * @ingroup File
 *
 * @param file_entry The pathname from the directory of the CSV file to be read.
 * @param hier_levels Number of levels that will be loaded.
 * @param sort_pathnames It indicates if the pathnames should be sorted. 
 *
 * @return The array of iftFiles.
 *
 * @exception The file entry <code><b>file_entry</b></code> is neither a directory nor a CSV file.
 */
iftFileSet *iftLoadFileSetFromDirOrCSV(const char *file_entry, size_t hier_levels, bool sort_pathnames);


/**
 * @brief      Writes a File Set as a CSV, one pathname per row.
 *
 * @param[in]  fset  File set to be written.
 * @param[in]  path  Output CSV file.
 * 
 * @author Samuka
 * @date   Jan, 27, 2017
 */
void iftWriteFileSetAsCSV(const iftFileSet *fset, const char *path);


/**
 * @brief Checks if the files with the pathnames from a given File Set exist, raising an exception if
 * some file does not exist.
 * @author Samuel Martins
 * @date Apr 4, 2016
 * @ingroup File 
 */
void iftCheckFileSetPathnames(const iftFileSet *fset);

/**
 * @brief Destroys an iftFile from memory, <b>BUT IT DOES NOT DELETE IT FROM THE DISK</b>.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param file The iftFile to be destroyed.
 */
void iftDestroyFile(iftFile **file);


/**
 * @brief Sorts the file list by pathname
 * @author Thiago Vallin Spina
 * @ingroup File
 *
 * @param files The file list.
 */
void iftSortFileSet(iftFileSet *files);


/**
 * @brief Sorts the file list in the directory by pathname
 * @author Thiago Vallin Spina
 * @ingroup File
 *
 * @param files The file list.
 */
void iftSortDir(iftDir *files);

/**
 * @brief Copies a fileset
 * @author Peixinho
 * @param Original fileset
 * @return Fileset copy
 * @date Jun, 2016
 */
iftFileSet *iftCopyFileSet(const iftFileSet* fileset);

/**
 * @brief Allocates an iftFileSet structure but does NOT allocate the elements of the file array.
 * @author Thiago Vallin Spina
 *
 * @param n The number of files.
 * @return The allocated structure
 */
iftFileSet *iftCreateFileSet(int nfiles);


/**
 * @brief Destroys a File Array (iftFileSet) from memory, <b>BUT IT DOES NOT DELETE IT FROM THE DISK</b>.
 * @author Samuel Martins
 * @date Oct 6th, 2015
 * @ingroup File
 *
 * @param file The iftFile to be destroyed.
 */
void iftDestroyFileSet(iftFileSet **farr);


/**
 * @brief Destroys an iftDir.
 * @author Samuel Martins
 * @ingroup File
 *
 * Destroys the iftDir <b><code>dir</code></b> from memory, BUT IT DOES NOT DELETE IT FROM THE DISK.
 * If its sub-directories also were listed/loaded, the function also destroys them.
 *
 * @param dir The iftDir to be destroyed.
 */
void iftDestroyDir(iftDir **dir);


/**
 * @brief Prints the information of an iftFile.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param f The iftFile to be printed.
 *
 * @warning It does not print the information of sub-directories from an iftFile dir.
 */
void iftPrintFileInfo(iftFile *f);


/**
 * @brief Prints the information of an iftDir.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param dir The iftDir to be printed.
 *
 * @attention It does not print the information of its sub-directories.
 */
void iftPrintDirInfo(iftDir *dir);


/**
 * @brief Prints all files and sub-directories from an iftFile directory as a tree.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param dir The iftDir to be printed.
 */
void iftPrintDirAsTree(iftDir *dir);


/**
 * @brief Checks if the <pathname> is a file (not dir) and exists on the disk.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname to be checked.
 * @return True if it is a file, false otherwise.
 */
bool iftFileExists(const char *pathname);


/**
 * @brief Checks if the <b><code>pathname</code></b> is a directory on the disk.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname to be checked.
 * @return True if it is a directory, false otherwise.
 */
bool iftDirExists(const char *pathname);


/**
 * @brief Joins two pathnames. It automatically treats the '/' (separation char) of the end of <code><b>pathname1</code></b> 
 * and the beginning of the <code><b>pathname2</code></b>.
 * @author Samuel Martins
 * @date Dec 13rd, 2015
 * @ingroup File
 *
 * @param pathname1 The leftmost pathname to be joined.
 * @param pathname2 The rightmost pathname to be joined.
 * @return The joined pathname.
 */
char *iftJoinPathnames(const char *pathname1, const char *pathname2);


/**
 * @brief Gets the parent dir from a file or directory.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname of the file/directory.
 * @return The parent directory from <pathname>.
 *
 * @warning It does not check if the parent dir exists.
 * @warning Return the parent dir WITHOUT THE SLASH AT THE END: e.g: /home/samuel --> parent_dir = /home
 */
char *iftParentDir(const char *pathname);


/**
 * @brief Creates all directories in a hierarchy on the disk if they do not exist.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname of the directory to be created on disk.
 *
 * @warning If the directory already exists, do nothing.
 * @note http://www.thegeekstuff.com/2012/06/c-directory/
 */
void iftMakeDir(const char *dir_path);


/**
 * @brief Creates a valid temporary file on <b>DISK</b>.
 * @author Samuel Martins
 * @date Dec 11, 2015
 * @ingroup File
 *
 * @param prefix Prefix from the tmp file. If it is NULL, no prefix is added.
 * @param suffix Suffix from the tmp file. If it is NULL, no prefix is added.
 * @return The temporary file pathname.
 * @sa iftTempPathname(), iftMakeTempDir()
 */
char *iftMakeTempFile(const char *prefix, const char *suffix);


/**
 * @brief Creates a valid temporary pathname. If a prefix is given (!= NULL), it is added into the temp pathname.
 * @author Samuel Martins
 * @date Dec 11th, 2015
 * @ingroup File
 *
 * @param prefix Prefix from the temporary pathname. If it is NULL, no prefix is added.
 * @param suffix Suffix from the temporary pathname. If it is NULL, no suffix is added.
 * @return The temporary pathname.
 *
 * @attention No file or directory is created on disk.
 * @sa iftMakeTempFile(), iftMakeTempDir()
 */
char *iftMakeTempPathname(const char *prefix, const char *suffix);


/**
 * @brief Creates a temporary directory on <b>DISK</b>. If a prefix is given (!= NULL), it is added into the tmp dir.
 * @author Samuel Martins
 * @date Dec 11th, 2015
 * @ingroup File
 *
 * @param prefix Prefix from the tmp directory. If it is NULL, no prefix is added.
 * @return The temporary directory pathname.
 * @sa iftMakeTempFile(), iftMakeTempPathname()
 */
char *iftMakeTempDir(const char *prefix);


/**
 * @brief Removes the file @p <b>pathname</b> from the <b>DISK</b>.
 * @author Samuel Martins
 * @date Oct 7th, 2015
 * @ingroup File
 *
 * @param pathname The file to be removed from the DISK.
 *
 * @warning If the <code><b>pathname</code></b> is a directory or it does not exist, nothing will be done.
 */
void iftRemoveFile(const char *pathname);


/**
 * @brief Removes the directory <code><b>dir_pathname</code></b> from the <b>DISK</b>.
 * @author Samuel Martins
 * @date Oct 7th, 2015
 * @ingroup File
 *
 * @param dir_pathname The directory to be removed from the DISK.
 *
 * @warning If the <code><b>dir_pathname</code></b> is a file or it does not exist, nothing will be done.
 */
void iftRemoveDir(const char *dir_pathname);


/**
 * @brief Gets the absolute path from a pathname.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname The pathname.
 * @return The absolute pathname from <code><b>pathname</code></b>.
 */
char *iftAbsPathname(const char *pathname);


/**
 * @brief Gets the filename from a pathname. If <code><b>suffix</code></b> is != NULL and != "", the suffix will be removed from the filename.
 * @author Samuel Martins
 * @date Dec 11, 2015
 * @ingroup File
 *
 * E.g: \n
 * ("/home/batman/Pictures/avatar.png", NULL) = "avatar.png" \n
 * ("/home/batman/Pictures/avatar.png", .png) = "avatar"
 * ("/home/batman/Pictures", .png) = "/home/batman/Pictures"
 * 
 * @param pathname The Input Pathname.
 * @param suffix The suffix to be removed from the Oathaname. If it is NULL or "", no suffix is removed.
 * @return The Filename from <b>pathname</b>.
 *
 * @note If the given <code><b>suffix</code></b> is not a suffix from the pathname, nothind will be removed.
 * @attention If <code><b>suffix</code></b> is != NULL and != "", the suffix will be removed from the filename.
 */
char *iftFilename(const char *pathname, const char *suffix);


/**
 * @brief Gets the Basename from a Pathname.
 * @author Samuel Martins
 * @date Mar 21, 2016
 * @ingroup File
 *
 * E.g: \n
 * ("/home/batman/Pictures/avatar.png") = "/home/batman/Pictures/avatar" \n
 * ("my_image.zip", .png) = "my_image"
 * 
 * @param  pathname Input Pathname
 * @return          The Basename from the Input Pathname
 */
char *iftBasename(const char *pathname);


/**
 * @brief Gets the Extension of file <code><b>pathname</code></b>. It there is no extension, a blank string "" will be returned.
 * @author Samuel Martins
 * @date Dec 11, 2015
 * @ingroup String
 *
 * @param str Pathname from the file.
 * @return The file extension or a blank string "", if there is no extension.
 *
 * @note The extension is returned getting the strings after the last '.' (including it)
 * @warning Only files with extensions "*.tar.gz", "*.tar.bz", and "*.tar.bz2" has the correct extension returned.
 */
const char *iftFileExt(const char *pathname);


/**
 * @brief Checks if the Pathname is from a valid Image Pathname (*.scn, *.pgm, *.ppm, *.png).
 * @author Samuel Martins
 * @date Dec 9th, 2015
 * @ingroup File
 *
 * @param img_pathname The Image Pathname to be checked.
 * @return True if the image pathname is valid, false otherwise.
 *
 * @exception img_pathname does not exist or is a directory.
 */
bool iftIsImageFile(const char *img_pathname);


/**
 * @brief Checks if the File Set only have supported images (*.scn, *.pgm, *.ppm, *.png), raising an error if false.
 * @author Samuel Martins
 * @date Aug 12, 2016
 * @ingroup File
 *
 * @param img_paths File Set with Image's pathnames.
 * @exception At least a pathname is not from an Image.
 */
void iftValidateFileSet(const iftFileSet *img_files);



/**
 * @brief      Filters the files from the file set <dst> by matching the filenames of each pathname from <src> with <dst>.
 * 
 * For each pathname from <src>, this functions tries to match its filename with the filename of one pathname in <dst>. \n
 * If matched, the pathname from <dst> is copied. Otherwise an exception will be raised.
 *
 * @param[in]  src   Source file set whose filenames will be filter the destination file set.
 * @param[in]  dst   Destination file set which will be filtered
 * @return           Filtered file set.
 * 
 * @author Samuka
 * @date Jan 27, 2017
 * @ingroup File
 */
iftFileSet *iftFilterFileSetByFilename(const iftFileSet *src, const iftFileSet *dst);


/**
 * @brief      Filters the files by indexes. 
 *
 * @param[in]  fset  File set to be filtered.
 * @param      idxs  Indexes from the files to be chosen.
 *
 * @return     Filtered file set.
 * 
 * @author Samuka
 * @date Jan 30, 2017
 * @ingroup File
 */
iftFileSet *iftFilterFileSetByIndexes(const iftFileSet *fset, const iftIntArray *idxs);


/**
 * @brief Check if a file is empty.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param fp File pointer to the file.
 * @return True if it is a directory, false otherwise.
 */
bool iftIsFileContentEmpty(FILE *fp);


/**
 * @brief Checks if the Image Pathname has a valid suported Extension (scn, pgm, ppm, png).
 * @author Samuel Martins 
 * @date Mar 3, 2016
 * @ingroup File
 *  
 * @param  img_pathname Image Pathname to be checked.
 * @return              True if it's valid, false otherwise;
 *
 * @warning It DOESN'T check if it exists a file with this pathname. 
 */
bool iftIsImagePathnameValid(const char *img_pathname);


/**
 * @brief Returns the content of a FILE with pathname <code><b>pathname</code></b> in a string.
 * @author Samuel Martins
 * @ingroup File
 *
 * @param pathname Pathname from the FILE to be read.
 * @return The FILE content as a String.
 *
 * @note All the content of the file is stored in an only string.
 */
char *iftReadFileAsString(const char *pathname);


/**
 * @brief Finds a string in a File. If the string exists, it returns true. Otherwise, returns false.
 * @author Samuel Martins
 * @date Dec 8th, 2015
 * @ingroup File
 *
 * @param pathname Pathname from the FILE to be read. 
 * @param str String to be searched in the FILE. 
 * @return True if the string exists in the File. False, otherwise.
 *
 * @exception pathname does not exist or is not a File.
 */
bool iftFindInFile(const char *pathname, const char *str);


/**
 * @brief Adds Escape Characters to '(', ')', and '/' (linux, mac) or '\' (windows) to a Pathname
 *        (it should work on Linux, Mac and Windows).
 * @author Samuel Martins
 * @date Dec 2nd, 2015
 * @ingroup File
 *
 * E.g (Linux): pathname = "/usr/local/(abc)"
 *              escaped  = "\/usr\/local\/\(abc\)"
 *
 * @param pathname Pathname to be used.
 * @return The escaped pathname.
 *
 * @warning Function limitation: The function just works to pathnames with at most 1000 characteres to escape.
 */
char *iftAddEscapeChars(const char *pathname);


/**
 * @brief Creates a dictionary from a File Set, where for each, the pair (key, value) corresponds to the
 * File Key and File Pathname.
 *
 * See iftFileKey() for more information about File Keys.
 * 
 * @param  fset The File Set to be used.
 * @return      The resulting Dictionary of Files.
 */
iftDict *iftFileSetToDict(const iftFileSet *fset);


/**
 * @brief Get the file index, according to the path.
 * @author Peixinho
 *
 * The image index is stored as part of the image file name.
 * The IFT default nomenclature is <label>_<sample>_<suffix>.<ext>
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return The image index.
 * @sa iftValidDataSetFilenameConvention
 */
int iftImageSampleId(char *path);


/**
 * @brief Get the file label, according to the path.
 * @author Peixinho
 *
 * The image label is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return The image label.
 * @sa iftValidDataSetFilenameConvention
 */
int iftImageSampleLabel(char *path);

/**
 * @brief Get the file label and id, according to the path.
 * @author Deangeli
 *
 * The image label is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return A int pointer which the first position (index 0) stores the image label and the second one (index 1) stores
 * the sample id.
 * @sa iftValidDataSetFilenameConvention
 */
int* iftGetImageSampleInfo(char *path);

/**
 * @brief Zip the Content of a Dir.
 * @author Samuel Martins
 * @date Jun 9, 2016
 * @ingroup File
 * 
 * @param dir_path     Directory's pathname.
 * @param out_zip_path Pathname from the Output Zip file.
 */
void iftZipDirContent(const char *dir_path, const char *out_zip_path);


/**
 * @brief Unzip a zip file to a directory.
 * @author Samuel Martins
 * @date Jun 9, 2016
 * @ingroup File
 * 
 * @param zip_path Zip file to be unzipped.
 * @param out_dir  Output directory where the zip file will be unzipped
 */
void iftUnzipFile(const char *zip_path, const char *out_dir);


/**
 * @brief Get the file suffix, if the path follows IFT's convention.
 * @author Thiago Vallin Spina
 *
 * The image suffix is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param path Image path.
 * @return The image suffix.
 * @sa iftValidDataSetFilenameConvention
 */
char* iftImageSuffix(char* path);


/**
 * @brief Get the file labels, according to the filed path under that directory.
 * @author Peixinho
 *
 * The image label is stored as part of the image file name.
 *
 * Ex.: 000002_000015_LeftLung.pgm
 *
 * @param dir Directory containing the files.
 * @return The image labels.
 * @sa iftValidDataSetFilenameConvention
 */
int* iftImageLabels(iftDir *dir);


/**
 * @brief Checks if the path's filename conforms to IFT's default nomenclature.
 *
 * @author Thiago Vallin Spina
 * @date Mar 09, 2016
 *
 * The IFT default nomenclature is <label>_<sample>_<suffix>.<ext>
 *
 * Ex1.: 000002_000015_LeftLung.pgm
 *
 * The suffix may also be empty to conform with the old convention, or be separated by multiple '_' characters
 *
 * Ex2.: 000002_000015.pgm
 * Ex3.: 000002_000015_Right_Lung.pgm
 *
 * @param path The path name.
 *
 * @note The old convention 000002_0000015.pgm is also respected.
 */
bool iftValidDataSetFilenameConvention(const char *path);

/**
 * @brief Verifies if all files in a FileSet conform to the specified dataset convention.
 *
 * @author Thiago Vallin Spina
 * @date Mar 09, 2016
 *
 * @param files The input files
 */
bool iftValidDataSetFilenames(const iftFileSet *files);

/**
 * @brief Finds the image file by label and returns its index in the list.
 * @author Thiago Vallin Spina
 *
 * @param files File list.
 * @param label Image label.
 *
 * @return The index of the file in the file list or NIL if it is not present.
 */
int iftFindFileByLabel(const iftFileSet *files, int label);

/**
 * @brief Finds the image file by label and sample id and returns its index in the list.
 *
 * @author Thiago Vallin Spina
 * @date Mar 09, 2016
 *
 * @param files File list.
 * @param label Image label.
 *
 * @return The index of the file in the file list or NIL if it is not present.
 */
int iftFindFileByLabelAndSample(const iftFileSet *files, int label, int sample);

/**
 * @brief Finds the file by the path's basename'sand returns its index in the list.
 * @author Thiago Vallin Spina
 *
 * @param files File list.
 * @param prefix The path's basename.
 *
 * @return The index of the file in the file list or NIL if it is not present.
 */
int iftFindFileByBasename(const iftFileSet *files, const char *bname, const char *suffix);

/**
 * @brief Filters subset <code><b>full_set</code></b> by the basenames of the paths in <code><b>subset</code></b> after removing their suffixes. NOTE: All
 * directories are sorted in-place afterwards.
 * @author Thiago Vallin Spina
 *
 * @param subset The subset used as filtering mask.
 * @param full_set The set of original files.
 * @param suffix_subset The suffix of files in <code><b>subset</code></b>
 * @param suffix_full_set The suffix of files in <code><b>full_set</code></b>
 *
 * @return The filtered pathnames.
 */
iftFileSet *iftFilterFilesBySubset(iftFileSet *subset, iftFileSet *full_set, const char *suffix_subset, const char *suffix_full_set);

/*****************************************************************/

/**
 * @brief Converts a FileSet to a CSV file structure.
 *
 * @author Thiago Vallin Spina
 * @date Mar 9, 2016
 *
 * @param files The input file set.
 * @return A CSV file where each row has a path from the file set.
 */
iftCSV *iftFileSetToCSV(const iftFileSet *files);


/**
 * @brief Links two file sets by indexing them into a dict by their filenames.
 *
 * It tries to index all pathnames from file set1 with file set2. \n\n
 * 
 * Ex: \n
 * file_set1 = ["./orig/0001.scn", "./orig/0002.scn", "./orig/0003.scn"] \n
 * file_set2 = ["./labels/0001.scn", "./labels/0002.scn", "./labels/0003.scn"] \n
 * out_dict = {"./orig/0001.scn": "./labels/0001.scn", "./orig/0002.scn": "./labels/0002.scn",
 *             "./orig/0003.scn": "./labels/0003.scn"}
 * 
 * @param  file_set1 First file set.
 * @param  file_set2 Second file set.
 * @return           Output dict.
 */
iftDict *iftLinkFileSetsByFilename(const iftFileSet *file_set1, const iftFileSet *file_set2);


/**
 * @brief Runs a program.
 * @param prog Program name.
 * @param args Arguments of the program.
 *
 * @note args can be built on-the-fly.
 */
void iftRunProgram(const char *prog, const char *args, ...);


#ifdef __cplusplus
}
#endif

#endif
