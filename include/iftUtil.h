/**
 * @file
 * @brief General utilities.
 * @author Samuel Martins
 * @date Jan 13, 2016
 *
 * @note Programs:
 * * @ref iftRegexExamples.c = Examples and a quick explanation about Regex
 * * @ref iftTestRegex.c = Checks if an expression/string matches with a Regex
 * * @ref iftTestCmdLineParser1.c, @ref iftTestCmdLineParser2.c = Show how to parse input command line
 *   arguments from the keyboard
 */

#ifndef IFT_UTIL_H
#define IFT_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif


#include <regex.h>
#include "iftCommon.h"
#include "iftDict.h"
#include "iftString.h"

typedef enum {
    IFT_LINUX, IFT_MACOS, IFT_WIN32, IFT_WIN64
} iftOS;


/**
 * @brief An Option from the Command Line.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 * @warning At least one name must be defined. Otherwise, an error will be raised in iftCreateCmdLineParser().
 */
typedef struct ift_cmd_line_opt {
    /** Short option name. IT MUST HAVE THE PREFIX '-' */
    char short_name[4];
    /** Long option name. IT MUST HAVE THE PREFIX '--' */
    char long_name[64];
    /** Tells if the option has an argument. */
    bool has_arg;
    /** Datatype of the Argument (if the option have it). */
    iftCDataType arg_type;
    /** Tells if the option is required. */
    bool required;
    /** Description of the option (optional). It is used in the iftPrintUsage(). */
    char help[512];
} iftCmdLineOpt;


/**
 * @brief The Command Line Parser.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 */
typedef struct ift_cmd_line_parser {
    /** Program name */
    char program_name[64];
    /** Number of options. */
    int n_opts;
    /** An array with the command line options. */
    iftCmdLineOpt *opts;
    /** Description of the program usage. It is used in the iftPrintUsage(). */
    char description[IFT_STR_DEFAULT_SIZE];
} iftCmdLineParser;

/**
 * @brief Keep information of progress of a task
 * @author Peixinho
 * @date Apr, 2017
 */
typedef struct ift_progress_bar {
    float start;
    float end;
    float cur;
    int size;
} iftProgressBar;

/**
 * Updates the current position of a progressbar object.
 * @param bar iftProgressBar object
 * @param step current position
 * @date Apr, 2017
 * @author Peixinho
 */
void iftProgressBarUpdate(iftProgressBar* bar, int cur);

/**
 * @brief Creates a progress bar
 * @param start starting point
 * @param end ending point
 * @return iftProgressBar object
 * @author Peixinho
 * @date Apr, 2017
 */
iftProgressBar* iftCreateProgressBar(float start, float end);


/**
 * @brief Destroys a progress bar
 * @param pbar Pointer to progress bar object
 * @return iftProgressBar object
 * @author Peixinho
 * @date Apr, 2017
 */
void iftDestroyProgressBar(iftProgressBar** pbar);

/**
 * @brief Creates a progress bar with specified size
 * @param start starting point
 * @param end ending point
 * @param size progress bar print size
 * @return iftProgressBar object
 * @author Peixinho
 * @date Apr, 2017
 */
iftProgressBar* iftCreateProgressBarSize(float start, float end, int size);

/**
 * @brief Checks if a string matches with a regular expression.
 * @author Samuel Martins
 * @date Dec 11st, 2015
 * @ingroup Regex
 * @note See a complete demo in @ref iftRegexExamples.c and @ref iftTestRegex.c
 *
 * @param str The string to be checked.
 * @param regex_pattern The regular expression.
 * @return True/False.
 *
 * 
 * @exception str is NULL.
 * @exception regex_pattern is NULL.
 * @exception regex_pattern does not compile.
 */
bool iftRegexMatch(const char *str, const char *regex_pattern, ...);


/**
 * @brief Creates a Command Line Parser from a set of command line options.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 *
 * @param description Description of the program usage (optional). It is used in iftPrintUsage().
 * @param Number of options.
 * @param Array of Command Line Options.
 *
 * @exception Number of options is negative.
 * @exception Command Line Option with no name defined.
 * @exception Option Short name without the prefix '-'.
 * @exception Option Short name is "-h".
 * @exception Option Long name without the prefix '--'.
 * @exception Option Long name is "--help".
 */
iftCmdLineParser *iftCreateCmdLineParser(const char *description, int n_opts, iftCmdLineOpt cmd_line_opts[]);


/**
 * @brief Destroys a Command Line Parser.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 */
void iftDestroyCmdLineParser(iftCmdLineParser **parser);


/**
 * @brief Parses an Input Command Line.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 * 
 * @param argc Size of the array <code><b>argv</b></code>.
 * @param argv String array.
 * @param parser The command line parser.
 * @return The dictionary with the passed arguments.
 * 
 * @note If no options is passed or the options "-h", "--h", or "--help" are passed, the function is going to call the
 * usage printing function iftPrintUsage().
 * @note See a complete demo in demo/Miscellaneous/iftTestCmdLineParser*.c
 */
iftDict *iftParseCmdLine(int argc, const char *argv[], iftCmdLineParser *parser);


/**
 * @brief Prints the Usage Message from the command line options.
 * @author Samuel Martins
 * @date Feb 15, 2016
 * @ingroup CommandLineParser
 * 
 * @param parser Command Line Parser that contains all command line options.
 *
 * @note This function is called by the iftParseCmdLine().
 */
void iftPrintUsage(iftCmdLineParser *parser);


/**
 * @brief getline implementation to ensure that the function will work on Linux and Windows
 *
 * Same parameters as getline from stdio.h.
 *
 * @author Thiago Vallin Spina
 * @date Aug 29, 2016
 *
 * @param lineptr The pointer to the line to be read
 * @param n The size of the line.
 * @param stream The file stream.
 *
 * @return The number of bytes in the string (until the NULL terminating character)
 *
 * @note obtained from http://stackoverflow.com/questions/735126/are-there-alternate-implementations-of-gnu-getline-interface
 * and fixed since there was a bug due to realloc being able to change the pointer address of bufptr, in which case
 * the character pointer *p would lose the original reference. Modification by TVS.
 *
 */
size_t iftGetLine(char **lineptr, size_t *n, FILE *stream);

/**
 * @brief Get the Operational System. Operation Systems supported: IFT_LINUX, IFT_MACOS, IFT_WIN32, IFT_WIN64
 * @author Samuka Martins
 * @date Jul 18, 2017
 * @return The enumeration of the Operation System
 */
iftOS iftGetOS(void);

/**
 * @brief Set the environment to run ANTs programs (available in ift_dir/externals/ANTs)
 * @note Supported by Linux and Mac OS
 *
 * @author Samuka Martins
 * @date Jul 18, 2017
 */
void iftSetANTsEnvironment(void);


#ifdef __cplusplus
}
#endif

#endif //IFT_UTIL_H
