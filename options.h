/**
* options.h
*
* This header file contents the definition of the Option structure and
* all possible options for the project. All options are passive and do no actions.
* They are only used in other functions.
*
*/

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <stdlib.h>

/**
* This structure describes an option.
*/
struct Option
{
	const char * name;
	const char * option;
	const char * shortoption;
	const bool hasvalue;
	const char * defaultvalue;
	const char * value;
	const char * description;

	// static unsigned int s_uiOptionsCounter;

	Option
		(
		const char * _name,
		const char * _option,
		const char * _shoption = NULL,
		const bool _hasval = true,
		const char * _defval = NULL,
		const char * _descr = NULL
		) :
		name(_name),
		option(_option),
		shortoption(_shoption),
		hasvalue(_hasval),
		defaultvalue(_defval),
		description(_descr)
	{
		// ++s_uiOptionsCounter;
	}
};

// unsigned int Option::s_uiOptionsCounter = 0;

/**
* The array with all possible options.
*/
Option Options[] =
{
	Option("Help", "help", "?", false, NULL, "shows the help message"),
	Option("List devices", "listdevices", "lsdev", false, NULL, "shows all supported platforms and devices and quits"),
	Option("Platform", "platform", "pf", true, "0", "defines the platform, which will be used for running OpenCL kernel. By default the first possible platform will be used."),
	Option("Device", "device", "dev", true, "0", "defines the device, which will be used for running OpenCL kernel. By default the first possible device on the platform chosen will be used"),
	Option("Memory usage", "memusage", "mu", true, "0.95", "defines the portion of maximal allowed device memory allocation size, which will be used as a reference for calibrating parameters. By default it is 0.95"),
	Option("Workitems", "workitems", "wi", true, "0", "defines the number of workitems per each workgroup. If this parameter specified, the workgroups parameter and window size should also be specified. By default the number of workitems will be defined according to the device parameters"),
	Option("Workgroups", "workgroups", "wg", true, "0", "defines the number of workgroups run along a device. You can override the optimal number of workgroups by specifying this parameter. By default the number of workgroups will be defined according to the device parameters"),
	Option("Window size", "windowsize", "ws", true, "0", "defines the number of blocks in a window. You can override the optimal window size by specifying this parameter. By default the window size will be defined according to the device parameters"),
	Option("Reference sequence file", "reference", "ref", true, NULL, "the name of the file containing the reference sequence in FASTA format"),
	Option("Queries database file", "queries", "db", true, "", "the name of the file containing the queries sequences in FASTQ format"),
	Option("Output file", "output", "out", true, "", "the name of the file, to which the output information will be printed. By default all output are be printed directly to console"),
	Option("Matrix file", "matrix", "mat", true, "", "the name of the file containing the nucleotide scoring matrix, for example BLOSUM62, BLOSUM50 etc. By default the matrix with diagonal containing the match score value and the other values equal to mismatch fee value is used"),
	Option("Match score", "matchscore", "mscore", true, "2", "the match score value. This value is used for the constructing the default scoring matrix. In the case, when both, the scoring matrix file and the match score value, are specified, the scoring matrix has bigger priority. By default is 2"),
	Option("Mismatch fee", "mismatchfee", "mfee", true, "-1", "the mismatch fee value. This value is used for the constructing the default scoring matrix. In the case, when both, the scoring matrix file and the match score value, are specified, the scoring matrix has bigger priority. By default is -1"),
	Option("Insertion fee", "insertfee", "ifee", true, "-1", "the insertion fee value. This fee will be applied in the case of extending the reference sequence with the nucleotide from query sequence. By default is -1"),
	Option("Removing fee", "removefee", "rfee", true, "-1", "the removing fee value. This fee will be applied in the case of removing the nucleotides from the reference sequence. By default is -1"),
	Option("No path", "nopath", "np", false, NULL, "makes paths to be eliminated. This option should only be used for the fast computation of similarity values."),
	Option(NULL, NULL)
};

#endif /* OPTIONS_H_ */
