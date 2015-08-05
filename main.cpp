#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <errno.h>

#include "fbcl.h"
extern "C" {
#include "file_handler.h"
#include "utils.h"
}

#ifdef WINDOWS
#include <iostream>
#include <windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif



int main(int argc, char **argv){

	int i = 0, j = 0;
	cl_device_id device_cpu, device_gpu, device;
	int valid_cpu=1, valid_gpu=1;
	ALGORITHM type;
	int workers, parallel_jobs;
	
	char * query_file;
	char * subject_db;
	char * out_file;
	
	char *symbols = (char*) "ACGT";
	char sub_matrix[36] = {
		0, 0, 0, 0, 0, 0,
		0, 2, -1, -1, -1, 0,
		0, -1, 2, -1, -1, 0,
		0, -1, -1, 2, -1, 0,
		0, -1, -1, -1, 2, 0,
		0, 0, 0, 0, 0, 0
	};
	
	/*if(argc < 2){
		printf("Error. Usage: PolitoSW.exe <type> <work-items>.\nExiting...\n");
		return(-1);
	}*/
	
	type = LINEAR;
	workers = 32;
	parallel_jobs = 16;
	
#ifdef WINDOWS
	//query_file = (char*) "C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/queries_db_few_long.txt";
	//subject_db = (char*) "C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/subjects_db_few_long.txt";
	query_file = (char*) "C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/RF_query.txt";
	subject_db = (char*) "C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/RF00001_short.txt";

	out_file = (char *) "C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/alignment_scores.txt";
#else
	query_file = (char*) "queries_db.txt";
	subject_db = (char*) "subjects_db.txt";
	out_file = (char *) "alignment_scores.txt";
#endif

	for(i = 1; i < argc; i++){
		if((strcmp(argv[i], "-w") == 0)||(strcmp(argv[i], "--workers") == 0)){				// workers different from default (=32)
			sscanf(argv[i+1], "%i", &workers);	
		}
		else if((strcmp(argv[i], "-p") == 0)||(strcmp(argv[i], "--parallel_jobs") == 0)){	// parallel_jobs different from default (=16)
			sscanf(argv[i+1], "%i", &parallel_jobs); 
		}
		else if((strcmp(argv[i], "-t") == 0)||(strcmp(argv[i], "--type") == 0)){	// type different from default (=0)
			sscanf(argv[i+1], "%i", (int*)&type); 
		}
		else if((strcmp(argv[i], "-s") == 0)||(strcmp(argv[i], "--subject") == 0)){	// subject_db different from default (="subject_db.txt")
			subject_db = (char*) malloc(100);
			sscanf(argv[i+1], "%s", subject_db); 
		}
		else if((strcmp(argv[i], "-q") == 0)||(strcmp(argv[i], "--query") == 0)){	// subject_db different from default (="subject_db.txt")
			query_file = (char*) malloc(100);
			sscanf(argv[i+1], "%s", query_file); 
		}
		else if((strcmp(argv[i], "-o") == 0)||(strcmp(argv[i], "--output") == 0)){	// subject_db different from default (="subject_db.txt")
			out_file = (char*) malloc(100);
			sscanf(argv[i+1], "%s", out_file); 
		}else if((strcmp(argv[i], "-h") == 0)||(strcmp(argv[i], "--help") == 0)){
			printf("PolitoSW allows to align DNA sequences with the Smith-Watermann algorithm, comparing a set of queries with a subject\n");
			printf("database, both read from input files. It then allows to output the alignment scores into an output file.\n");
			printf("This application is optimized to run with OpenCL devices and parallelizes many alignments at once;\n");
			printf("it also allows to choose between different implementations of the algorithm, from the linear gap model to the affine one,\n");
			printf("and the Dynamic Gap Selector (DGS) optimized version.\n\n");
			printf("Usage:\n");
			printf("\t-h, --help:           Print this helpful help.\n");
			printf("\t-t, --type:           0=linear, 1=affine, 2=DGS. [default=0]\n");
			printf("\t-w, --worker:         Number of work-item per alignment. Must be a multiple of 32. [default=32]\n");
			printf("\t-p, --parallel_jobs:  Number of parallel alignments. [default=16]\n");
			printf("\t-s, --subject:        Subject filename from which read the reference sequences. [default=subjects.txt]\n");
			printf("\t-q, --query:          Query filename from which read the query sequences. [default=queries.txt\n");
			printf("\t-o, --output:         Output file into which save the alignment scores. [default=alignment_scores.txt]\n\n");
			return 1;
		}
		else{
			//printf("Error, wrong usage. Type -h or --help to see the help.\n");
			//return 1;
			;
		}
	}
	
	
	//sscanf(argv[1], "%i", (int*)&type); 
	//sscanf(argv[2], "%i", &workers);
	
	// if #WI < sub_matrix then sub_matrix will not be completely copied in local memory
	// thus providing wrong results.
	// See polito_locale.cl
	//if(workers < 36){ printf("Sorry, workers must be >= 36.\nExiting..\n"); return(-1); }
		
	//Get device
	if (!fbclGetFirstCPU(&device_cpu)){
		printf("ERR: valid CPU not found!\n");
		valid_cpu = 0;
	}
	if (!fbclGetFirstGPU(&device_gpu)){
		printf("ERR: valid GPU not found!\n");
		valid_gpu = 0;
	}

	if (valid_gpu){
		device = device_gpu;
	}
	else if (valid_cpu){
		device = device_cpu;
	}
	else{
		getchar();
		return(1);
	}	
	
	size_t wg_size;
	
	cl_int err;
	fbclGetMaxWorkGroupSize(device, &wg_size, &err, NULL);
	//printf("Max workgroup size: %ld\n", wg_size);
	
	if(workers > wg_size){
		printf("\n** Error: workers must be <= %ld.\nExiting...\n", wg_size);
		#ifdef WINDOWS
			getchar();
		#endif
		return -1;
	}

	printf("Running parameters:\n type=%d, workers=%d, parallel=%d, subject_db=%s queries=%s output_file=%s\n", (int)type, workers, parallel_jobs, subject_db, query_file, out_file);

#ifdef WINDOWS
	Sleep(1);
#else
	sleep(1);
#endif

	

	//device = device_cpu;
	//originalTestSW(device, workers, type);
	
	//testSW(device, workers, type);
	if((parallel_jobs == 0)&&(type == LINEAR)){
		clAlignSWSequentialChar(device, workers, type, query_file, subject_db, out_file, symbols, sub_matrix);
	}else{
		if(parallel_jobs == 0)
			parallel_jobs=1;
		//clAlignSWParallelSingle(device, workers, type, query_file, subject_db, out_file, symbols, sub_matrix);
		clAlignSWParallelChar(device, parallel_jobs, workers, type, query_file, subject_db, out_file, symbols, sub_matrix);
	}
	printf("\n");
	

#ifdef WINDOWS
	getchar();
#endif
	return(0);
}
/**********************************************************************************************************************************
 *											Original version of testSW																 *
 **********************************************************************************************************************************/

void originalTestSW(cl_device_id device, int worker, ALGORITHM type_algorithm){
	char *source_code;
	size_t source_size;
	//int type_algorithm = 1; //0: linear, 1: affine, 2: DGS_SW
	//int worker=32;
	clock_t start_timer, end_timer;
	double execution_time;

	//SW Test (by Wikipedia)
	char *symbols = (char*) "ACGT";
	// Amminoacid substitution matrix
	// DNA sub matrix, match= 2, unmatch= -1
	int sub_matrix[36] = {
		0, 0, 0, 0, 0, 0,
		0, 2, -1, -1, -1, 0,
		0, -1, 2, -1, -1, 0,
		0, -1, -1, 2, -1, 0,
		0, -1, -1, -1, 2, 0,
		0, 0, 0, 0, 0, 0
	};
	char *reference = (char*)"CATGTTTCCACTTACAGATCCTTCAAAAAGAGTGTTTCAAAACTGCTCTATGAAAAGGAATGTTCAACTCTGTGAGTTAAATAAAAGCATCAAAAAAAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTTCCATTTTCTCTATAAGCCTCAAAGCTGTCCAAATGTCCACTTGCAGATACTACAAAAAGAGTGTTTCAAAAGTGCTCAATGAAAAGGAATGTTCAGCTCTGTGAGTTAAATGCAAACATCACAAATAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGGGAAGATAATTCCGTGTCCAGCGAAGGCTTCAAAGCTTTCAAAATATCCACTTGCAAATTCTACAAAAAGAGTGTTTCAAAGCTGCTTTATCAAAAGAAAGTTTCAACTCTGTGAGTTGAATGTGCACATCACAAAGAAGTTTCTGAGAATGCCTTCAGTCTGGTTTTTATGTGAAGATATTCCCTTTTCCAACGAAAGCCTCGAAGCTGTCCAAATATCCACTTGTAAGTGCTGCAAAAAGAGTGTTTCAAAACTGCTACAGCAAAAGAAAGGTTTATCTCTGTGAGTTGAGTAGACACATCAAGAAGAAATTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTTCCTTTGTCACCATAGGCCTCCAAGCCCTCCAAATGTCCACTTGCAGATGCTACAAAAAGAGTGTTTCAAAACTGCTGTATGAAAAGAAATGCTCAAATCTGTGAGATAAATGCATACATCACAAAGAAGTCTTTGAGAATGCTTCTGTCTAGTTTTTATGTTAAGATATTTCCTATTTCACCATACGTCTCAACGCACACAAAATGTACACTTGCAGATGCTACAAAGAGAGTGTTTCAAAACTTGTAGATCAAAACAAGTGTTCAACTTTGTGAGTTGAGGACACACATCTGAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTCCCGTTTCCAGCGAAAGCCCCAAAACTATCCAAATATCCACTTGCACATTCTACAAAAAGAGTGTTTCAAATCTGCTCTATCAAAATAAAGGTTCAACTCTGTGAGTTGACTACACACATCACAAAGAAGTTTCTAAGAATGCTTCTGTCTGGTTTTTATGGGAAGATATTTCCTTTTTCAACATAGGCCTTGCAGCATCTACAAAAAGAGTTTTTCAAAACTCCTCTAAGAAAAGGAATGTTCAACTCCATGAGTTTAATGCAAAGATCACAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTAACCTGAAGACAGTTCCGTTTCCAGTGAAGGCTTGAAAGCTGTCCAAATATCCACATGCAAATTCTCAAAACGAGTGTTTAAAAGCTGCTCTATCACTAGAAAGTTTCACCTCTGTGAGCTGAATGCACACAGAAGTTTCTGAGAATGCTTCTGTCTGGTTTTTATGTGAAGATATTCCCGTTTCCAACCAAAGCCTCAAAGCTGTCCAAATATCCATTTGCAGATCCTACAGGGAGAGTGTTTCAAAACTGCTCTATAAAAAGAAAGGTTCTACTCTGTGAGTGGAGTACACACATCACAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAAATAGTTCGTTCTCCAAAATAGTCCTCAAAGCGCTCCAAATGTCCACTTGCAGATTCTACAAAAAGAGTGTTTCAAAACTGCTCTATGAAAAGGAATGTTCAACTCTGTGAGCTGAATGCAACCATCACAAAGAATTTTCTGAGAATGATTCTAATTATTATGTGAAGATATTCCCGTTTCGAACGAAGGCCTCAAAGCGCTCCAGATGTCCACTTGCAGATTCTACAAAAAGAGTGTTTCAAAACTGCCCTATGGAAAGGAGGATTCAACTCTATGAGTTGAATGCACACATCACAAATGGTTTTTGCGAATGCTTCTGTCTGCTTTTTATGTGAAGATATTTCCTTTTTCACCATAGGCCTCAAAGAAATCCAAATATCCACTTGCAGATACTACAGAAAGATTGTTTCAAAACTTCTGTTTCGACTCTGTGAGTTGAATGCACACATCATAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTTTGTGAAGATATTCCCGTTTCCAAGGAAAGCCTCAAAGCTGTCCTAATATCCACTTTTGAATTTTACAAGTAGAGTGTTTGAAAACTGCTCTGTCAAAAGAAAGGTTTAACTCTGTGTGTTGAGGGCTCACATCACAAAGAAGTTTCTGAGAATGCTTCCATCTATTTTTTATGTTAAGATATTTCCTTTTTCACCATAGGCCTCAAAGCACACCAAATGTCCACTTGCATATTGTGCAAAGAGTGTTTCAAAACTGTGCTATCAAAAGAAGCGTTCAACACTGTGAGTTGAACGCACACATCACAAAGAAGTTTCTGAGCATGCTTCTGTATAGTTTTTATGTGAAGATATTCCCGTATCCGAAGAAGGCCTCCAGGTGGTCTAAATATTTGCACATTCTACTAAAAGAGTGTTTAAAAGCTGCTGTATGATAAAGTATGTTCGAATTTGTGAGTTGAATGCACACATCACAAAGAAGTTTCTGAGAATGCTTCTGTCTAATATTTATGTGAAGATATTACGGTTTCCAACGAAGTCCTCAAAACCGTCCAAATATCCACTTGCAGATTTTATAAAAAGAGTGTTTCAAAGCTGCTCTATCAAAAGAAAGCTTCAACACTGTGAGTGGAATGCACACATCACAAAGTAGTTTCTGAAAATGCTGCTGTCTAGTTTTTATGTGAAGATATTCCCATTTCCAATGAAGGCCACAAAGCCGTCCAAATATCCACTTGCA";
	char *query = (char*) "ATTCATCTCACACATTTGAACATTTCTTTGATTGAATGGAAATACTCTTTTTGCATAATCTGCAAATGGATATTTGTGATATGTACTTTCATCTCACAGAGTTGAACCGGCCTATGGAGAAAAAGGAAATATCTTCAGATAAACTTCACATAAAAAATAGACAGAAGCATTGTGAGAAACCATGGTGAAATAGGAAATATCGTCACATAAAAAATAGTTTGGAAACACCCTTTTTGTAGAATCTACAAATGGCAGTTTTCTAAACAGACCTTTTGCAGAATCTACAAAACCGGTCTTTTTGTACAATCTGCAAAGGGATATTTCCAGTTTGGGGACAGTCTTTTTGTATAATGTGCAAAAGGCATATGGTGAAATAGGAAATGTCTTCACATAAAATGTGCTTTCACCACAAGGAGTTGAGCCTCTCTTTTGAGATTTGGACCTTGCTTTTCATTGAGCAGTTTGGCACTAGACAGAAGCATTTTGAGAAAACTCTTTGTGATGTGTTTACGGTTTGAGACCTATGGTGAAAAAGTAATATGTAATGTGTGCATTTATCTTACAGTGTTAAACCTTGGGGAACTTCTTTGTGATGTCTCCATTCATCTGACAGAGGAACTTCTTTGTGATGTCTCCATTCATCTGACAAAGGGATATTTTTGAGACATTTGAAGCCTAGAGTGATCTTTTGATTGAGAAGTATGGAAATGGTGGTCTTTTAGAATCTGGAAAGGGATATTTCTTAGCCCTTTGAGGATCTGGAAAGGGATATTTCTTAGCCCTTTGAGGCCTGGGATATTTCTTAGCCCTTTGAGGCCTATGGTGAGAGTTTGGAAACACTCTTTCTGTGACATCTGTAAATGGGGAAACACTCTTTCTGTGACATCTGTAAATGGATATCTTTCTTTTTATTGAGCAGTTTGGATACAGTCTTTTTATTTGTGAGCCCTTTATTGCCTATGGTGAAATAGGTGTGCATTCATCTCACATAGTTGAAACTTTCTTTGGCATTCATCTCACATAGTTGAAACTTTCTTGGGATTGCTCACAGAATTGAAACTATCCTTTGATTGAGGAGTTATGGTGAAAGAGAAATATCTTCCCATAAAAACTAGAATTCCAAGAAATTTCTTTGTGATGTGCCCATTCATCATGGATGTTTGGATTGCTTTGAGGCCTATGTTGAGAGAACGGTTTGATGCCTATGGTGAAAAAGAAATATCTCTGAAAAACTTCTCTGTGATGTGAGTATTCATGTCATTCTTTGTGATGTTTGCATTCATCTCACATATTTGAGTGTATTCGTCTCACAGAGTTGAACCTTTCTTTGCAGATGTGTGGTTTCATCTCACAGAGTTGAACCTTTCTCTGCAGAGGGATATTTGGGAGCCATTTAAGGACTATCTATGCTGAAAAAGGAAATATCTTCACATATAAGCT";

	//char *reference = "MWTRLAFCCWALALVSGWTNFSPMAPSLNFSFRLFPEASPGALGRLAVVPRSGEEEAVGSKVERLGRTFR";
	//char *query = "MWTRLTFCCWALALVSGWTNFQPMAPSLNFSFRLFPEASPGALGRLAVPPRSGEEEAVGSKVERLGRTFR";

	fbclVariadicGeneral(NULL, NULL, &worker, 3, worker, strlen(query) + 1, strlen(reference) + 1);

	int *buffer_reference, *buffer_query, size_reference, size_query;
	buffer_reference = get_buffer(reference, symbols, worker, &size_reference);
	buffer_query = get_buffer(query, symbols, worker, &size_query);

	BLOCK_AUX *buffer_aux_row, *buffer_aux_column;
	buffer_aux_row = get_buffer_aux(size_reference);
	buffer_aux_column = get_buffer_aux(size_query);

	// DNA sub matrix, match= 2, unmatch= -1
	int dna[] = {
		0, 0, 0, 0, 0, 0,
		0, 2, -1, -1, -1, 0,
		0, -1, 2, -1, -1, 0,
		0, -1, -1, 2, -1, 0,
		0, -1, -1, -1, 2, 0,
		0, 0, 0, 0, 0, 0
	};



	INFO info;
	info.query_lenght = size_query;
	info.reference_lenght = size_reference;
	info.block_size = worker;
	info.block_x = info.query_lenght / info.block_size;
	info.block_y = info.reference_lenght / info.block_size;
	info.penality1 = -1;
	info.penality2 = -10;
	info.symbol_number = strlen(symbols);

	BLOCK *kernel_array;
	kernel_array = (BLOCK*)malloc(sizeof(BLOCK)*info.block_size);
	for (int i = 0; i < info.block_size; i++){
		kernel_array[i].diagonal = 0;
		kernel_array[i].h_max = 0;
		kernel_array[i].left = 0;
		kernel_array[i].up = 0;
		kernel_array[i].e = 0;
		kernel_array[i].f = 0;
		kernel_array[i].e_diagonal = 0;
		kernel_array[i].f_diagonal = 0;
		kernel_array[i].e_aux = 0;
		kernel_array[i].f_aux = 0;
		kernel_array[i].gQry = info.penality1;
		kernel_array[i].gSbj = info.penality1;
	}

	//RUN SW TEST
	switch (type_algorithm)
	{
	case 0:
		#ifdef WINDOWS
			source_code = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw.cl", &source_size);
		#else
			source_code = fbclGetSourceFromFile("politosw.cl", &source_size);
		#endif
		printf("Algorith used for the alignment: SW-Linear \n");
		break;
	case 1:
		source_code = fbclGetSourceFromFile("politosw_affine.cl", &source_size);
		printf("Algorith used for the alignment: SW-Affine \n");
		break;
	case 2:
		source_code = fbclGetSourceFromFile("politosw_dgssw.cl", &source_size);
		printf("Algorith used for the alignment: SW-DGS \n");
		break;
	default:
		break;
	}

	start_timer = clock();
	fbclRun(device, info.block_size, source_code, source_size, (char*)"politosw", YES, &execution_time,
	(char*)"rrarar",
		sizeof(int) * info.reference_lenght, (void*)buffer_reference,
		sizeof(int) * info.query_lenght, (void*)buffer_query,
		sizeof(BLOCK_AUX) * info.reference_lenght, (void*)buffer_aux_row,
		sizeof(int) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)sub_matrix,
		sizeof(BLOCK) * info.block_size, (void*)kernel_array,
		sizeof(INFO), (void*)&info
		);
	end_timer = clock();

	free(source_code);

	int all_max = 0;

	printf("\n");
	for (int i = 0; i < info.block_size; i++){
		printf("Kernel %d: %d\n", i, kernel_array[i].h_max);
		if (all_max < kernel_array[i].h_max)
			all_max = kernel_array[i].h_max;
	}

	printf("\nResult: %d\n", all_max);
	printf("Took %lf", execution_time);
}



/**********************************************************************************************************************************
*											Modded version of testSW															  *
**********************************************************************************************************************************/

void testSW(cl_device_id _device, int worker, ALGORITHM type_algorithm){

	//SW Test (by Wikipedia)
	char *symbols = (char*) "ACGT";
	int sub_matrix[36] = {
		0, 0, 0, 0, 0, 0,
		0, 2, -1, -1, -1, 0,
		0, -1, 2, -1, -1, 0,
		0, -1, -1, 2, -1, 0,
		0, -1, -1, -1, 2, 0,
		0, 0, 0, 0, 0, 0
	};
	
	int i;
	char *source_code, *source_code_local;
	size_t source_size, source_size_local;
	int all_max;
	double execution_time, avg_time;
	double exec_time=0.0, tot_exec_time=0.0;
	
	#ifdef WINDOWS_PROFILE
		clock_t start_timer, end_timer, tot_start_timer, tot_end_timer;
		LARGE_INTEGER frequency, start, end, tot_start, tot_end;
		QueryPerformanceFrequency(&frequency);
	#else
		struct timeval tv_start, tv_end, tv_tot_start, tv_tot_end;
	#endif
	
	int ref_lines, ref_size, ref_max_length;
	int qry_lines, qry_size, qry_max_length;
	
	SEQUENCE *reference, *query;
#ifdef WINDOWS
	reference = readIndexedSequencesFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/subjects.txt" , &ref_lines, &ref_size, &ref_max_length);
	query = readIndexedSequencesFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/queries.txt", &qry_lines, &qry_size, &qry_max_length);
#else
	reference = readIndexedSequencesFromFile("subjects.txt", &ref_lines, &ref_size, &ref_max_length);
	query = readIndexedSequencesFromFile("queries.txt", &qry_lines, &qry_size, &qry_max_length);
#endif
	
	//printf("\nSubjects:\n");
	//printSequence(reference, ref_lines);
	
	//printf("\nQueries:\n");
	//printSequence(query, qry_lines);
	
	
	// testing: bypassing file reading (windows sucks)
	(*reference).seq = (char*) "CATGTTTCCACTTACAGATCCTTCAAAAAGAGTGTTTCAAAACTGCTCTATGAAAAGGAATGTTCAACTCTGTGAGTTAAATAAAAGCATCAAAAAAAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTTCCATTTTCTCTATAAGCCTCAAAGCTGTCCAAATGTCCACTTGCAGATACTACAAAAAGAGTGTTTCAAAAGTGCTCAATGAAAAGGAATGTTCAGCTCTGTGAGTTAAATGCAAACATCACAAATAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGGGAAGATAATTCCGTGTCCAGCGAAGGCTTCAAAGCTTTCAAAATATCCACTTGCAAATTCTACAAAAAGAGTGTTTCAAAGCTGCTTTATCAAAAGAAAGTTTCAACTCTGTGAGTTGAATGTGCACATCACAAAGAAGTTTCTGAGAATGCCTTCAGTCTGGTTTTTATGTGAAGATATTCCCTTTTCCAACGAAAGCCTCGAAGCTGTCCAAATATCCACTTGTAAGTGCTGCAAAAAGAGTGTTTCAAAACTGCTACAGCAAAAGAAAGGTTTATCTCTGTGAGTTGAGTAGACACATCAAGAAGAAATTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTTCCTTTGTCACCATAGGCCTCCAAGCCCTCCAAATGTCCACTTGCAGATGCTACAAAAAGAGTGTTTCAAAACTGCTGTATGAAAAGAAATGCTCAAATCTGTGAGATAAATGCATACATCACAAAGAAGTCTTTGAGAATGCTTCTGTCTAGTTTTTATGTTAAGATATTTCCTATTTCACCATACGTCTCAACGCACACAAAATGTACACTTGCAGATGCTACAAAGAGAGTGTTTCAAAACTTGTAGATCAAAACAAGTGTTCAACTTTGTGAGTTGAGGACACACATCTGAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTCCCGTTTCCAGCGAAAGCCCCAAAACTATCCAAATATCCACTTGCACATTCTACAAAAAGAGTGTTTCAAATCTGCTCTATCAAAATAAAGGTTCAACTCTGTGAGTTGACTACACACATCACAAAGAAGTTTCTAAGAATGCTTCTGTCTGGTTTTTATGGGAAGATATTTCCTTTTTCAACATAGGCCTTGCAGCATCTACAAAAAGAGTTTTTCAAAACTCCTCTAAGAAAAGGAATGTTCAACTCCATGAGTTTAATGCAAAGATCACAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTAACCTGAAGACAGTTCCGTTTCCAGTGAAGGCTTGAAAGCTGTCCAAATATCCACATGCAAATTCTCAAAACGAGTGTTTAAAAGCTGCTCTATCACTAGAAAGTTTCACCTCTGTGAGCTGAATGCACACAGAAGTTTCTGAGAATGCTTCTGTCTGGTTTTTATGTGAAGATATTCCCGTTTCCAACCAAAGCCTCAAAGCTGTCCAAATATCCATTTGCAGATCCTACAGGGAGAGTGTTTCAAAACTGCTCTATAAAAAGAAAGGTTCTACTCTGTGAGTGGAGTACACACATCACAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTATGTGAAAATAGTTCGTTCTCCAAAATAGTCCTCAAAGCGCTCCAAATGTCCACTTGCAGATTCTACAAAAAGAGTGTTTCAAAACTGCTCTATGAAAAGGAATGTTCAACTCTGTGAGCTGAATGCAACCATCACAAAGAATTTTCTGAGAATGATTCTAATTATTATGTGAAGATATTCCCGTTTCGAACGAAGGCCTCAAAGCGCTCCAGATGTCCACTTGCAGATTCTACAAAAAGAGTGTTTCAAAACTGCCCTATGGAAAGGAGGATTCAACTCTATGAGTTGAATGCACACATCACAAATGGTTTTTGCGAATGCTTCTGTCTGCTTTTTATGTGAAGATATTTCCTTTTTCACCATAGGCCTCAAAGAAATCCAAATATCCACTTGCAGATACTACAGAAAGATTGTTTCAAAACTTCTGTTTCGACTCTGTGAGTTGAATGCACACATCATAAAGAAGTTTCTGAGAATGCTTCTGTCTAGTTTTTTTGTGAAGATATTCCCGTTTCCAAGGAAAGCCTCAAAGCTGTCCTAATATCCACTTTTGAATTTTACAAGTAGAGTGTTTGAAAACTGCTCTGTCAAAAGAAAGGTTTAACTCTGTGTGTTGAGGGCTCACATCACAAAGAAGTTTCTGAGAATGCTTCCATCTATTTTTTATGTTAAGATATTTCCTTTTTCACCATAGGCCTCAAAGCACACCAAATGTCCACTTGCATATTGTGCAAAGAGTGTTTCAAAACTGTGCTATCAAAAGAAGCGTTCAACACTGTGAGTTGAACGCACACATCACAAAGAAGTTTCTGAGCATGCTTCTGTATAGTTTTTATGTGAAGATATTCCCGTATCCGAAGAAGGCCTCCAGGTGGTCTAAATATTTGCACATTCTACTAAAAGAGTGTTTAAAAGCTGCTGTATGATAAAGTATGTTCGAATTTGTGAGTTGAATGCACACATCACAAAGAAGTTTCTGAGAATGCTTCTGTCTAATATTTATGTGAAGATATTACGGTTTCCAACGAAGTCCTCAAAACCGTCCAAATATCCACTTGCAGATTTTATAAAAAGAGTGTTTCAAAGCTGCTCTATCAAAAGAAAGCTTCAACACTGTGAGTGGAATGCACACATCACAAAGTAGTTTCTGAAAATGCTGCTGTCTAGTTTTTATGTGAAGATATTCCCATTTCCAATGAAGGCCACAAAGCCGTCCAAATATCCACTTGCA";

	//(*reference).seq = "CATGTTTCCACTTACAGATCCTTCAAAAAGAGTGTTTCAAAACTGCTCTA";
	(*reference).index = 1;
	(*reference).seq_length = strlen( (*reference).seq);

	(*query).seq = (char*) "ATTCATCTCACACATTTGAACATTTCTTTGATTGAATGGAAATACTCTTTTTGCATAATCTGCAAATGGATATTTGTGATATGTACTTTCATCTCACAGAGTTGAACCGGCCTATGGAGAAAAAGGAAATATCTTCAGATAAACTTCACATAAAAAATAGACAGAAGCATTGTGAGAAACCATGGTGAAATAGGAAATATCGTCACATAAAAAATAGTTTGGAAACACCCTTTTTGTAGAATCTACAAATGGCAGTTTTCTAAACAGACCTTTTGCAGAATCTACAAAACCGGTCTTTTTGTACAATCTGCAAAGGGATATTTCCAGTTTGGGGACAGTCTTTTTGTATAATGTGCAAAAGGCATATGGTGAAATAGGAAATGTCTTCACATAAAATGTGCTTTCACCACAAGGAGTTGAGCCTCTCTTTTGAGATTTGGACCTTGCTTTTCATTGAGCAGTTTGGCACTAGACAGAAGCATTTTGAGAAAACTCTTTGTGATGTGTTTACGGTTTGAGACCTATGGTGAAAAAGTAATATGTAATGTGTGCATTTATCTTACAGTGTTAAACCTTGGGGAACTTCTTTGTGATGTCTCCATTCATCTGACAGAGGAACTTCTTTGTGATGTCTCCATTCATCTGACAAAGGGATATTTTTGAGACATTTGAAGCCTAGAGTGATCTTTTGATTGAGAAGTATGGAAATGGTGGTCTTTTAGAATCTGGAAAGGGATATTTCTTAGCCCTTTGAGGATCTGGAAAGGGATATTTCTTAGCCCTTTGAGGCCTGGGATATTTCTTAGCCCTTTGAGGCCTATGGTGAGAGTTTGGAAACACTCTTTCTGTGACATCTGTAAATGGGGAAACACTCTTTCTGTGACATCTGTAAATGGATATCTTTCTTTTTATTGAGCAGTTTGGATACAGTCTTTTTATTTGTGAGCCCTTTATTGCCTATGGTGAAATAGGTGTGCATTCATCTCACATAGTTGAAACTTTCTTTGGCATTCATCTCACATAGTTGAAACTTTCTTGGGATTGCTCACAGAATTGAAACTATCCTTTGATTGAGGAGTTATGGTGAAAGAGAAATATCTTCCCATAAAAACTAGAATTCCAAGAAATTTCTTTGTGATGTGCCCATTCATCATGGATGTTTGGATTGCTTTGAGGCCTATGTTGAGAGAACGGTTTGATGCCTATGGTGAAAAAGAAATATCTCTGAAAAACTTCTCTGTGATGTGAGTATTCATGTCATTCTTTGTGATGTTTGCATTCATCTCACATATTTGAGTGTATTCGTCTCACAGAGTTGAACCTTTCTTTGCAGATGTGTGGTTTCATCTCACAGAGTTGAACCTTTCTCTGCAGAGGGATATTTGGGAGCCATTTAAGGACTATCTATGCTGAAAAAGGAAATATCTTCACATATAAGCT";

	(*query).index = 1;
	(*query).seq_length = strlen((*query).seq);
	
	
	// Returns the minimun value (third parameter) between the arguments after '3'
	//fbclVariadicGeneral(NULL,NULL, &worker, 3, worker, strlen(query)+1, strlen(reference)+1);

	// Eseguo il padding a 16 sul minimo tra worker e la lunghezza di query e reference
	fbclVariadicGeneral16(NULL, NULL, &worker, 3, worker, strlen(query[0].seq) + 1, strlen(reference[0].seq) + 1);
	


	
	//printf("Reference: (original size=%ld) %s\n", strlen(reference[0].seq), reference[0].seq);
	//printf("Query: (original size=%ld), %s\n", strlen(query[0].seq), query[0].seq);
	
	

	

	// Initialization phase:
	// Get device LOCAL_MEMORY amount (hard constraint: algorithm must not exceed that value, otherwise kernel initialization fails)
	cl_int error = 0;
	cl_ulong device_local_memory = 0;
	fbclGetDeviceLocalMemory(_device, &device_local_memory, &error);
	if (error) printf("** Error getting device local memory. **\n");
	printf("\n\n\n[# INFO: Device Local memory: %ld bytes #]\n\n", device_local_memory);

	//RUN SW TEST
	printf("\n\n");
	

	switch (type_algorithm)
	{
		case 0:
		#ifdef WINDOWS
			source_code = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw.cl", &source_size);
			source_code_local = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw_local.cl", &source_size_local);
		#else
			source_code = fbclGetSourceFromFile("politosw.cl", &source_size);
			source_code_local = fbclGetSourceFromFile("politosw_local.cl", &source_size_local);

		#endif
			
			
			printf("Algorithm used for the alignment: SW-Linear + SW-Linear_Local \n");
			break;
		case 1:
			source_code = fbclGetSourceFromFile("politosw_affine.cl", &source_size);
			printf("Algorithm used for the alignment: SW-Affine \n");
			break;
		case 2:
			source_code = fbclGetSourceFromFile("politosw_dgssw.cl", &source_size);
			printf("Algorithm used for the alignment: SW-DGS \n");
			break;
		default:
			break;
	}

	

	
	// Kernel execution on all queries and subjects
	
	
	// [1]: separated single-run instances of the old kernel
	//	on the first query and subject ONLY
	// Each run is executed completely in a sequential manner
	
	
	int size_reference, size_query;
	int *buffer_reference = get_buffer(reference[0].seq, symbols, worker, &size_reference);
	int *buffer_query = get_buffer(query[0].seq, symbols, worker, &size_query);

	INFO info;
	info.query_lenght = size_query;
	info.reference_lenght = size_reference;
	info.block_size = worker;
	info.block_x = info.query_lenght / info.block_size;
	info.block_y = info.reference_lenght / info.block_size;
	info.penality1 = -1;
	info.penality2 = -10;
	info.symbol_number = strlen(symbols);

	BLOCK_AUX *buffer_aux_row = get_buffer_aux(size_reference);
	BLOCK * kernel_array = get_block(info.block_size, info.penality1);

	

	//printf("\n\nBuffer Reference: ");
	//for (i = 0; i < size_reference; i++) printf("%d", buffer_reference[i]);
	//printf("\nBuffer Query: ");
	//for (i = 0; i < size_query; i++) printf("%d", buffer_query[i]);
	
	
	//printf("\n\nQuery: %s\n", query[0].seq);
	//printf("Reference: %s\n", reference[0].seq);
	//printf("source: %s\n", source_code);

	printf("\n\nAll-sequential version for %d runs\n\n", KERNEL_RUNS);
	avg_time = 0.0;
	TOT_PROFILE_START();
	
	for(i=0; i<KERNEL_RUNS; i++){
		
		PROFILE_START();
		fbclRun(_device, info.block_size, source_code, source_size, (char*)"politosw", NO,	&execution_time,
			(char*)"rrarar",
			sizeof(int) * info.reference_lenght, (void*)buffer_reference,
			sizeof(int) * info.query_lenght, (void*)buffer_query,
			sizeof(BLOCK_AUX) * info.reference_lenght, (void*)buffer_aux_row,
			sizeof(int) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)sub_matrix,
			sizeof(BLOCK) * info.block_size, (void*)kernel_array,
			sizeof(INFO), (void*)&info
		);
			
		PROFILE_STOP();
		
		printf("\nRun %d:\n", i);
		all_max = getKernelMax(kernel_array, info.block_size);

		printf("  Result: %d\n", all_max);
		printf("  Running time: %lf ms\n", exec_time);
		printf("  [OCL profiling event] Took %lf ms\n", execution_time);
		
		avg_time += execution_time;
	}
	TOT_PROFILE_STOP();
	
	printf("\nTotal execution time for the %d runs: %lf ms\n", KERNEL_RUNS, tot_exec_time );
	printf(" Average kernel runtime for %d runs is %lf ms\n", KERNEL_RUNS, avg_time/KERNEL_RUNS);
	printf("--> Total runtime for all-sequential (original version) kernel runs: %lf ms\n", avg_time);
	

	

	
	
	free(reference);
	free(query);
	
	free(buffer_reference);
	free(buffer_query);
	
	free(buffer_aux_row);
	free(kernel_array);
	
	free(source_code);
	free(source_code_local);
	return;
}


/*********************************************************************************************************************************************************
 *						Sequentially align a set of references taken from query_file on a subject db taken from subject_db.								 *
 *										Prints the alignment results and the time taken on the out_file													 *
 *********************************************************************************************************************************************************/


void clAlignSWSequentialChar(cl_device_id device, int worker, ALGORITHM alg_type, char * query_file, char * subject_db, char * out_file, char * symbols, char * sub_matrix){
	int i;
	char *source_code, *source_code_local;
	size_t source_size, source_size_local;
	int all_max;
	double execution_time, avg_time;
	double exec_time = 0.0, tot_exec_time = 0.0;

#ifdef WINDOWS_PROFILE
	clock_t start_timer, end_timer, tot_start_timer, tot_end_timer;
	LARGE_INTEGER frequency, start, end, tot_start, tot_end;
	QueryPerformanceFrequency(&frequency);
#else
	struct timeval tv_start, tv_end, tv_tot_start, tv_tot_end;
#endif


	int ref_entries, ref_filesize, ref_max_length;
	int qry_entries, qry_filesize, qry_max_length;

	SEQUENCE *reference, *query;

	reference = readIndexedSequencesFromFile(subject_db, &ref_entries, &ref_filesize, &ref_max_length);
	query = readIndexedSequencesFromFile(query_file, &qry_entries, &qry_filesize, &qry_max_length);

	//sortSequencesByLength(reference, ref_entries, 'd');
	//sortSequencesByLength(query, qry_entries, 'd');


	printf("Reference lines: %d\n", ref_entries);
	printf("Query lines: %d\n", qry_entries);

	// Alignment results: ref_entries rows, qry_entries columns
	int ** scores = (int **)malloc(sizeof(int*) * ref_entries);
	float ** align_times = (float **)malloc(sizeof(float*) * ref_entries);
	for (i = 0; i < ref_entries; i++){
		scores[i] = (int *)calloc(qry_entries, sizeof(int));
		align_times[i] = (float *)calloc(qry_entries, sizeof(float));
	}

	// Getting local memory size
	cl_int error = 0;
	cl_ulong device_local_memory = 0;
	fbclGetDeviceLocalMemory(device, &device_local_memory, &error);
	if (error) printf("** Error getting device local memory. **\n");
	printf("\n\n\n[# INFO: Device Local memory: %ld bytes #]\n\n", device_local_memory);


	printf("\n\n");


	switch (alg_type)
	{
	case 0:
#ifdef WINDOWS
		//source_code_old = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw.cl", &source_size_old);
		source_code = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw_local.cl", &source_size);
#else
		//source_code_old = fbclGetSourceFromFile("politosw.cl", &source_size_old);
		source_code = fbclGetSourceFromFile("politosw_local.cl", &source_size);
#endif


		printf("Algorithm used for the alignment: SW-Linear \n");
		break;
	case 1:
		source_code = fbclGetSourceFromFile("politosw_affine.cl", &source_size);
		printf("Algorithm used for the alignment: SW-Affine \n");
		break;
	case 2:
		source_code = fbclGetSourceFromFile("politosw_dgssw.cl", &source_size);
		printf("Algorithm used for the alignment: SW-DGS \n");
		break;
	default:
		break;
	}

	
	cl_context _context;
	cl_command_queue _queue;
	cl_program _program;
	cl_kernel _kernel;

	fbclInitKernel(device, &_context, &_queue, source_code, source_size, &_program, &_kernel, (char*)"politosw", YES);



	/*************************************************		Execution	 *****************************************************/
	size_t size_char_query, size_char_reference;
	INFO info;
	info.penality1 = -1;
	info.penality2 = -10;
	info.symbol_number = strlen(symbols);

	BLOCK_AUX_LINEAR_SHORT *buffer_aux_row_linear;
	BLOCK_LINEAR_SHORT * kernel_array_linear;

	char * char_buffer_reference;
	char * char_buffer_query;

	// alignment specific

	int original_workers = worker;
	avg_time = 0.0;
	TOT_PROFILE_START();
	int j = 0;
	for (i = 0; i < ref_entries; i++){		// for all the subjects
		for (j = 0; j < qry_entries; j++){		// for all the queries
			
			worker = original_workers;

			// Setup alignment-specific data
			// worker is a multiple of 16
			fbclVariadicGeneral16(NULL, NULL, &worker, 3, worker, strlen(query[j].seq) + 1, strlen(reference[i].seq) + 1);

			
			size_t tmp1, tmp2;
			char_buffer_reference = getCharBuffer(reference[i].seq, reference[i].seq_length, symbols, worker, &size_char_reference);
			//size_char_reference = (int) tmp1;
			
			char_buffer_query = getCharBuffer(query[j].seq, query[j].seq_length, symbols, worker, &size_char_query);
			//size_char_query = (int) tmp2;
			
			info.block_size = worker;

			//char_buffer_{reference,query} are multiple of 16

			buffer_aux_row_linear = get_buffer_aux_linear(size_char_reference);
			kernel_array_linear = get_block_linear(worker);


			info.query_lenght = size_char_query;
			info.reference_lenght = size_char_reference;
			info.block_x = info.query_lenght / info.block_size;
			info.block_y = info.reference_lenght / info.block_size;
			
			
			// Compute internal parameters to efficiently copy from global to local memory
			short ref_complete_blocks = 0, ref_pad_block_copiers = 0, ref_single_pad_copiers = 0;
			divideEvenlyBlocks(info.reference_lenght, (short)worker, (short)16, &ref_complete_blocks, &ref_pad_block_copiers, &ref_single_pad_copiers);
			
			short qry_complete_blocks = 0, qry_pad_block_copiers = 0, qry_single_pad_copiers = 0;
			divideEvenlyBlocks(info.query_lenght, (short)worker, (short)16, &qry_complete_blocks, &qry_pad_block_copiers, &qry_single_pad_copiers);
			
			short sub_complete_copiers = 0, sub_complete_elements = 0, sub_pad_elements = 0;
			short sub_size = sizeof(char_sub_matrix) / sizeof(char_sub_matrix[0]);
			divideEvenly(sub_size, (short)worker, &sub_complete_copiers, &sub_complete_elements, &sub_pad_elements);
			


			CONTROL_SHORT ctrl;
			ctrl.ref_complete_blocks = ref_complete_blocks;
			ctrl.ref_pad_block_copiers = ref_pad_block_copiers;
			ctrl.ref_single_pad_copiers = ref_single_pad_copiers;
			
			ctrl.qry_complete_blocks = qry_complete_blocks;
			ctrl.qry_pad_block_copiers = qry_pad_block_copiers;
			ctrl.qry_single_pad_copiers = qry_single_pad_copiers;
			
			ctrl.sub_complete_copiers = sub_complete_copiers;
			ctrl.sub_complete_elements = sub_complete_elements;
			ctrl.sub_pad_elements = sub_pad_elements;


			printf("\n\n");
			printf("ref index: %d, qry index: %d, workers=%d\n", reference[i].index, query[j].index, worker);
			//printf("DEBUG: ref_l = %d, ref_complete_blocks = %d, ref_pad_block_copiers = %d, ref_single_pad_copiers = %d\n", info.reference_lenght, ctrl.ref_complete_blocks, ctrl.ref_pad_block_copiers, ctrl.ref_single_pad_copiers);
			//printf("DEBUG: qry_l = %d, qry_complete_blocks = %d, qry_pad_block_copiers = %d, qry_single_pad_copiers = %d\n", info.query_lenght, ctrl.qry_complete_blocks, ctrl.qry_pad_block_copiers, ctrl.qry_single_pad_copiers);
			//printf("DEBUG: sub_l = %d, sub_complete_copiers = %d, sub_complete_elements = %d, sub_pad_elements = %d\n\n", sub_size, ctrl.sub_complete_copiers, ctrl.sub_complete_elements, ctrl.sub_pad_elements);
			
			

			// Before starting the computation, we should compute the internal memory request for the run and evaluate if it fits
			// into local memory. If not, we should select a kernel with less local-memory parameters

			int local_memory_base_usage = size_char_reference + size_char_reference + sub_size;
			int local_memory_medium_usage = local_memory_base_usage + sizeof(BLOCK_LINEAR_SHORT)*info.block_size;
			int local_memory_max_usage = local_memory_medium_usage + sizeof(BLOCK_AUX_LINEAR_SHORT)*size_char_reference;

			printf("Local memory usage: base = %d, medium = %d, max = %d\n", local_memory_base_usage, local_memory_medium_usage, local_memory_max_usage);
			int score = 0;

			PROFILE_START();
			switch (alg_type){
			case LINEAR:
				fbclRunSequential(device, &_context, &_queue, &_program, &_kernel, info.block_size, source_code, source_size, (char*)"politosw", NO, &execution_time,
					(char*)"rrrrrwlllll",

					sizeof(char) * size_char_reference, (void*)char_buffer_reference,
					sizeof(char) * size_char_query, (void*)char_buffer_query,
					sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)char_sub_matrix,
					sizeof(INFO), (void*)&info,
					sizeof(CONTROL_SHORT), (void*)&ctrl,
					sizeof(int), (void*)&score,

					sizeof(char) * size_char_reference, NULL,
					sizeof(char) * size_char_query, NULL,
					sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), NULL,
					sizeof(BLOCK_AUX_LINEAR_SHORT) * info.reference_lenght, NULL,
					sizeof(BLOCK_LINEAR_SHORT) * info.block_size, NULL

					);
				break;
			default:
				printf("** Error: Invalid algorithm type (%d).\n", alg_type);
				break;
			}
			PROFILE_STOP();

			
			//all_max = getKernelMaxLinear(kernel_array_linear, info.block_size);
			//scores[i][j] = all_max;
			scores[i][j] = score;
			align_times[i][j] = execution_time;

			printf("  Result: %d\n", scores[i][j]);
			printf("  Running time: %lf ms\n", exec_time);
			printf("  [OCL profiling event] Took %lf ms\n", execution_time);

			avg_time += execution_time;

			free(buffer_aux_row_linear);
			free(kernel_array_linear);

			free(char_buffer_reference);
			free(char_buffer_query);


		}
	}

	TOT_PROFILE_STOP();



	printf("\n\n\nTotal execution time for the %d alignments: %lf ms\n", ref_entries*qry_entries, tot_exec_time);
	printf("--> Total runtime for all-sequential kernel alignments: %lf ms\n", avg_time);
	printf("(average on %d alignments: %lf ms)\n", ref_entries*qry_entries, avg_time / (ref_entries*qry_entries));



	/************************************************* Execution ended: destroy OpenCL structures and save output *****************************************************/
	fbclDestroy(&_context, &_program, &_kernel, &_queue);

	
	FILE * fp;
	if ((fp = fopen(out_file, "w")) == NULL){
		printf("Error opening ouput file. Not saving results\n");
	}
	char buff[MAX_LINE];
	sprintf(buff, "PolitoSW: aligned query file %s on subject database %s.\n\n", query_file, subject_db);
	fprintf(fp, "%s", buff);
	
	printf("\n\nFinal scores:\n");
	for (i = 0; i<ref_entries; i++){
		for (j = 0; j<qry_entries; j++){

			
			sprintf(buff, "[r:%3d (l=%4d) ** q:%3d (l=%4d)] = %6d - %10.5f ms\n", reference[i].index, reference[i].seq_length, query[j].index, query[j].seq_length, scores[i][j], align_times[i][j]);
			printf("%s", buff);

			fprintf(fp, "%s", buff);
		}
	}

	fclose(fp);

	free(source_code);
	free(source_code_local);
	
	free(reference);
	free(query);
	
	free(scores);
	free(align_times);
	return;
}


/*********************************************************************************************************************************************************
 *											Test for the parallel version of the SW alignment. Aligns only LINEAR_DIMENSION / workers couples at a time                      *
 *********************************************************************************************************************************************************/

int clAlignSWParallelSingle(cl_device_id device, int max_worker, ALGORITHM alg_type, char * query_file, char * subject_db, char * out_file, char * symbols, char * sub_matrix){
	int i = 0, j = 0;
	char *source_code;
	size_t source_size;
	int all_max;
	double execution_time, avg_time;
	double exec_time = 0.0, tot_exec_time = 0.0;

#ifdef WINDOWS_PROFILE
	clock_t start_timer, end_timer, tot_start_timer, tot_end_timer;
	LARGE_INTEGER frequency, start, end, tot_start, tot_end;
	QueryPerformanceFrequency(&frequency);
#else
	struct timeval tv_start, tv_end, tv_tot_start, tv_tot_end;
#endif


	size_t wg_size = 0;
	int err = 0;
	// Better device fitting:
	//int worker = WORKERS_PER_GROUP;
	int worker = max_worker;
	fbclGetMaxWorkGroupSize(device, &wg_size, &err, NULL);

	printf("max work group size: %ld\n", wg_size);

	int ref_entries, ref_filesize, ref_max_length;
	int qry_entries, qry_filesize, qry_max_length;

	SEQUENCE *reference, *query;

	reference = readIndexedSequencesFromFile(subject_db, &ref_entries, &ref_filesize, &ref_max_length);
	query = readIndexedSequencesFromFile(query_file, &qry_entries, &qry_filesize, &qry_max_length);
	
	//printf("Debug: reference sequence:\n");
	//printSequence(reference, ref_entries);
	//printf("\n\n\nDebug: query sequence:\n");
	//printSequence(query, qry_entries);

	if(reference == NULL){
		printf("Error reading reference file.\nExiting..\n");
		return -1;
	}
	if (query == NULL){
		printf("Error reading query file.\nExiting..\n");
		return -1;
	}

	//sortSequencesByLength(reference, ref_entries, 'd');		// 'd' = descending order
	//sortSequencesByLength(query, qry_entries, 'd');


	printf("Reference lines: %d\n", ref_entries);
	printf("Query lines: %d\n", qry_entries);
	
	
	// Alignment results: ref_entries rows, qry_entries columns
	int ** scores = (int **)malloc(sizeof(int*) * ref_entries);
	float ** align_times = (float **)malloc(sizeof(float*) * ref_entries);
	for (i = 0; i < ref_entries; i++){
		scores[i] = (int *)calloc(qry_entries, sizeof(int));
		align_times[i] = (float *)calloc(qry_entries, sizeof(float));
	}

	// Getting local memory size
	cl_int error = 0;
	cl_ulong device_local_memory = 0;
	fbclGetDeviceLocalMemory(device, &device_local_memory, &error);
	if (error) printf("** Error getting device local memory. **\n");
	printf("\n\n\n[# INFO: Device Local memory: %ld bytes #]\n\n", device_local_memory);


	printf("\n\n");


	switch (alg_type)
	{
	case 0:
#ifdef WINDOWS
		source_code = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw_parallelSingle.cl", &source_size);
#else
		source_code = fbclGetSourceFromFile("politosw_parallelSingle.cl", &source_size);
#endif
		printf("Algorithm used for the alignment: SW-Linear \n");
		break;
	case 1:
		source_code = fbclGetSourceFromFile("politosw_affine.cl", &source_size);
		printf("Algorithm used for the alignment: SW-Affine \n");
		break;
	case 2:
		source_code = fbclGetSourceFromFile("politosw_dgssw.cl", &source_size);
		printf("Algorithm used for the alignment: SW-DGS \n");
		break;
	default:
		break;
	}


	cl_context _context;
	cl_command_queue _queue;
	cl_program _program;
	cl_kernel _kernel;

	avg_time = 0.0;
	TOT_PROFILE_START();


	fbclInitKernel(device, &_context, &_queue, source_code, source_size, &_program, &_kernel, (char*)"politosw", YES);

	/************************************************* Parallel execution *****************************************************/

	
	GLOBAL_INFO info;
	info.penality1 = -1;
	info.penality2 = -10;
	info.symbol_number = strlen(symbols);
	info.block_size = worker;

	BLOCK_AUX_LINEAR_SHORT *buffer_aux_row_linear;
	BLOCK_LINEAR_SHORT * kernel_array_linear;

	size_t size_char_query = 1, size_char_reference = 1;
	char * char_buffer_reference = NULL;
	char * char_buffer_query = NULL;// = (char*)malloc(sizeof(char));

	int max_reference_length = 0, max_query_length = 0;
	int local_memory_base_usage = 0;
	int local_memory_medium_usage = 0;
	int local_memory_max_usage = 0;
	// alignment specific



	

	// Setup alignment-specific data
	
	// Parallel test: align the first query with the first 10 subjects (references)
	//int jobs = MIN(ref_entries, LINEAR_DIMENSION/WORKERS_PER_GROUP);
	int jobs = MIN(ref_entries, LINEAR_DIMENSION/max_worker);
	printf("-- jobs = %d\n", jobs);
	
	CONTROL_PARALLEL * ctrl = (CONTROL_PARALLEL *)calloc(jobs, sizeof(CONTROL_PARALLEL));

	if(ctrl == NULL){
		perror("Malloc error: ");
		return -1;
	}
	
	
	for(i = 0; i < jobs; i++){
		
		//printf("*********** i = %d ***********\n", i);
		unsigned long offset = 0;
		size_t tmp_size = 0;
		size_t final_dest_size = 0;

		char * tmp_buff = getCharBuffer(reference[i].seq, reference[i].seq_length, symbols, worker, &tmp_size);
		
		if(char_buffer_reference == NULL){
			//printf("Ref buffer is null: let's first fill it\n");
			char_buffer_reference = getCharBuffer(reference[i].seq, reference[i].seq_length, symbols, worker, &final_dest_size);
		}else{
			//printf("Realloc'ing ref global buffer\n");
			char_buffer_reference = reallocGlobalBuffer(char_buffer_reference, size_char_reference, tmp_buff, tmp_size, &final_dest_size, &offset);
		}
		//perror("Malloc error: ");
		//printf("tmp: workers=%d, size=%ld, buffer: \n", worker, tmp_size);
		//for(j=0; j<tmp_size; j++){
		//	printf("%d", tmp_buff[j]);
		//}
		//printf("\nref: size=%ld, buffer:\n", final_dest_size);
		//for(j=0; j<final_dest_size; j++){
		//	printf("%d", char_buffer_reference[j]);
		//}
		//printf("\n");
		size_char_reference = final_dest_size;
		
		if(tmp_size > max_reference_length) max_reference_length = tmp_size;
		
		//printf("reference_length[%d] = %ld\n", i, tmp_size);
		
		
		ctrl[i].reference_length = tmp_size;
		ctrl[i].block_y = ctrl[i].reference_length / info.block_size;
		ctrl[i].ref_id = reference[i].index;
		ctrl[i].ref_offset = offset;
		//printf("ref %d offset = %ld\n", i, ctrl[i].ref_offset);
		
		if((offset % 16) != 0 ){
			printf("ERROR: offset/16=%lf\n", (double)offset/16);
			return(-1);
		}
		
		/*
		printf("**DEBUG:\ni=%d:\n seq=%s\n index=%d\n", i, reference[i].seq, reference[i].index);
		printf("\nnew reference[%d]: length=%d, offset=%ld\nseq= ", i, ctrl[i].reference_lenght, ctrl[i].ref_offset);
		for (j = 0; j < ctrl[i].reference_lenght; j++) printf("%d", tmp_buff[j]);

		printf("\nchar_buffer_reference: length=%ld\nseq= ", size_char_reference);
		for (j = 0; j < size_char_reference; j++){
			if(char_buffer_reference[j] == 0)
				printf("0");
			else
				printf("%d", char_buffer_reference[j]);
		}
		printf("\nreference at offset %ld:\n", ctrl[i].ref_offset);
		for(j=0; j<ctrl[i].reference_lenght; j++) printf("%d", char_buffer_reference[ctrl[i].ref_offset + j]);
		
		//printf("\n\n");
		 */

		free(tmp_buff);
		tmp_buff = NULL;
		offset = 0;
		tmp_buff = getCharBuffer(query[0].seq, query[0].seq_length, symbols, worker, &tmp_size);
		if(char_buffer_query == NULL){
			char_buffer_query = getCharBuffer(query[0].seq, query[0].seq_length, symbols, worker, &final_dest_size);
		}else{
			//printf("Realloc'ing qry global buffer\n");
			char_buffer_query = reallocGlobalBuffer(char_buffer_query, size_char_query, tmp_buff, tmp_size, &final_dest_size, &offset);
		}
		size_char_query = final_dest_size;
		
		if(tmp_size > max_query_length) max_query_length = tmp_size;
		
		
		free(tmp_buff);
		tmp_buff = NULL;
		
		ctrl[i].query_length = tmp_size;
		ctrl[i].block_x = ctrl[i].query_length / info.block_size;
		ctrl[i].qry_id = query[0].index;
		ctrl[i].qry_offset = (unsigned long)offset;
		//printf("qry %d offset = %ld\n", i, ctrl[i].qry_offset);
		if((offset % 16) != 0 ){
			printf("ERROR: offset/16=%lf\n", (double)offset/16);
			return(-1);
		}

		// Compute internal parameters to efficiently copy from global to local memory
		short ref_complete_blocks = 0, ref_pad_block_copiers = 0, ref_single_pad_copiers = 0;
		divideEvenlyBlocks(ctrl[i].reference_length, (short)worker, (short)16, &ref_complete_blocks, &ref_pad_block_copiers, &ref_single_pad_copiers);
		ctrl[i].ref_complete_blocks = ref_complete_blocks;
		ctrl[i].ref_pad_block_copiers = ref_pad_block_copiers;
		ctrl[i].ref_single_pad_copiers = ref_single_pad_copiers;

		short qry_complete_blocks = 0, qry_pad_block_copiers = 0, qry_single_pad_copiers = 0;
		divideEvenlyBlocks(ctrl[i].query_length, (short)worker, (short)16, &qry_complete_blocks, &qry_pad_block_copiers, &qry_single_pad_copiers);
		ctrl[i].qry_complete_blocks = qry_complete_blocks;
		ctrl[i].qry_pad_block_copiers = qry_pad_block_copiers;
		ctrl[i].qry_single_pad_copiers = qry_single_pad_copiers;

		short sub_complete_copiers = 0, sub_complete_elements = 0, sub_pad_elements = 0;
		short sub_size = sizeof(char_sub_matrix) / sizeof(char_sub_matrix[0]);
		divideEvenly(sub_size, (short)worker, &sub_complete_copiers, &sub_complete_elements, &sub_pad_elements);
		ctrl[i].sub_complete_copiers = sub_complete_copiers;
		ctrl[i].sub_complete_elements = sub_complete_elements;
		ctrl[i].sub_pad_elements = sub_pad_elements;		
		
		/*
		printf("\n\n");
		printf("ref index: %d, qry index: %d, workers=%d\n", reference[i].index, query[0].index, worker);
		printf("DEBUG: ref_l = %d, ref_complete_blocks = %d, ref_pad_block_copiers = %d, ref_single_pad_copiers = %d, ref_offset = %ld\n", ctrl[i].reference_lenght, ctrl[i].ref_complete_blocks, ctrl[i].ref_pad_block_copiers, ctrl[i].ref_single_pad_copiers, ctrl[i].ref_offset);
		printf("DEBUG: qry_l = %d, qry_complete_blocks = %d, qry_pad_block_copiers = %d, qry_single_pad_copiers = %d, qry_offset = %ld\n", ctrl[i].query_lenght, ctrl[i].qry_complete_blocks, ctrl[i].qry_pad_block_copiers, ctrl[i].qry_single_pad_copiers, ctrl[i].qry_offset);
		printf("DEBUG: sub_l = %d, sub_complete_copiers = %d, sub_complete_elements = %d, sub_pad_elements = %d\n\n", sub_size, ctrl[i].sub_complete_copiers, ctrl[i].sub_complete_elements, ctrl[i].sub_pad_elements);
		*/
		
		local_memory_base_usage = reference[i].seq_length + query[0].seq_length + sub_size;
		local_memory_medium_usage = local_memory_base_usage + sizeof(BLOCK_LINEAR_SHORT)*info.block_size;
		local_memory_max_usage = local_memory_medium_usage + sizeof(BLOCK_AUX_LINEAR_SHORT)*max_reference_length;

		//printf("Local memory usage: base = %d, medium = %d, max = %d\n\n\n\n", local_memory_base_usage, local_memory_medium_usage, local_memory_max_usage);
		
	}
	//printf("\n\nDEBUG: max_ref_length = %d\n\n", max_reference_length);

	//buffer_aux_row_linear = get_buffer_aux_linear(max_reference_length);
	kernel_array_linear = get_block_linear(worker*jobs);

	printf("\n\nRunning kernel..\n");
	
	/*
	printf("Debug informations:\n");
	for(i=0; i<jobs; i++){
		printf("\n\n\n**i=%d:\nctrl[%d].reference_length: %d\n", i, i, ctrl[i].reference_lenght);
		printf("*ctrl[%d].block_y: %d\n", i, ctrl[i].block_y);
		
		printf("*ctrl[%d].query_length: %d\n", i, ctrl[i].query_lenght);
		printf("*ctrl[%d].block_x: %d\n", i, ctrl[i].block_x);
		
		printf("*ctrl[%d].ref_id: %ld\n", i, ctrl[i].ref_id);
		printf("****ctrl[%d].ref_offset: %ld\n", i, ctrl[i].ref_offset);
		printf("*ctrl[%d].ref_complete_blocks: %d\n", i, ctrl[i].ref_complete_blocks);
		printf("*ctrl[%d].ref_pad_block_copiers: %d\n", i, ctrl[i].ref_pad_block_copiers);
		printf("*ctrl[%d].ref_single_pad_copiers: %d\n", i, ctrl[i].ref_single_pad_copiers);
		printf("*ctrl[%d].ref_pad_block_copiers: %d\n", i, ctrl[i].ref_pad_block_copiers);
		
		printf("*ctrl[%d].qry_id: %ld\n", i, ctrl[i].qry_id);
		printf("****ctrl[%d].qry_offset: %ld\n", i, ctrl[i].qry_offset);
		printf("*ctrl[%d].qry_complete_blocks: %d\n", i, ctrl[i].qry_complete_blocks);
		printf("*ctrl[%d].qry_pad_block_copiers: %d\n", i, ctrl[i].qry_pad_block_copiers);
		printf("*ctrl[%d].qry_single_pad_copiers: %d\n", i, ctrl[i].qry_single_pad_copiers);
		printf("*ctrl[%d].qry_pad_block_copiers: %d\n", i, ctrl[i].qry_pad_block_copiers);
		
		printf("*ctrl[%d].sub_complete_copiers: %d\n", i, ctrl[i].sub_complete_copiers);
		printf("*ctrl[%d].sub_complete_elements: %d\n", i, ctrl[i].sub_complete_elements);
		printf("*ctrl[%d].sub_pad_elements: %d\n", i, ctrl[i].sub_pad_elements);
		
		printf("*ctrl[%d].score: %d\n", i, ctrl[i].score);
	}*/
	
	
	//CL_SEQUENCE * s = (CL_SEQUENCE *) calloc(jobs, sizeof(CL_SEQUENCE));
	int * sc = (int*) calloc(jobs, sizeof(int));
	//cl_ulong * sc1 = (cl_ulong*) calloc(jobs*info.block_size, sizeof(cl_ulong));
	int * lengths = (int*) calloc(jobs, sizeof(int));
	err = 0;
	execution_time = 0.0;

	printf("******************************** DEBUG ********************************\n");
	printf("Dimensions::\n  size_char_reference dimension=%ld\n  size_char_query dimension=%ld\n", sizeof(char) * size_char_reference, sizeof(char) * size_char_query);
	printf("  char_sub_matrix=%ld\n  info=%ld\n  control=%ld\n  sc=%ld\n", sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), sizeof(GLOBAL_INFO), sizeof(CONTROL_PARALLEL)*jobs, sizeof(cl_int)*jobs);
	printf("Sum of global sizes=%ld\n\n", sizeof(char)*size_char_reference + sizeof(char)*size_char_query + sizeof(char)*(info.symbol_number + 2)*(info.symbol_number + 2) + sizeof(GLOBAL_INFO) + sizeof(CONTROL_PARALLEL)*jobs + sizeof(cl_int)*jobs);
	
	printf("Local sizes::\n  max_reference_length=%ld\n  max_query_length=%ld\n  sub_matrix=%ld\n", sizeof(char) * max_reference_length, sizeof(char) * max_query_length, sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2));
	printf("  kernel_aux=%ld\n  kernel_matrix=%ld\n\n", sizeof(BLOCK_AUX_LINEAR_SHORT) * max_reference_length, sizeof(BLOCK_LINEAR_SHORT) * info.block_size);
	printf("Sum of local sizes=%ld\n\n\n", sizeof(char)*max_reference_length + sizeof(char)*max_query_length + sizeof(char)*(info.symbol_number + 2)*(info.symbol_number + 2) + sizeof(BLOCK_AUX_LINEAR_SHORT)* max_reference_length + sizeof(BLOCK_LINEAR_SHORT)*info.block_size);
	
	
	
	printf("****************************** RUN KERNEL *****************************\n");
	PROFILE_START();
	switch (alg_type){
	case LINEAR:
		err = fbclRunParallelSingle(device, &_context, &_queue, &_kernel, info.block_size, (char*)"politosw", NO, &execution_time,
				(char*)"rrrraalllll",
				
				// Global Parameters
				sizeof(char) * size_char_reference, (void*)char_buffer_reference,
				sizeof(char) * size_char_query, (void*)char_buffer_query,
				sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)char_sub_matrix,
				sizeof(GLOBAL_INFO), (void*)&info,
				
				sizeof(CONTROL_PARALLEL)*jobs, (void*)ctrl,
				sizeof(cl_int)*jobs, (void*)sc,
				
				// Local Parameters
				sizeof(char) * max_reference_length, NULL,
				sizeof(char) * max_query_length, NULL,
				sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), NULL,
				sizeof(BLOCK_AUX_LINEAR_SHORT) * max_reference_length, NULL,
				sizeof(BLOCK_LINEAR_SHORT) * info.block_size, NULL
				
			/*Global parameters
				sizeof(char) * size_char_reference, (void*)char_buffer_reference,
				sizeof(char) * size_char_query, (void*)char_buffer_query,
				sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)char_sub_matrix,
				sizeof(GLOBAL_INFO), (void*)&info,
				sizeof(CONTROL_PARALLEL)*jobs, (void*)ctrl,
			
			// Local parameters
				sizeof(char) * max_reference_length, NULL,
				sizeof(char) * max_query_length, NULL,
				sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), NULL,
				sizeof(BLOCK_AUX_LINEAR_SHORT) * max_reference_length, NULL,
				sizeof(BLOCK_LINEAR_SHORT) * info.block_size, NULL*/

			);
		
		break;
	default:
		printf("** Error: Invalid algorithm type (%d).\n", alg_type);
		break;
	}
	PROFILE_STOP();
	
	printf("\nKernel finished! err = %d\n", err);
	
	/*
	if(err){
		printf("Error with fbclRunParallelSingle..\nFreeing memory and exiting\n");
		
		free(buffer_aux_row_linear);
		free(kernel_array_linear);

		free(char_buffer_reference);
		free(char_buffer_query); 
		fbclDestroy(&_context, &_program, &_kernel, &_queue);

		
		free(source_code);
		printf("source\n");
		free(reference);
		printf("ref\n");
		free(query);
		printf("qry\n");
		free(scores);
		printf("scores\n");
		free(align_times);
		printf("align\n");
		return(-1);
	}*/

	//all_max = getKernelMaxLinear(kernel_array_linear, info.block_size);
	//scores[i][j] = all_max;
	
	// In clAlignSWParallelSingle we evaluate the first column only
	for(i=0; i < jobs; i++){
		//scores[i][0] = sc[i];
		scores[i][0] = ctrl[i].score;
		
		printf("  ** Result %d: %d **\n", i, ctrl[i].score);
		//printf("  ** Length %d: %d **\n", i, lengths[i]);
		printf("  Running time: %lf ms\n", exec_time);
		printf("  [OCL profiling event] Took %lf ms\n\n\n\n", execution_time);
	}

	align_times[0][0] = execution_time;
	
	TOT_PROFILE_STOP();



	//printf("\n\n##Max group id: %d\n", ctrl[0].score);
	/*
	for(i=0; i<jobs; i++){
		for(j=0; j<worker; j++){
			//printf("Kernel[%d][%d]: %d\n", i, j, kernel_array_linear[i*worker + j].h_max);
			printf("Kernel[%d][%d]: %ld\n", i, j, sc1[i*worker + j]);
		}
		printf("\n\n\n");
	}*/
	
	/*
	cl_int ii, jj;
	for(ii=0; ii<jobs; ii++){
		cl_int l = s[ii].seq_length;
		printf("Sequence %d: length = %d, l=%d\n", ii, l, ctrl[ii].reference_lenght);
		for(jj=0; jj<ctrl[ii].reference_lenght; jj++){
			printf("%d", s[ii].seq[jj]);
		}
		printf("\n\n");
	}*/
	
	//printf("%ld\n", sizeof(CL_SEQUENCE));


	/************************************************* Execution ended: destroy OpenCL structures and save output *****************************************************/
	
	fbclDestroy(&_context, &_program, &_kernel, &_queue);

	/*
	
	FILE * fp;
	if ((fp = fopen(out_file, "w")) == NULL){
		printf("Error opening ouput file. Not saving results\n");
	}

	printf("\n\nFinal scores:\n");
	for (i = 0; i<ref_entries; i++){
		for (j = 0; j<qry_entries; j++){

			char buff[MAX_LINE];
			sprintf(buff, "[r:%3d (l=%4d) ** q:%3d (l=%4d)] = %6d - %10.5f ms\n", reference[i].index, reference[i].seq_length, query[j].index, query[j].seq_length, scores[i][j], align_times[i][j]);
			//printf("%s", buff);

			fprintf(fp, buff);
		}
	}
	printf("\n");
	fclose(fp);	
	*/



	//free(buffer_aux_row_linear);
	//free(kernel_array_linear);
	
	free(ctrl);
	free(sc);
	
	free(char_buffer_reference);
	free(char_buffer_query);
	
	free(source_code);
	for(i=0; i<ref_entries; i++) free(reference[i].seq);
	free(reference);
	for(i=0; i<qry_entries; i++) free(query[i].seq);
	free(query);
	
	for (i = 0; i < ref_entries; i++){
		free(scores[i]);
		free(align_times[i]);
	}
	free(scores);
	free(align_times);

	return 0;
	
}



/*********************************************************************************************************************************************************
*					Align in a parallelized way a set of references taken from query_file on a subject db taken from subject_db.						 *
*										Prints the alignment results and the time taken on the out_file													 *
*********************************************************************************************************************************************************/

int clAlignSWParallelChar(cl_device_id device, int parallel_jobs, int max_worker, ALGORITHM alg_type, char * query_file, char * subject_db, char * out_file, char * symbols, char * sub_matrix){
	int i = 0, j = 0, h = 0;
	char *source_code;
	size_t source_size;
	int all_max;
	double execution_time, avg_time;
	double exec_time = 0.0, tot_exec_time = 0.0;

#ifdef WINDOWS_PROFILE
	clock_t start_timer, end_timer, tot_start_timer, tot_end_timer;
	LARGE_INTEGER frequency, start, end, tot_start, tot_end;
	QueryPerformanceFrequency(&frequency);
#else
	struct timeval tv_start, tv_end, tv_tot_start, tv_tot_end;
#endif


	size_t wg_size = 0;
	int err = 0;
	// Better device fitting:
	//int worker = WORKERS_PER_GROUP;
	int worker = max_worker;
	fbclGetMaxWorkGroupSize(device, &wg_size, &err, NULL);

	printf("max work group size: %ld\n", wg_size);

	int ref_entries, ref_filesize, ref_max_length;
	int qry_entries, qry_filesize, qry_max_length;

	SEQUENCE *reference, *query;

	reference = readIndexedSequencesFromFile(subject_db, &ref_entries, &ref_filesize, &ref_max_length);
	query = readIndexedSequencesFromFile(query_file, &qry_entries, &qry_filesize, &qry_max_length);
	
	//printf("Debug: reference sequence:\n");
	//printSequence(reference, ref_entries);
	//printf("\n\n\nDebug: query sequence:\n");
	//printSequence(query, qry_entries);

	if(reference == NULL){
		printf("Error reading reference file.\nExiting..\n");
		return -1;
	}
	if (query == NULL){
		printf("Error reading query file.\nExiting..\n");
		return -1;
	}

	//sortSequencesByLength(reference, ref_entries, 'd');		// 'd' = descending order
	//sortSequencesByLength(query, qry_entries, 'd');


	printf("Reference lines: %d\n", ref_entries);
	printf("Query lines: %d\n", qry_entries);
	
	
	// Alignment results: ref_entries rows, qry_entries columns
	int ** scores = (int **)malloc(sizeof(int*) * ref_entries);
	for (i = 0; i < ref_entries; i++){
		scores[i] = (int *)calloc(qry_entries, sizeof(int));
	}

	// Getting local memory size
	cl_int error = 0;
	cl_ulong device_local_memory = 0;
	fbclGetDeviceLocalMemory(device, &device_local_memory, &error);
	if (error) printf("** Error getting device local memory. **\n");
	printf("\n\n\n[# INFO: Device Local memory: %ld bytes #]\n\n", device_local_memory);


	printf("\n\n");


	switch (alg_type)
	{
	case LINEAR:
#ifdef WINDOWS
		source_code = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw_parallel.cl", &source_size);
#else
		source_code = fbclGetSourceFromFile("politosw_parallel.cl", &source_size);
#endif
		printf("Algorithm used for the alignment: SW-Linear \n");
		break;
	case AFFINE:
#ifdef WINDOWS
		source_code = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw_affine_parallel.cl", &source_size);
#else
		source_code = fbclGetSourceFromFile("politosw_affine_parallel.cl", &source_size);
#endif
		printf("Algorithm used for the alignment: SW-Affine \n");
		break;
	case DGS:
#ifdef WINDOWS
		source_code = fbclGetSourceFromFile("C:/Users/stefano/Dropbox/openCL/politosw/openCL_test1/politosw_dgssw_parallel.cl", &source_size);
#else
		source_code = fbclGetSourceFromFile("politosw_dgssw_parallel.cl", &source_size);
#endif
		printf("Algorithm used for the alignment: SW-DGS \n");
		break;
	default:
		break;
	}


	cl_context _context;
	cl_command_queue _queue;
	cl_program _program;
	cl_kernel _kernel;

	avg_time = 0.0;
	TOT_PROFILE_START();


	fbclInitKernel(device, &_context, &_queue, source_code, source_size, &_program, &_kernel, (char*)"politosw", YES);

	/************************************************* Parallel execution *****************************************************/

	
	GLOBAL_INFO info;
	info.penality1 = -1;
	info.penality2 = -10;
	info.symbol_number = strlen(symbols);
	info.block_size = worker;
	info.groups_per_line = LINEAR_DIMENSION / worker;
	
	//BLOCK_LINEAR_SHORT * kernel_array_linear;

	size_t size_char_query = 1, size_char_reference = 1;
	char * char_buffer_reference = NULL;
	char * char_buffer_query = NULL;// = (char*)malloc(sizeof(char));

	int max_reference_length = 0, max_query_length = 0, sub_size = 0;
	int local_memory_base_usage = 0;
	int local_memory_medium_usage = 0;
	int local_memory_max_usage = 0;
	// alignment specific

	
	// Setup alignment-specific data
	
	int tot_alignments = ref_entries * qry_entries;
	
	info.total_alignments = tot_alignments;
	info.parallel_jobs = parallel_jobs;

	
	//int jobs = LINEAR_DIMENSION / max_worker;				// horizontal groups
	int jobs = parallel_jobs;
	
	
	CONTROL_PARALLEL * ctrl = (CONTROL_PARALLEL *)calloc(tot_alignments, sizeof(CONTROL_PARALLEL));

	if(ctrl == NULL){
		perror("Malloc error: ");
		return -1;
	}

	CONTROL_PARALLEL * tmp_ctrl_ref = (CONTROL_PARALLEL *)calloc(ref_entries, sizeof(CONTROL_PARALLEL));
	CONTROL_PARALLEL * tmp_ctrl_qry = (CONTROL_PARALLEL *)calloc(qry_entries, sizeof(CONTROL_PARALLEL));

	for(i = 0; i < ref_entries; i++){
		unsigned long offset = 0;
		size_t tmp_size = 0;
		size_t final_dest_size = 0;

		char * tmp_buff = getCharBuffer(reference[i].seq, reference[i].seq_length, symbols, worker, &tmp_size);

		if (char_buffer_reference == NULL){
			//printf("Ref buffer is null: let's first fill it\n");
			char_buffer_reference = getCharBuffer(reference[i].seq, reference[i].seq_length, symbols, worker, &final_dest_size);
		}
		else{
			//printf("Realloc'ing ref global buffer\n");
			char_buffer_reference = reallocGlobalBuffer(char_buffer_reference, size_char_reference, tmp_buff, tmp_size, &final_dest_size, &offset);
		}
		
		size_char_reference = final_dest_size;

		if (tmp_size > max_reference_length) max_reference_length = tmp_size;

		/* //DEBUG: print char buffer

		printf("**DEBUG:\ni=%d:\n seq=%s\n index=%d\n", i, reference[i].seq, reference[i].index);
		printf("\nnew reference[%d]: length=%d, offset=%ld\nseq= ", i, ctrl[i].reference_padded_lenght, ctrl[i].ref_offset);
		for (j = 0; j < ctrl[i].reference_padded_lenght; j++) printf("%d", tmp_buff[j]);

		printf("\nchar_buffer_reference: length=%ld\nseq= ", size_char_reference);
		for (j = 0; j < size_char_reference; j++){
		if(char_buffer_reference[j] == 0)
		printf("0");
		else
		printf("%d", char_buffer_reference[j]);
		}
		printf("\nreference at offset %ld:\n", ctrl[i].ref_offset);
		for(j=0; j<ctrl[i].reference_padded_lenght; j++) printf("%d", char_buffer_reference[ctrl[i].ref_offset + j]);

		//printf("\n\n");
		*/

		free(tmp_buff);
		tmp_buff = NULL;


		tmp_ctrl_ref[i].reference_length = reference[i].seq_length;
		tmp_ctrl_ref[i].reference_padded_length = tmp_size;
		tmp_ctrl_ref[i].block_y = tmp_size / info.block_size;
		tmp_ctrl_ref[i].ref_id = reference[i].index;
		tmp_ctrl_ref[i].ref_offset = (cl_long)offset;
		//printf("ref %d offset = %ld\n", i, tmp_ctrl_ref[i].ref_offset);

		if ((offset % 16) != 0){
			printf("ERROR: offset/16=%lf\n", (double)offset / 16);
			return(-1);
		}

		// Compute internal parameters to efficiently copy from global to local memory
		short ref_complete_blocks = 0, ref_pad_block_copiers = 0, ref_single_pad_copiers = 0;
		divideEvenlyBlocks(tmp_ctrl_ref[i].reference_padded_length, (short)worker, (short)16, &ref_complete_blocks, &ref_pad_block_copiers, &ref_single_pad_copiers);
		tmp_ctrl_ref[i].ref_complete_blocks = ref_complete_blocks;
		tmp_ctrl_ref[i].ref_pad_block_copiers = ref_pad_block_copiers;
		tmp_ctrl_ref[i].ref_single_pad_copiers = ref_single_pad_copiers;

		short sub_complete_copiers = 0, sub_complete_elements = 0, sub_pad_elements = 0;
		sub_size = sizeof(char_sub_matrix) / sizeof(char_sub_matrix[0]);
		divideEvenly(sub_size, (short)worker, &sub_complete_copiers, &sub_complete_elements, &sub_pad_elements);
		tmp_ctrl_ref[i].sub_complete_copiers = sub_complete_copiers;
		tmp_ctrl_ref[i].sub_complete_elements = sub_complete_elements;
		tmp_ctrl_ref[i].sub_pad_elements = sub_pad_elements;


		
	}

	for(j = 0; j < qry_entries; j++){
		unsigned long offset = 0;
		size_t tmp_size = 0;
		size_t final_dest_size = 0;

		char * tmp_buff = getCharBuffer(query[j].seq, query[j].seq_length, symbols, worker, &tmp_size);

		if (char_buffer_query == NULL){
			char_buffer_query = getCharBuffer(query[j].seq, query[j].seq_length, symbols, worker, &final_dest_size);
		}
		else{
			//printf("Realloc'ing qry global buffer\n");
			char_buffer_query = reallocGlobalBuffer(char_buffer_query, size_char_query, tmp_buff, tmp_size, &final_dest_size, &offset);
		}
		size_char_query = final_dest_size;

		if (tmp_size > max_query_length) max_query_length = tmp_size;


		free(tmp_buff);
		tmp_buff = NULL;

		tmp_ctrl_qry[j].query_length = query[j].seq_length;
		tmp_ctrl_qry[j].query_padded_length = tmp_size;
		tmp_ctrl_qry[j].block_x = tmp_size / info.block_size;
		tmp_ctrl_qry[j].qry_id = query[j].index;
		tmp_ctrl_qry[j].qry_offset = (long)offset;
		//printf("qry %d offset = %ld\n", i, tmp_ctrl_qry[j].qry_offset);
		if ((offset % 16) != 0){
			printf("ERROR: offset/16=%lf\n", (double)offset / 16);
			return(-1);
		}

		// Compute internal parameters to efficiently copy from global to local memory

		short qry_complete_blocks = 0, qry_pad_block_copiers = 0, qry_single_pad_copiers = 0;
		divideEvenlyBlocks(tmp_ctrl_qry[j].query_padded_length, (short)worker, (short)16, &qry_complete_blocks, &qry_pad_block_copiers, &qry_single_pad_copiers);
		tmp_ctrl_qry[j].qry_complete_blocks = qry_complete_blocks;
		tmp_ctrl_qry[j].qry_pad_block_copiers = qry_pad_block_copiers;
		tmp_ctrl_qry[j].qry_single_pad_copiers = qry_single_pad_copiers;
	}

	
	// Fill the MAIN CONTROL structure
	int al = 0;
	
	for(i = 0; i < ref_entries; i++){
		for(j = 0; j < qry_entries; j++){
			
			ctrl[al].reference_length = tmp_ctrl_ref[i].reference_length;
			ctrl[al].reference_padded_length = tmp_ctrl_ref[i].reference_padded_length;
			ctrl[al].block_y = tmp_ctrl_ref[i].block_y;
			

			ctrl[al].ref_id = tmp_ctrl_ref[i].ref_id;
			ctrl[al].ref_offset = tmp_ctrl_ref[i].ref_offset;
			ctrl[al].ref_complete_blocks = tmp_ctrl_ref[i].ref_complete_blocks;
			ctrl[al].ref_pad_block_copiers = tmp_ctrl_ref[i].ref_pad_block_copiers;
			ctrl[al].ref_single_pad_copiers = tmp_ctrl_ref[i].ref_single_pad_copiers;
			
			ctrl[al].sub_complete_copiers = tmp_ctrl_ref[i].sub_complete_copiers;
			ctrl[al].sub_complete_elements = tmp_ctrl_ref[i].sub_complete_elements;
			ctrl[al].sub_pad_elements = tmp_ctrl_ref[i].sub_pad_elements;
			

			ctrl[al].query_length = tmp_ctrl_qry[j].query_length;
			ctrl[al].query_padded_length = tmp_ctrl_qry[j].query_padded_length;
			
			
			ctrl[al].qry_id = tmp_ctrl_qry[j].qry_id;
			ctrl[al].qry_offset = tmp_ctrl_qry[j].qry_offset;
			ctrl[al].qry_complete_blocks = tmp_ctrl_qry[j].qry_complete_blocks;
			ctrl[al].qry_pad_block_copiers = tmp_ctrl_qry[j].qry_pad_block_copiers;
			ctrl[al].qry_single_pad_copiers = tmp_ctrl_qry[j].qry_single_pad_copiers;
			ctrl[al].block_x = tmp_ctrl_qry[j].block_x;


			

			
			/*
			printf("\n\n");
			printf("ref index: %d, qry index: %d, workers=%d\n", reference[i].index, query[0].index, worker);
			printf("DEBUG: ref_l = %d, ref_complete_blocks = %d, ref_pad_block_copiers = %d, ref_single_pad_copiers = %d, ref_offset = %ld\n", ctrl[al].reference_padded_lenght, ctrl[al].ref_complete_blocks, ctrl[al].ref_pad_block_copiers, ctrl[al].ref_single_pad_copiers, ctrl[al].ref_offset);
			printf("DEBUG: qry_l = %d, qry_complete_blocks = %d, qry_pad_block_copiers = %d, qry_single_pad_copiers = %d, qry_offset = %ld\n", ctrl[al].query_padded_length, ctrl[al].qry_complete_blocks, ctrl[al].qry_pad_block_copiers, ctrl[al].qry_single_pad_copiers, ctrl[al].qry_offset);
			printf("DEBUG: sub_l = %d, sub_complete_copiers = %d, sub_complete_elements = %d, sub_pad_elements = %d\n\n", sub_size, ctrl[al].sub_complete_copiers, ctrl[al].sub_complete_elements, ctrl[al].sub_pad_elements);
			*/

			local_memory_base_usage = reference[i].seq_length + query[j].seq_length + sub_size;
			local_memory_medium_usage = local_memory_base_usage + sizeof(BLOCK_LINEAR_SHORT)*info.block_size;
			local_memory_max_usage = local_memory_medium_usage + sizeof(BLOCK_AUX_LINEAR_SHORT)*max_reference_length;

			//printf("Local memory usage: base = %d, medium = %d, max = %d\n\n\n\n", local_memory_base_usage, local_memory_medium_usage, local_memory_max_usage);
			
			al++;
		}
	}
	//printf("\n\nDEBUG: max_ref_length = %d\n\n", max_reference_length);
	
	
	//buffer_aux_row_linear = get_buffer_aux_linear(max_reference_length);
	//kernel_array_linear = get_block_linear(worker*jobs);

	printf("\n\nRunning kernel..\n");
	
	/*
	printf("Debug informations:\n");
	for(i=0; i<tot_alignments; i++){
		printf("\n\n\n**i=%d:\nctrl[%d].reference_length: %d\n", i, i, ctrl[i].reference_length);
		printf("*ctrl[%d].block_y: %d\n", i, ctrl[i].block_y);
		
		printf("*ctrl[%d].query_length: %d\n", i, ctrl[i].query_padded_length);
		printf("*ctrl[%d].block_x: %d\n", i, ctrl[i].block_x);
		
		printf("*ctrl[%d].ref_id: %ld\n", i, ctrl[i].ref_id);
		printf("*ctrl[%d].ref_offset: %ld\n", i, ctrl[i].ref_offset);
		printf("*ctrl[%d].ref_complete_blocks: %d\n", i, ctrl[i].ref_complete_blocks);
		printf("*ctrl[%d].ref_pad_block_copiers: %d\n", i, ctrl[i].ref_pad_block_copiers);
		printf("*ctrl[%d].ref_single_pad_copiers: %d\n", i, ctrl[i].ref_single_pad_copiers);
		printf("*ctrl[%d].ref_pad_block_copiers: %d\n", i, ctrl[i].ref_pad_block_copiers);
		
		printf("*ctrl[%d].qry_id: %ld\n", i, ctrl[i].qry_id);
		printf("*ctrl[%d].qry_offset: %ld\n", i, ctrl[i].qry_offset);
		printf("*ctrl[%d].qry_complete_blocks: %d\n", i, ctrl[i].qry_complete_blocks);
		printf("*ctrl[%d].qry_pad_block_copiers: %d\n", i, ctrl[i].qry_pad_block_copiers);
		printf("*ctrl[%d].qry_single_pad_copiers: %d\n", i, ctrl[i].qry_single_pad_copiers);
		printf("*ctrl[%d].qry_pad_block_copiers: %d\n", i, ctrl[i].qry_pad_block_copiers);
		
		printf("*ctrl[%d].sub_complete_copiers: %d\n", i, ctrl[i].sub_complete_copiers);
		printf("*ctrl[%d].sub_complete_elements: %d\n", i, ctrl[i].sub_complete_elements);
		printf("*ctrl[%d].sub_pad_elements: %d\n", i, ctrl[i].sub_pad_elements);
		
		printf("*ctrl[%d].score: %d\n", i, ctrl[i].score);

		//getchar();
	}*/
	
	

	err = 0;
	execution_time = 0.0;
	size_t score_matrix_size = sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2);
	size_t size_block = 0, size_block_aux = 0;

	printf("******************************** DEBUG ********************************\n");
	printf("Dimensions::\n  size_char_reference dimension=%ld\n  size_char_query dimension=%ld\n", sizeof(char) * size_char_reference, sizeof(char) * size_char_query);
	printf("  char_sub_matrix=%ld\n  info=%ld\n  control=%ld\n\n", score_matrix_size, sizeof(GLOBAL_INFO), sizeof(CONTROL_PARALLEL)*tot_alignments);
	printf("Sum of global sizes=%ld\n\n", sizeof(char)*size_char_reference + sizeof(char)*size_char_query + score_matrix_size + sizeof(GLOBAL_INFO) + sizeof(CONTROL_PARALLEL)*tot_alignments);
	
	printf("Local sizes::\n  max_reference_length=%ld\n  max_query_length=%ld\n  sub_matrix=%ld\n", sizeof(char) * max_reference_length, sizeof(char) * max_query_length, score_matrix_size);
	switch(alg_type){
		case LINEAR:
			size_block = sizeof(BLOCK_LINEAR_SHORT);
			size_block_aux = sizeof(BLOCK_AUX_LINEAR_SHORT);
			break;
		case AFFINE:
			size_block = sizeof(BLOCK_AFFINE);
			size_block_aux = sizeof(BLOCK_AUX_AFFINE);
			break;
		case DGS:
			size_block = sizeof(BLOCK_DGS);
			size_block_aux = sizeof(BLOCK_AUX_DGS);
			break;
		default: break;
	}
	printf("  kernel_aux=%ld\n  kernel_matrix=%ld\n\n", size_block_aux * max_reference_length, size_block * info.block_size);
	printf("Sum of local sizes=%ld\n\n\n", sizeof(char)*max_reference_length + sizeof(char)*max_query_length + score_matrix_size + size_block_aux * max_reference_length + size_block*info.block_size);
	


	printf(" * Total alignments = %d\n\n", info.total_alignments);
	
	//getchar();
	
	//BLOCK_AUX_AFFINE * buffer_aux_affine = get_buffer_aux_affine(max_reference_length);
	//BLOCK_AFFINE * buffer_affine = get_buffer_affine(info.block_size);

	printf("****************************** RUN KERNEL *****************************\n");
	PROFILE_START();
	switch (alg_type){
	case LINEAR:
		err = fbclRunParallel(device, &_context, &_queue, &_kernel, parallel_jobs, info.block_size, (char*)"politosw", YES, &execution_time,
			(char*)"rrrralllll",
			
			// Global parameters
			sizeof(char) * size_char_reference, (void*)char_buffer_reference,
			sizeof(char) * size_char_query, (void*)char_buffer_query,
			sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)char_sub_matrix,
			sizeof(GLOBAL_INFO), (void*)&info,
				
			sizeof(CONTROL_PARALLEL)*tot_alignments, (void*)ctrl,
			
			// Local parameters
			sizeof(char) * max_reference_length, NULL,
			sizeof(char) * max_query_length, NULL,
			sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), NULL,
			size_block_aux * max_reference_length, NULL,
			size_block * info.block_size, NULL

			);
		
		break;
	case AFFINE:
		err = fbclRunParallel(device, &_context, &_queue, &_kernel, parallel_jobs, info.block_size, (char*)"politosw", YES, &execution_time,
			(char*)"rrrralllll",

			// Global parameters
			sizeof(char) * size_char_reference, (void*)char_buffer_reference,
			sizeof(char) * size_char_query, (void*)char_buffer_query,
			sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)char_sub_matrix,
			sizeof(GLOBAL_INFO), (void*)&info,

			sizeof(CONTROL_PARALLEL)*tot_alignments, (void*)ctrl,
			//sizeof(BLOCK_AUX_AFFINE)*max_reference_length, (void*)buffer_aux_affine,
			//sizeof(BLOCK_AFFINE)*info.block_size, (void*)buffer_affine,

			// Local parameters
			sizeof(char) * max_reference_length, NULL,
			sizeof(char) * max_query_length, NULL,
			sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), NULL,
			size_block_aux * max_reference_length, NULL,									// TOO BIG to keep in local; GLOBAL
			size_block * info.block_size, NULL

			);
			break;
	case DGS:
		err = fbclRunParallel(device, &_context, &_queue, &_kernel, parallel_jobs, info.block_size, (char*)"politosw", YES, &execution_time,
			(char*)"rrrralllll",

			// Global parameters
			sizeof(char) * size_char_reference, (void*)char_buffer_reference,
			sizeof(char) * size_char_query, (void*)char_buffer_query,
			sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), (void*)char_sub_matrix,
			sizeof(GLOBAL_INFO), (void*)&info,

			sizeof(CONTROL_PARALLEL)*tot_alignments, (void*)ctrl,
			//sizeof(BLOCK_AUX_AFFINE)*max_reference_length, (void*)buffer_aux_affine,
			//sizeof(BLOCK_AFFINE)*info.block_size, (void*)buffer_affine,

			// Local parameters
			sizeof(char) * max_reference_length, NULL,
			sizeof(char) * max_query_length, NULL,
			sizeof(char) * (info.symbol_number + 2) * (info.symbol_number + 2), NULL,
			size_block_aux * max_reference_length, NULL,									
			size_block * info.block_size, NULL

			);
		break;
	default:
		printf("** Error: Invalid algorithm type (%d).\n", alg_type);
		break;
	}
	PROFILE_STOP();
	
	printf("\nKernel finished! err = %d\n", err);
	
	
	if(err){
		printf("Error with fbclRunParallel..\nFreeing memory and exiting\n");
		
		free(char_buffer_reference);
		free(char_buffer_query); 
		fbclDestroy(&_context, &_program, &_kernel, &_queue);

		
		free(source_code);
		printf("source\n");
		free(reference);
		printf("ref\n");
		free(query);
		printf("qry\n");
		free(scores);
		printf("scores\n");
		return(-1);
	}


	al = 0;
	
	TOT_PROFILE_STOP();


	/************************************************* Execution ended: destroy OpenCL structures and save output *****************************************************/
	
	fbclDestroy(&_context, &_program, &_kernel, &_queue);

	
	
	FILE * fp;
	if ((fp = fopen(out_file, "w")) == NULL){
		printf("Error opening ouput file. Not saving results\n");
	}

	/*
	for (i = 0; i<tot_alignments; i++){
		printf("\n\n\n**i=%d:\nctrl[%d].reference_length: %d\n", i, i, ctrl[i].reference_length);
		printf("*ctrl[%d].block_y: %d\n", i, ctrl[i].block_y);

		printf("*ctrl[%d].query_length: %d\n", i, ctrl[i].query_padded_length);
		printf("*ctrl[%d].block_x: %d\n", i, ctrl[i].block_x);

		printf("*ctrl[%d].ref_id: %ld\n", i, ctrl[i].ref_id);
		printf("*ctrl[%d].ref_offset: %ld\n", i, ctrl[i].ref_offset);
		printf("*ctrl[%d].ref_complete_blocks: %d\n", i, ctrl[i].ref_complete_blocks);
		printf("*ctrl[%d].ref_pad_block_copiers: %d\n", i, ctrl[i].ref_pad_block_copiers);
		printf("*ctrl[%d].ref_single_pad_copiers: %d\n", i, ctrl[i].ref_single_pad_copiers);
		printf("*ctrl[%d].ref_pad_block_copiers: %d\n", i, ctrl[i].ref_pad_block_copiers);

		printf("*ctrl[%d].qry_id: %ld\n", i, ctrl[i].qry_id);
		printf("*ctrl[%d].qry_offset: %ld\n", i, ctrl[i].qry_offset);
		printf("*ctrl[%d].qry_complete_blocks: %d\n", i, ctrl[i].qry_complete_blocks);
		printf("*ctrl[%d].qry_pad_block_copiers: %d\n", i, ctrl[i].qry_pad_block_copiers);
		printf("*ctrl[%d].qry_single_pad_copiers: %d\n", i, ctrl[i].qry_single_pad_copiers);
		printf("*ctrl[%d].qry_pad_block_copiers: %d\n", i, ctrl[i].qry_pad_block_copiers);

		printf("*ctrl[%d].sub_complete_copiers: %d\n", i, ctrl[i].sub_complete_copiers);
		printf("*ctrl[%d].sub_complete_elements: %d\n", i, ctrl[i].sub_complete_elements);
		printf("*ctrl[%d].sub_pad_elements: %d\n", i, ctrl[i].sub_pad_elements);

		printf("*ctrl[%d].score: %d\n", i, ctrl[i].score);

		//getchar();
	}*/

	printf("\n\nFinal scores:\n");
	


	char buff[MAX_LINE];
	char type[MAX_LINE];
	switch(alg_type){
		case LINEAR:
			sprintf(type, "LINEAR");
			break;
		case AFFINE:
			sprintf(type, "AFFINE");
			break;
		case DGS:
			sprintf(type, "DGS");
			break;
		default:
			sprintf(type, "(unknown)");
			break;
	}
	sprintf(buff, "PolitoSW [type: %s] run on %s and %s.\n\n", type, subject_db, query_file);
	printf("%s", buff);
	
	fprintf(fp, "%s", buff);
	
	for(i=0; i<tot_alignments; i++){
		
		#ifdef WINDOWS
		sprintf_s(buff, "[r:%3ld (l=%4d) ** q:%3ld (l=%4d)] = %6d - %10.5f ms  --  [ idx = %3d: group (%2d, %2d), iter %2d ]\n", (unsigned long)ctrl[i].ref_id, (int)ctrl[i].reference_length, (unsigned long)ctrl[i].qry_id, (int)ctrl[i].query_length, (int)ctrl[i].score, execution_time, (int)ctrl[i].idx, (int)ctrl[i].gid_x, (int)ctrl[i].gid_y, (int)ctrl[i].iter);
		#else
		sprintf(buff, "[r:%3ld (l=%4d) ** q:%3ld (l=%4d)] = %6d - %10.5f ms  --  [ idx = %3d: group (%2d, %2d), iter %2d ]\n", (unsigned long)ctrl[i].ref_id, (int)ctrl[i].reference_length, (unsigned long)ctrl[i].qry_id, (int)ctrl[i].query_length, (int)ctrl[i].score, execution_time, (int)ctrl[i].idx, (int)ctrl[i].gid_x, (int)ctrl[i].gid_y, (int)ctrl[i].iter);
		#endif
		printf("%s", buff);
		
		fprintf(fp, "%s\n\n", buff);

	}
	printf("\n\nTotal execution time = %lf ms\n", tot_exec_time);
	printf("\n");
	fclose(fp);	
	
	
	

	//free(buffer_affine);
	//free(buffer_aux_affine);

	free(tmp_ctrl_ref);
	free(tmp_ctrl_qry);
		
	free(ctrl);
	
	free(char_buffer_reference);
	free(char_buffer_query);
	
	free(source_code);
	for(i=0; i<ref_entries; i++) free(reference[i].seq);
	free(reference);
	for(i=0; i<qry_entries; i++) free(query[i].seq);
	free(query);
	
	for (i = 0; i < ref_entries; i++){
		free(scores[i]);
	}
	free(scores);

	return 0;
}


