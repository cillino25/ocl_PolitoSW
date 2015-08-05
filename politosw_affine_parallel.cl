
typedef struct _block_aux_affine {
	short h;
	short e;
	short f;
	short e_diagonal;
	short f_diagonal;
	short e_aux;
	short f_aux;
	/*cl_int h;
	cl_int e;
	cl_int f;
	cl_int e_diagonal;
	cl_int f_diagonal;
	cl_int e_aux;
	cl_int f_aux;*/
} BLOCK_AUX_AFFINE;

typedef struct _block_affine {
	int up;
	int left;
	int diagonal;
	int h_max;
	int e;
	int f;
	int e_diagonal;
	int f_diagonal;
	int e_aux;
	int f_aux;
} BLOCK_AFFINE;

typedef struct _global_info {
	int symbol_number;
	int block_size;
	int penality1;
	int penality2;
	
	int groups_per_line;
	int total_alignments;
	int parallel_jobs;
} GLOBAL_INFO;

typedef struct _control_parallel {
	int	reference_length;
	int reference_padded_length;
	int	query_length;
	int query_padded_length;
	int  block_x;
	int  block_y;

	ulong ref_id;
	ulong ref_offset;
	int   ref_complete_blocks;
	int   ref_pad_block_copiers;
	int   ref_single_pad_copiers;

	ulong qry_id;
	ulong qry_offset;
	int   qry_complete_blocks;
	int   qry_pad_block_copiers;
	int   qry_single_pad_copiers;

	short sub_complete_copiers;
	short sub_complete_elements;
	short sub_pad_elements;
	
	int gid_x;
	int gid_y;
	int iter;
	int idx;
	
	int score;
} CONTROL_PARALLEL;


void charBlockCopy16(__constant char * source, long offset, __local char * dest, ushort local_size, ushort k, short complete_blocks, short pad_block_copiers, short single_pad_copiers){
	int i = 0;
	long off = offset / 16;		// Allowed since host code gives 16-wide aligned sequences
	
	// Block copy (each worker can copy a 16 wide value)
	for(i=0; i < complete_blocks; i++){
		long base = k*complete_blocks + i;
		
		char16 t = vload16(base + off, source);
		
		dest[base * 16 +  0] = t.s0;
		dest[base * 16 +  1] = t.s1;
		dest[base * 16 +  2] = t.s2;
		dest[base * 16 +  3] = t.s3;
		dest[base * 16 +  4] = t.s4;
		dest[base * 16 +  5] = t.s5;
		dest[base * 16 +  6] = t.s6;
		dest[base * 16 +  7] = t.s7;
		dest[base * 16 +  8] = t.s8;
		dest[base * 16 +  9] = t.s9;
		dest[base * 16 + 10] = t.sA;
		dest[base * 16 + 11] = t.sB;
		dest[base * 16 + 12] = t.sC;
		dest[base * 16 + 13] = t.sD;
		dest[base * 16 + 14] = t.sE;
		dest[base * 16 + 15] = t.sF;

	}
	
	
	// Block padding
	// Only block_copiers can copy 16 wide values
	if(k < pad_block_copiers){	
		long base = local_size * complete_blocks + k;
		
		char16 t = vload16(base + off, source);
		
		dest[base * 16 +  0] = t.s0;
		dest[base * 16 +  1] = t.s1;
		dest[base * 16 +  2] = t.s2;
		dest[base * 16 +  3] = t.s3;
		dest[base * 16 +  4] = t.s4;
		dest[base * 16 +  5] = t.s5;
		dest[base * 16 +  6] = t.s6;
		dest[base * 16 +  7] = t.s7;
		dest[base * 16 +  8] = t.s8;
		dest[base * 16 +  9] = t.s9;
		dest[base * 16 + 10] = t.sA;
		dest[base * 16 + 11] = t.sB;
		dest[base * 16 + 12] = t.sC;
		dest[base * 16 + 13] = t.sD;
		dest[base * 16 + 14] = t.sE;
		dest[base * 16 + 15] = t.sF;
	}
	
	// Padding // useless if multiple of 16
	// 
	if((k - pad_block_copiers < single_pad_copiers)&&(k - pad_block_copiers > 0)){
		long loc_base = local_size * complete_blocks * 16 + pad_block_copiers * 16;
		dest[loc_base + k - pad_block_copiers] = source[off + loc_base + k - pad_block_copiers];
	}
}

__kernel void politosw(
	__constant char *gl_reference, 
	__constant char *gl_query, 
	__constant char *gl_sub_matrix, 
	__global GLOBAL_INFO *info,
	
	__global CONTROL_PARALLEL *ctrl,
	//__global BLOCK_AUX_AFFINE *aux_row,  
	//__global BLOCK_AFFINE *kernel_matrix,

	__local char * reference,
	__local char * query,
	__local char * sub_matrix,
	__local BLOCK_AUX_AFFINE *aux_row,  
	__local BLOCK_AFFINE *kernel_matrix
	)
	{

	// lascio i valori interni a INT solo per compatibilità con schede che non supportano il char
	int valid, finish;

	int loc_size = get_local_size(0);
	int gid_x = get_group_id(0);
	int gid_y = get_group_id(1);
	//ushort num_groups_x = info->groups_per_line;
	int num_groups_x = get_num_groups(0);
	int num_groups_y = get_num_groups(1);
	
	int y = get_local_id(1);
	
	int k = get_local_id(0);
	int k1=(k+1) % loc_size;


	int total_alignments = info->total_alignments;
	int jobs = info->parallel_jobs;
	
	int al_index = 0;
	int al_iter = 0;

	while(1){
		loc_size = get_local_size(0);
		
		al_index = al_iter*jobs + gid_y*num_groups_x + gid_x;
		
		if(al_index >= total_alignments) break;


		int x, y, i, ib, bx, by, x1, y1;
		i=0;
		ib=0-k;
		x=0-k;
		y=k;	
		bx=by=0;

		// Copy reference into local memory
		int a = ctrl[al_index].ref_complete_blocks;
		int b = ctrl[al_index].ref_pad_block_copiers;
		int c = ctrl[al_index].ref_single_pad_copiers;
		ulong offset = ctrl[al_index].ref_offset;
		charBlockCopy16(gl_reference, offset, reference, loc_size, k, a, b, c);


		// Copy query into local memory
		a = ctrl[al_index].qry_complete_blocks;
		b = ctrl[al_index].qry_pad_block_copiers;
		c = ctrl[al_index].qry_single_pad_copiers;
		offset = ctrl[al_index].qry_offset;

		charBlockCopy16(gl_query, offset, query, loc_size, k, a, b, c);


		// Copy sub_matrix into local memory
		a = ctrl[al_index].sub_complete_copiers;
		b = ctrl[al_index].sub_complete_elements;
		c = ctrl[al_index].sub_pad_elements;
	
		if(k < a){
			for(i=0; i < b; i++){
				sub_matrix[k*b + i] = gl_sub_matrix[k*b + i];
			}
		}

		if(k == (loc_size - 1)){
			short base = b * a;
			for(i=0; i < c; i++){
				sub_matrix[base + i] = gl_sub_matrix[base + i];
			}
		}

		// Copy INFO structure in private memory
		int loc_block_x, loc_block_y, loc_block_size;
		int d1, d2;
		int n;
		n = info->symbol_number;
		d1 = info->penality1;
		d2 = info->penality2;
		loc_block_x = ctrl[al_index].block_x;
		loc_block_y = ctrl[al_index].block_y;
		loc_block_size = loc_size;
	
		int h, h_max, r, q, s;
		int up, left, left1, diagonal;
		int e, f, e_aux, f_aux;

		int state=0;

		finish=0;
		h_max=0;
		
		
		barrier(CLK_LOCAL_MEM_FENCE);
		//barrier(CLK_GLOBAL_MEM_FENCE);
		
		
		while(!finish)
		{	
			valid=1;
			if (ib<0 || state==3)
				valid=0;

			left1 = kernel_matrix[k].left;
			//Per l'affine controllare dove salvare e e f, anche al cambio del blocco/riga
		
			h=0;
			e=0;
			f=0;
			e_aux = kernel_matrix[k].e_aux;
			f_aux = kernel_matrix[k].f_aux;
			
			if(valid) {
				r= reference[y];
				q= query[x];
				if (q>0 && q<=n && r>0 && r<=n){				

					s=sub_matrix[r*(n+2)+q];		    
				
					left	=	kernel_matrix[k].left;
					up		=	kernel_matrix[k].up;
					diagonal=	kernel_matrix[k].diagonal;

					e = kernel_matrix[k].e + d2;
					f = kernel_matrix[k].f + d2;
				
					e = max(left+d1, e);
					e = max(0, e);
					f = max(up+d1, f);
					f = max(0, f);

					h=max(h,diagonal+s);
					h=max(h,kernel_matrix[k].e_diagonal+s);
					h=max(h,kernel_matrix[k].f_diagonal+s);


					h_max = max(h_max, h);			
				}	
			
			}
		
		
			// Aggiornamento degli indici		
			if (state!=3){
				ib+=1;
				x1=x;
				y1=y;

				if (ib==loc_block_size){
					//CAMBIO BLOCCO
					state=1;
					by+=1;
					ib=0;
					if (by==loc_block_y){				
						//CAMBIO RIGA
						state=2;
						by=0;
						bx+=1;
						if (bx==loc_block_x){
							//FINE
							state=3;					
						}
					}
					x = bx*loc_block_size;
					y = by*loc_block_size + k;					
				}
				else{
					state = 0;
					x++;
				}
			}
		
								
			//barrier(CLK_LOCAL_MEM_FENCE);	
			barrier(CLK_GLOBAL_MEM_FENCE);	
		
			if(valid){
				if (state==0){
					// Verso altri
					if (y1==(loc_block_y * loc_block_size - 1));
					else{
						kernel_matrix[k1].left=h;
						kernel_matrix[k1].e=e;

						kernel_matrix[k1].f_aux = f;
						kernel_matrix[k1].e_aux = e;
					}
				
					//Verso se stessi
					kernel_matrix[k].up=h;		
					kernel_matrix[k].diagonal=left1;
					kernel_matrix[k].f=f;

					kernel_matrix[k].e_diagonal = e_aux;
					kernel_matrix[k].f_diagonal = f_aux;

					//Salvataggio				
				}
				else if (state==1 ){
					// Verso altri
					kernel_matrix[k1].left=h;
					kernel_matrix[k1].e=e;
					kernel_matrix[k1].f_aux = f;
					kernel_matrix[k1].e_aux = e;
				
					//Verso se stessi
					kernel_matrix[k].up=aux_row[y].h;		
					kernel_matrix[k].diagonal=aux_row[y-1].h;
					kernel_matrix[k].f=aux_row[y].f;

					kernel_matrix[k].e_diagonal = aux_row[y].e_aux;
					kernel_matrix[k].f_diagonal = aux_row[y].f_aux;

					//Salvataggio
					aux_row[y1].h = h;
					aux_row[y1].f = f;

					aux_row[y1].e_aux = e_aux;
					aux_row[y1].f_aux = f_aux;

				}
				else if (state==2 ){
					// Verso altri
					if (y1 == (loc_block_y * loc_block_size - 1));
					else{
						kernel_matrix[k1].left=h;
						kernel_matrix[k1].e=e;
						kernel_matrix[k1].f_aux = f;
						kernel_matrix[k1].e_aux = e;
					}

					// Verso se stessi
					kernel_matrix[k].up=aux_row[y].h;		
					kernel_matrix[k].f=aux_row[y].f;
				

					if (y==0){
						kernel_matrix[k].diagonal = 0;
						kernel_matrix[k].left = 0;
						kernel_matrix[k].e = 0;
						kernel_matrix[k].e_diagonal = 0;
						kernel_matrix[k].f_diagonal = 0;
					}
					else{
						kernel_matrix[k].diagonal=aux_row[y-1].h;
						kernel_matrix[k].e_diagonal = aux_row[y].e_aux;
						kernel_matrix[k].f_diagonal = aux_row[y].f_aux;
					}

					//Salvataggio
					aux_row[y1].h = h;
					aux_row[y1].f = f;		
					aux_row[y1].e_aux = e_aux;
					aux_row[y1].f_aux = f_aux;
				}
				else if(valid && state==3){
					kernel_matrix[k1].left=h;
					kernel_matrix[k1].e=e; 
				}
			}			
		
		
			barrier(CLK_LOCAL_MEM_FENCE);	
			barrier(CLK_GLOBAL_MEM_FENCE);
		
			i+=1;
			if(i == (loc_block_x * loc_block_y * loc_block_size + loc_block_size - 1))
				finish=1;
		}
	
		kernel_matrix[k].h_max = h_max;

		barrier(CLK_LOCAL_MEM_FENCE);
		barrier(CLK_GLOBAL_MEM_FENCE);

		// In-loco search of maximum alignment score
		while(loc_size > 2){
			loc_size = loc_size / 2;
			barrier(CLK_LOCAL_MEM_FENCE);
			if(k < loc_size){
				kernel_matrix[k].h_max = max(kernel_matrix[k].h_max, kernel_matrix[k + loc_size - 1].h_max);	// should be replaced with SELECT, more efficient
			}
		}
	
	
		// Copy back data from local to global memory
		if(k==0){
			atomic_xchg(&ctrl[al_index].score, kernel_matrix[0].h_max);
			atomic_xchg(&ctrl[al_index].gid_x, gid_x);
			atomic_xchg(&ctrl[al_index].idx, al_index);	
			atomic_xchg(&ctrl[al_index].gid_y, gid_y);
			atomic_xchg(&ctrl[al_index].iter, al_iter);
			
		}
		
		
		// end alignment operations
		
		al_iter++;
		barrier(CLK_LOCAL_MEM_FENCE);


	}

}
