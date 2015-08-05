
typedef struct _block_aux_linear {
	short h;
} BLOCK_AUX_LINEAR_SHORT;

typedef struct _block_linear {
	short up;
	short left;
	short diagonal;
	short h_max;
} BLOCK_LINEAR_SHORT;

typedef struct _info {
	int	reference_lenght;
	int	query_lenght;
	int symbol_number; 
	int block_size;
	int block_x;
	int block_y;	
	int penality1;
	int penality2;
} INFO;


typedef struct _control {
	short ref_complete_blocks;
	short ref_pad_block_copiers;
	short ref_single_pad_copiers;
	short qry_complete_blocks;
	short qry_pad_block_copiers;
	short qry_single_pad_copiers;
	short sub_complete_copiers;
	short sub_complete_elements;
	short sub_pad_elements;
} CONTROL_SHORT;


void charBlockCopy16(__constant char * source, __local char * dest, ushort local_size, ushort k, short complete_blocks, short pad_block_copiers, short single_pad_copiers){
	int i = 0;
	// Block copy (each worker can copy a 16 wide value)
	for(i=0; i < complete_blocks; i++){
		char16 t = vload16(k * complete_blocks + i, source);
		
		dest[k * complete_blocks * 16 + i * 16 +  0] = t.s0;
		dest[k * complete_blocks * 16 + i * 16 +  1] = t.s1;
		dest[k * complete_blocks * 16 + i * 16 +  2] = t.s2;
		dest[k * complete_blocks * 16 + i * 16 +  3] = t.s3;
		dest[k * complete_blocks * 16 + i * 16 +  4] = t.s4;
		dest[k * complete_blocks * 16 + i * 16 +  5] = t.s5;
		dest[k * complete_blocks * 16 + i * 16 +  6] = t.s6;
		dest[k * complete_blocks * 16 + i * 16 +  7] = t.s7;
		dest[k * complete_blocks * 16 + i * 16 +  8] = t.s8;
		dest[k * complete_blocks * 16 + i * 16 +  9] = t.s9;
		dest[k * complete_blocks * 16 + i * 16 + 10] = t.sA;
		dest[k * complete_blocks * 16 + i * 16 + 11] = t.sB;
		dest[k * complete_blocks * 16 + i * 16 + 12] = t.sC;
		dest[k * complete_blocks * 16 + i * 16 + 13] = t.sD;
		dest[k * complete_blocks * 16 + i * 16 + 14] = t.sE;
		dest[k * complete_blocks * 16 + i * 16 + 15] = t.sF;

	}
	
	
	// Block padding
	// Only block_copiers can copy 16 wide values
	if(k < pad_block_copiers){	
		
		char16 t = vload16(local_size * complete_blocks + k, source);
		dest[local_size * complete_blocks * 16 + k * 16 +  0] = t.s0;
		dest[local_size * complete_blocks * 16 + k * 16 +  1] = t.s1;
		dest[local_size * complete_blocks * 16 + k * 16 +  2] = t.s2;
		dest[local_size * complete_blocks * 16 + k * 16 +  3] = t.s3;
		dest[local_size * complete_blocks * 16 + k * 16 +  4] = t.s4;
		dest[local_size * complete_blocks * 16 + k * 16 +  5] = t.s5;
		dest[local_size * complete_blocks * 16 + k * 16 +  6] = t.s6;
		dest[local_size * complete_blocks * 16 + k * 16 +  7] = t.s7;
		dest[local_size * complete_blocks * 16 + k * 16 +  8] = t.s8;
		dest[local_size * complete_blocks * 16 + k * 16 +  9] = t.s9;
		dest[local_size * complete_blocks * 16 + k * 16 + 10] = t.sA;
		dest[local_size * complete_blocks * 16 + k * 16 + 11] = t.sB;
		dest[local_size * complete_blocks * 16 + k * 16 + 12] = t.sC;
		dest[local_size * complete_blocks * 16 + k * 16 + 13] = t.sD;
		dest[local_size * complete_blocks * 16 + k * 16 + 14] = t.sE;
		dest[local_size * complete_blocks * 16 + k * 16 + 15] = t.sF;
	}
	
	// Padding // useless if multiple of 16
	// 
	/*if(k - pad_block_copiers < single_pad_copiers){
		short base = local_size * complete_blocks * 16 + pad_block_copiers * 16;
		dest[base + k - pad_block_copiers] = source[base + k - pad_block_copiers];
	}*/
	if((k - pad_block_copiers < single_pad_copiers)&&(k - pad_block_copiers > 0)){
		short base = local_size * complete_blocks * 16 + pad_block_copiers * 16;
		dest[base + k - pad_block_copiers] = source[base + k - pad_block_copiers];
	}
}



// Optimized kernel parameters
__kernel void politosw(
	__constant char *gl_reference, 
	__constant char *gl_query, 
	__constant char * gl_sub_matrix, 
	__constant INFO *info,
	__constant CONTROL_SHORT *ctrl,
	__global int *al_score,
	//__global float * out_buffer,
	
	__local char * reference,
	__local char * query,
	__local char * sub_matrix, 
	__local BLOCK_AUX_LINEAR_SHORT *aux_row,  
	__local BLOCK_LINEAR_SHORT *kernel_matrix
	 
	)
	
	{




	uchar valid, finish;
	
	ushort k = get_local_id(0);
	ushort k1=(k+1) % info->block_size;

	ushort loc_size = info->block_size;

	int x, y, i, ib, bx, by, x1, y1;
	i=0;
	ib=0-k;
	x=0-k;
	y=k;	
	bx=by=0;


	
	// Copy reference into local memory
	
	short a = ctrl->ref_complete_blocks;
	short b = ctrl->ref_pad_block_copiers;
	short c = ctrl->ref_single_pad_copiers;

	charBlockCopy16(gl_reference, reference, loc_size, k, a, b, c);


	// Copy query into local memory

	a = ctrl->qry_complete_blocks;
	b = ctrl->qry_pad_block_copiers;
	c = ctrl->qry_single_pad_copiers;

	charBlockCopy16(gl_query, query, loc_size, k, a, b, c);


	// Copy sub_matrix into local memory

	a = ctrl->sub_complete_copiers;
	b = ctrl->sub_complete_elements;
	c = ctrl->sub_pad_elements;
	
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
	

	

	// Test
	/*
	for(i=0; i < ctrl->n_ref; i++ ){
		out_buffer[k * t + i] = (float)loc_reference[k * t + i] + (float)k/1000;
	}*/

	// Copy INFO structure in private memory
	int loc_block_x, loc_block_y, loc_block_size;
	char d1, d2;
	uchar n;
	n = (uchar) info->symbol_number;
	d1 = (char) info->penality1;
	d2 = (char) info->penality2;
	loc_block_x = info->block_x;
	loc_block_y = info->block_y;
	loc_block_size = info->block_size;

	
	int h, h_max;
	
	uchar r, q;	// reference and query amminoacids
	char s;		// substitution matrix value
	
	int up, left, left1, diagonal;

	/* 
		state: 0 ok, 1 cambio blocco, 2 cambio riga, 3 fine	
	*/
	uchar state=0;
	
	finish = 0;
	
	h_max=0;
	i=0;

	barrier(CLK_LOCAL_MEM_FENCE);
	//barrier(CLK_GLOBAL_MEM_FENCE);


	/********************   SW Algorithm   **************************/
	
	while(finish == 0)
	{	
		valid = 1;
		if (ib<0 || state==3)
			valid = 0;

		left1 = kernel_matrix[k].left;
		
		h=0;

		if(valid == 1) {
			r = reference[y];
			q = query[x];
						
			if (q > 0 && q <= (int)n && r>0 && r <= (int)n){
				

				s = (char)sub_matrix[ r * (n + 2) + q ];
				
				left	=	kernel_matrix[k].left;
				up		=	kernel_matrix[k].up;
				diagonal=	kernel_matrix[k].diagonal;

				h=max(h, up + d1);
				h=max(h, left + d1);
				h=max(h, diagonal + s);

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
				x = bx * loc_block_size;
				y = by * loc_block_size + k;					
			}
			else{
				state = 0;
				x += 1;
			}
		}
						
		barrier(CLK_LOCAL_MEM_FENCE);	
		//barrier(CLK_GLOBAL_MEM_FENCE);

		if(valid == 1){
			if (state == 0){
				kernel_matrix[k].up=h;		
				kernel_matrix[k].diagonal=left1;
			
				kernel_matrix[k1].left=h;
			}
			else if (state == 1){
				aux_row[y1].h=h;
				kernel_matrix[k].up=aux_row[y].h;
				kernel_matrix[k].diagonal=aux_row[y-1].h;
						
				kernel_matrix[k1].left=h;
			}
			else if (state == 2){
				aux_row[y1].h=h;				
				kernel_matrix[k].up=aux_row[y].h;			
				if (y==0) kernel_matrix[k].diagonal=0;
				else kernel_matrix[k].diagonal=aux_row[y-1].h;				
							
				if (y==(loc_block_y * loc_block_size - 1));
				else kernel_matrix[k1].left=h;				
			}
			else if((valid == 1) && state==3){
				kernel_matrix[k1].left=h;
			}
		}		
		barrier(CLK_LOCAL_MEM_FENCE);	
		//barrier(CLK_GLOBAL_MEM_FENCE);

		i+=1;
		if(i==(loc_block_x * loc_block_y * loc_block_size + loc_block_size - 1))
			finish = 1;
	}
	kernel_matrix[k].h_max = h_max;
	
	// In-loco search of maximum alignment score

	
	while(loc_size > 2){
		loc_size = loc_size / 2;
		barrier(CLK_LOCAL_MEM_FENCE);
		if(k < loc_size){
			kernel_matrix[k].h_max = max(kernel_matrix[k].h_max, kernel_matrix[k + loc_size - 1].h_max);	// should be replaced with SELECT, more efficient
		}
	}

	

	// Copy back data from local to global memory
	if(k == 0)
		*al_score = (int)kernel_matrix[0].h_max;

}
