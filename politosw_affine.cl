
typedef struct _block_aux {
	int h;
	int e;
	int f;
	int e_diagonal;
	int f_diagonal;
	int e_aux;
	int f_aux;
	int gSbj;
	int gQry;
} BLOCK_AUX;

typedef struct _block {
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
	int gSbj;
	int gQry;
    #ifdef DEBUG 
		BLOCK_DEBUG debug;
	#endif
} BLOCK;

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

__kernel void politosw(
	__global int *reference, 
	__global int *query, 
	__global BLOCK_AUX *aux_row, 
	__global BLOCK_AUX *aux_column,
	__global int *sub_matrix, 
	__global BLOCK *kernel_matrix, 
	__global INFO *info)
	{
	int valid, finish;
	/* 
		k:	id del kernel
		k1:	id del kernel sucessivo 
	*/
	// int k=get_local_id(0);
	int k = get_global_id(0);
	int k1=(k+1) % info->block_size;
	#ifdef DEBUG
	kernel_matrix[k].debug.k=k;
	kernel_matrix[k].debug.k1=k1;
	#endif
	/* 
		i:		iterazione
		ib:		iterazione sul blocco

		x,y:	coordinate della matrice 
		bx,by:	coordinate del blocco
	*/
	int x, y, i, ib, bx, by, x1, y1;
	i=0;
	ib=0-k;
	x=0-k;
	y=k;	
	bx=by=0;

	/*
		n:	numero di simboli
		d1:	penalit� 1
		d2: penalit� 2		
	*/
	int d1, d2, n;
	n=info->symbol_number;
	d1=info->penality1;
	d2=info->penality2;
	
	/*
		h:	valore della matrice
		r:	valore del reference
		q:	valore della query
		s:	valore di sostituzione

		up,left,diagonal: valori per il calcolo di h
		e,f: valori usati per l'affine gap
	*/
	int h, h_max, r, q, s;
	int up, left, left1, diagonal;
	int e, f, e_aux, f_aux;

	/* 
		state: 0 ok, 1 cambio blocco, 2 cambio riga, 3 fine	
	*/
	int state=0;
	finish=0;
	h_max=0;
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
			#ifdef DEBUG	
			kernel_matrix[k].debug.all_h[i]=h;
			#endif
		}
		
		#ifdef debug		
		kernel_matrix[k].debug.all_up[i]=		kernel_matrix[k].up;
		kernel_matrix[k].debug.all_diagonal[i]=	kernel_matrix[k].diagonal;
		kernel_matrix[k].debug.all_left[i]=		kernel_matrix[k].left;
		kernel_matrix[k].debug.all_e[i]=		e;
		kernel_matrix[k].debug.all_f[i]=		f;
		#endif

		// Aggiornamento degli indici		
		if (state!=3){
			ib+=1;
			x1=x;
			y1=y;

			if (ib==info->block_size){
				//CAMBIO BLOCCO
				state=1;
				by+=1;
				ib=0;
				if (by==info->block_y){				
					//CAMBIO RIGA
					state=2;
					by=0;
					bx+=1;
					if (bx==info->block_x){
						//FINE
						state=3;					
					}
				}
				x=bx*info->block_size;
				y=by*info->block_size + k;					
			}
			else{
				state=0;
				x+=1;
			}
		}
		
		#ifdef DEBUG			
		kernel_matrix[k].debug.all_state[i]=state;
		#endif
								
		barrier(CLK_LOCAL_MEM_FENCE);	
		barrier(CLK_GLOBAL_MEM_FENCE);	
		if(valid){
			if (state==0){
				// Verso altri
				if (y1==(info->block_y*info->block_size-1));
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
				if (y1==(info->block_y*info->block_size-1));
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
		
		/*if(valid){
			if (state==0){
				kernel_matrix[k].up=h;
				kernel_matrix[k].f=f;
				kernel_matrix[k].diagonal=left1;
				
				if (y==(info->block_y*info->block_size-1));
				else {
					kernel_matrix[k1].left=h;
					kernel_matrix[k1].e=e;
				}
				//kernel_matrix[k1].left=h;
				//kernel_matrix[k1].e=e;
			}
			else if (state==1){
				aux_row[y1].h=h;
				aux_row[y1].e=e;
				aux_row[y1].f=f;
				kernel_matrix[k].up=aux_row[y].h;
				kernel_matrix[k].f=aux_row[y].f;
				kernel_matrix[k].diagonal=aux_row[y-1].h;
						
				kernel_matrix[k1].left=h;
				kernel_matrix[k1].e=e;
			}
			else if (state==2 ){
				aux_row[y1].h=h;
				aux_row[y1].e=e;
				aux_row[y1].f=f;				
				kernel_matrix[k].up=aux_row[y].h;
				kernel_matrix[k].f=aux_row[y].f;

				if (y==0) {
					kernel_matrix[k].diagonal=0;
					kernel_matrix[k].e = 0;
					kernel_matrix[k].left = 0;
				}
				else {
					kernel_matrix[k].diagonal=aux_row[y-1].h;
				}				
							
				if (y==(info->block_y*info->block_size-1));
				else {
					kernel_matrix[k1].left=h;
					kernel_matrix[k1].e=e;
				}
			}
			//else if(valid && state==3){ //Ridondante
			else if(state==3){
				kernel_matrix[k1].left=h;
				kernel_matrix[k1].e=e;
			}
		}
		*/		
		barrier(CLK_LOCAL_MEM_FENCE);	
		barrier(CLK_GLOBAL_MEM_FENCE);
		
		#ifdef DEBUG
		kernel_matrix[k].debug.all_up_after[i]=kernel_matrix[k].up;
		kernel_matrix[k].debug.all_diagonal_after[i]=kernel_matrix[k].diagonal;
		kernel_matrix[k].debug.all_left_after[i]=kernel_matrix[k].left;
		kernel_matrix[k].debug.all_e_after[i]=kernel_matrix[k].e;
		kernel_matrix[k].debug.all_f_after[i]=kernel_matrix[k].f;
		#endif

		i+=1;
		if(i==(info->block_x*info->block_y*info->block_size+info->block_size-1))
			finish=1;
	}
	#ifdef DEBUG
	kernel_matrix[k].debug.i=i;
	#endif
	kernel_matrix[k].h_max = h_max;
}