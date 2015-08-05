/// 9 int -> 36 Byte
// Reference Lenght * 36 Byte 
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

/// 12 int -> 48 Byte 
// Workers * 48 Byte
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
} BLOCK;

/// 8 int -> 32 Byte
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

__kernel void politosw(
	__global int *reference, 
	__global int *query, 
	__global BLOCK_AUX *aux_row, 
	__global int *sub_matrix, 
	__global BLOCK *kernel_matrix, 
	__global INFO *info)
	{
	int valid, finish;
	/* 
		k:	id del kernel
		k1:	id del kernel sucessivo 
	*/
	int k = get_local_id(0);
	int k1=(k+1) % info->block_size;

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
		d1:	penalità 1
		d2: penalità 2		
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

		up,left,digonal: valori per il calcolo di h
	*/
	int h, h_max, r, q, s;
	int up, left, left1, diagonal;

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
		
		h=0;
		if(valid) {
			r= reference[y];
			q= query[x];				
			if (q>0 && q<=n && r>0 && r<=n){				
				s=sub_matrix[r*(n+2)+q];		    
				
				left	=	kernel_matrix[k].left;
				up		=	kernel_matrix[k].up;
				diagonal=	kernel_matrix[k].diagonal;

				h=max(h,up+d1);
				h=max(h,left+d1);
				h=max(h,diagonal+s);

				h_max = max(h_max, h);			
			}		
		}

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
						
		barrier(CLK_LOCAL_MEM_FENCE);	
		barrier(CLK_GLOBAL_MEM_FENCE);
		
		if(valid){
			if (state==0){
				kernel_matrix[k].up=h;
				kernel_matrix[k].diagonal=left1;
			
				kernel_matrix[k1].left=h;
			}
			else if (state==1 ){
				aux_row[y1].h=h;
				kernel_matrix[k].up=aux_row[y].h;
				kernel_matrix[k].diagonal=aux_row[y-1].h;
						
				kernel_matrix[k1].left=h;
			}
			else if (state==2 ){
				aux_row[y1].h=h;	
				kernel_matrix[k].up=aux_row[y].h;			
				if (y==0){
					kernel_matrix[k].diagonal=0;
				}else{
					kernel_matrix[k].diagonal=aux_row[y-1].h;				
				}
							
				if (y != (info->block_y*info->block_size-1)){
					kernel_matrix[k1].left=h;
				}
			}
			else if(valid && state==3){
				kernel_matrix[k1].left=h;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);	
		barrier(CLK_GLOBAL_MEM_FENCE);

		i+=1;
		if(i==(info->block_x*info->block_y*info->block_size+info->block_size-1))
			finish=1;
	}
	kernel_matrix[k].h_max = h_max;
}
