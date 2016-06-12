//====================================================================================================100
//		UPDATE
//====================================================================================================100

//    2006.03   Rob Janiczek
//        --creation of prototype version
//    2006.03   Drew Gilliam
//        --rewriting of prototype version into current version
//        --got rid of multiple function calls, all code in a  
//         single function (for speed)
//        --code cleanup & commenting
//        --code optimization efforts   
//    2006.04   Drew Gilliam
//        --added diffusion coefficent saturation on [0,1]
//		2009.12 Lukasz G. Szafaryn
//		-- reading from image, command line inputs
//		2010.01 Lukasz G. Szafaryn
//		--comments

//====================================================================================================100
//	DEFINE / INCLUDE
//====================================================================================================100

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "define.c"
#include "graphics.c"
#include "resize.c"
#include "timer.c"

FILE *fileResults;

/**Returns the Mean Squared Error based on the distance between the resulting golden image and the approximation
 * @params Ne, number of elements
 * @params img1, golden image
 * @params img2, approximation image
 **/
float image_distance(int Ne, fp* img1, fp* img2){
    
    int i;
    float mse = 0;
    
    for (i = 0; i < Ne; i++){
        mse += pow((int)img1[i] - (int)img2[i],2);
    }
    
    return mse/Ne;
}

int main(int argc, char *argv []){

    //VARIABLES
    // time
    long long time0;
    long long time1;
    long long time2;
    long long time3;
    long long time4;
    long long time5;
    long long time6;
    long long time7;
    long long time8;
    long long time9;
    long long time10;
    
    time0 = get_time();

    // inputs image, input paramenters
    fp* image_ori;																// originalinput image
    int image_ori_rows;
    int image_ori_cols;
    long image_ori_elem;

    // inputs image, input paramenters
    char   *filename = 0;
    fp* image;                  //input image
    long Nr,Nc;			// IMAGE nbr of rows/cols/elements
    long Ne;                    //number of elements

    // algorithm parameters
    int niter;			// nbr of iterations
    fp lambda;			// update step size

    // size of IMAGE
    int r1,r2,c1,c2;		// row/col coordinates of uniform ROI
    long NeROI;			// ROI nbr of elements

    // ROI statistics
    fp meanROI, varROI, q0sqr;					//local region statistics

    // surrounding pixel indicies
    int *iN,*iS,*jE,*jW;    

    // center pixel value
    fp Jc;

    // directional derivatives
    fp *dN,*dS,*dW,*dE;
    
    // calculation variables
    fp tmp,sum,sum2;
    fp G2,L,num,den,qsqr,D;

    // diffusion coefficient
    fp *c; 
    fp cN,cS,cW,cE;
    
    // counters
    long i,j;    // image row/col
    long k;      // image single index    

    // number of threads
    int threads;

    //Approximation: vars added
    char *config_filename = 0;
    char line[1024]; 
    int numApprox = 0;
    int row_top_idx, row_bottom_idx, col_left_idx, col_right_idx;
    fp* previmg;                //copy of the original image
    fp* golden;
    int delta = 0;
    int *loops;
    float *lambdas;
    float *thresholds_pix_changed;
    float *thresholds_prct_changed;
    float *img_areas;
    double  timing;
    
    time1 = get_time();

    //GET INPUT PARAMETERS

    if(argc != 6){
        printf("ERROR: wrong number of arguments\n");
	return 0;
    } else {
	filename = argv[1];
        config_filename = argv[2];
	Nr = atoi(argv[3]);				// it is 502 in the original image: Rows
        Nc = atoi(argv[4]);				// it is 458 in the original image: Columns
	threads = atoi(argv[5]);                        // Number of threads
    }
        
/*
    
    */
    // Reading Configuration File
        
    FILE *infile;
    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: no such file (%s)\n", filename);
        return 0;
    }
    
    FILE *infileconfig;
    if ((infileconfig = fopen(config_filename, "r")) == NULL) {
        fprintf(stderr, "Error: no such file (%s)\n", config_filename);
        return 0;
    }
        
    while (fgets(line, 1024, infileconfig) != NULL) {
        if (line[0] == '#') {
            continue;
        } else { 
            if (!strcmp(line, "\n") || strcmp(line,"\r\n")){//strcmp verifies line has more than only \n
                if (strtok(line, " \t\n") != 0) 
                    numApprox++;
            } else {
                continue;
            }
        }
    }
    rewind(infileconfig);
        
    printf("\nNumber of computations: %d\n",numApprox);
        
    if (numApprox < 2) {
        printf("\nError: You need to specify one configuration set for the 'golden' and at least one for the approximation\n");
        return 0;
    }
   
    loops = (int*)malloc(numApprox * sizeof(int));
    lambdas = (float*)malloc(numApprox * sizeof(float));
    thresholds_pix_changed = (float*)malloc(numApprox * sizeof(float));
    thresholds_prct_changed = (float*)malloc(numApprox * sizeof(float));
    img_areas = (float*)malloc(numApprox * sizeof(float)); 
    
    i = 0;

    while (fgets(line, 1024, infileconfig) != NULL) {
        if (strtok(line, " \t\n") == NULL) continue; 
        
        if (line[0] == '#') {
            continue;
        } else  if ( i < numApprox ) { 
            int loop = atoi(strtok(NULL, " ,\t\n")); 
            if (loop <= 0){
                printf("Error: Parameter 'Number of iterations' must be greater than 0\n");
                return 0;
            }
            loops[i] = loop;

            float lambda = atof(strtok(NULL, " ,\t\n")); 
            if (lambda < 0  || lambda > 1){
                printf("Error: Parameter 'Lambda' must be between 0 and 1\n");
                return 0;
            }
            lambdas[i] = lambda;
            
            float pixc = atof(strtok(NULL, " ,\t\n")); 
            if (pixc < 0  || pixc > 1){
                printf("Error: Parameter 'Percentage of changed pixel in image' must be between 0 and 1\n");
                return 0;
            }
            thresholds_pix_changed[i] = pixc;
            
            float change = atof(strtok(NULL, " ,\t\n")); 
            if (change < 0  || change > 1){
                printf("Error: Parameter 'Percentage of change in image' must be between 0 and 1\n");
                return 0;
            }
            thresholds_prct_changed[i] = change;

            float area = atof(strtok(NULL, " ,\t\n")); 
            if (area < 0  || area > 1){
                printf("Error: Parameter 'Percentage of image area to compute' must be between 0 and 1\n");
                return 0;
            }
            img_areas[i] = area;
            i++;
        }    
    }
        
    fclose(infileconfig);
    
    omp_set_num_threads(threads);

    time2 = get_time();

    //READ IMAGE (SIZE OF IMAGE HAS TO BE KNOWN)
    /*
    image_ori_rows = 502;   // Fixed 1704-502
    image_ori_cols = 458;   // Fixed 2272-458
    */
    
    // Alternative, reading from the PGM file
    i = 0;
    
    while (fgets(line, 1024, infile) != NULL && i < 3) {
        if (line[0] == '#') {
            continue;
        } else {
            //The rows and columns are specified in the line # 1
            if (i == 1){
                const char* val1 = strtok(line, " \t");
                const char* val2 = strtok(NULL, " \t");
                image_ori_cols = atoi(val1);
                image_ori_rows = atoi(val2);
                printf("Original #Cols # = %d\n", image_ori_cols);
                printf("Original #Rows # = %d\n", image_ori_rows);
            } else {
                strtok(line, " \t\n");
            }
            i++;
        }
    }
    rewind(infile);
    
    image_ori_elem = image_ori_rows * image_ori_cols;   //Matrix

    image_ori = (fp*)malloc(sizeof(fp) * image_ori_elem); //Allocate memory for matrix

    read_graphics(filename, image_ori, image_ori_rows, image_ori_cols, 1);

    time3 = get_time();     // Reading graphics time

    //RESIZE IMAGE (ASSUMING COLUMN MAJOR STORAGE OF image_orig)

    Ne = Nr*Nc;

    image = (fp*)malloc(sizeof(fp) * Ne);
    previmg = (fp*)malloc(sizeof(fp) * Ne);
    golden = (fp*)malloc(sizeof(fp) * Ne);
    
    resize(image_ori, image_ori_rows, image_ori_cols, image, previmg, Nr, Nc, 1);

    time4 = get_time();

    //SETUP
    r1 = 0;                                                 // top row index of ROI, Regions of Interests?
    r2 = Nr - 1;                                            // bottom row index of ROI
    c1 = 0;                                                 // left column index of ROI
    c2 = Nc - 1;                                            // right column index of ROI

    // ROI image size    
    NeROI = (r2-r1+1)*(c2-c1+1);				// number of elements in ROI, ROI size
    
    // allocate variables for surrounding pixels
    iN = malloc(sizeof(int*)*Nr) ;				// north surrounding element
    iS = malloc(sizeof(int*)*Nr) ;				// south surrounding element
    jW = malloc(sizeof(int*)*Nc) ;				// west surrounding element
    jE = malloc(sizeof(int*)*Nc) ;				// east surrounding element
    
    // allocate variables for directional derivatives       
    dN = malloc(sizeof(fp)*Ne) ;				// north direction derivative
    dS = malloc(sizeof(fp)*Ne) ;				// south direction derivative
    dW = malloc(sizeof(fp)*Ne) ;				// west direction derivative
    dE = malloc(sizeof(fp)*Ne) ;				// east direction derivative

    // allocate variable for diffusion coefficient
    c  = malloc(sizeof(fp)*Ne) ;				// diffusion coefficient
        
    // N/S/W/E indices of surrounding pixels (every element of IMAGE)
    // #pragma omp parallel
    for (i=0; i<Nr; i++) {
        iN[i] = i-1;						// holds index of IMAGE row above
        iS[i] = i+1;						// holds index of IMAGE row below
    }
    // #pragma omp parallel
    for (j=0; j<Nc; j++) {
        jW[j] = j-1;						// holds index of IMAGE column on the left
        jE[j] = j+1;						// holds index of IMAGE column on the right
    }
    // N/S/W/E boundary conditions, fix surrounding indices outside boundary of IMAGE
    iN[0]    = 0;						// changes IMAGE top row index from -1 to 0
    iS[Nr-1] = Nr-1;						// changes IMAGE bottom row index from Nr to Nr-1 
    jW[0]    = 0;						// changes IMAGE leftmost column index from -1 to 0
    jE[Nc-1] = Nc-1;						// changes IMAGE rightmost column index from Nc to Nc-1

    //COMPUTATION
    
    int loop = 0;
    int m = 0;
    
    // Result file
    char buffer [255];
    sprintf(buffer,"../results_%lld.txt",time0);
        
    fileResults = fopen(buffer, "w");
        
    if (fileResults == NULL){
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(fileResults, "\n┌———————————————————————————————————————————————— GLOBAL CONFIGURATIONS —————————————————————————————————————————————————┐");
    fprintf(fileResults, "\n│%36s%-35s%-1s%16d%35s"," ","Number of Rows","=", Nr,"│");
    fprintf(fileResults, "\n│%36s%-35s%-1s%16d%35s"," ","Number of Columns","=", Nc,"│");
    fprintf(fileResults, "\n│%36s%-35s%-1s%16d%35s"," ","Number of Elements","=", Ne,"│");
    fprintf(fileResults, "\n│%36s%-35s%-1s%16d%35s"," ","Number of Threads","=", threads,"│");
    fprintf(fileResults, "\n└————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————┘");
    
    for (m = 0; m < numApprox; m++) {
        
        fprintf(fileResults, "\n\n\n\n ××××××××××××××××××××××××××××××  SETTING N° %d  ×××××××××××××××××××××××××××××",m+1);
        fprintf(fileResults, "\n×%15s%-31s%-1s%16d%14s"," ","Number of Loops","=", loops[m],"×");
        fprintf(fileResults, "\n×%15s%-31s%-1s%16f%14s"," ","Lambda","=",lambdas[m],"×");
        fprintf(fileResults, "\n×%15s%-31s%-1s%16f%14s"," ","Max. Change in pixel","=", thresholds_prct_changed[m],"×");
        fprintf(fileResults, "\n×%15s%-31s%-1s%16f%14s"," ","Max. Pixels changed","=", thresholds_pix_changed[m],"×");            
        fprintf(fileResults, "\n×%15s%-31s%-1s%16f%14s"," ","Image area to compute","=", img_areas[m],"×");            
        fprintf(fileResults, "\n ×××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××\n\n");
        
        time5 = get_time();
        //Reseting and scaling image based on the initial values. FROM 0-255 TO 0-1 AND EXTRACT
        int p;
        for (p=0; p < Ne; p++) {
            image[p] = exp(previmg[p]/255);     // exponentiate input IMAGE and copy to output image
        }
        time6 = get_time();
        
        // Are we going to tun the algorithm in the full image? or only in the 80% of it, for example.
        //Row boundaries, later we will loop the image from row_top_idx to row_bottom_idx    
        row_top_idx = lround(0.5*Nr*(1-img_areas[m]));
        row_bottom_idx = Nr - row_top_idx;

        //Col/row boundaries, later we will loop the image from col_left_idx to col_right_idx
        col_left_idx = lround(0.5*Nc*(1-img_areas[m]));
        col_right_idx = Nc - col_left_idx;

        fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");
        fprintf(fileResults,"\t\t\t\t\t\tComputing SRAD in Rows[%d,%d] Columns[%d,%d]\n", row_top_idx, row_bottom_idx, col_left_idx, col_right_idx);
        fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");
        
        loop = 0;
        timing = omp_get_wtime();

        
        do {
            delta = 0;
            // ROI statistics for entire ROI (single number for ROI)
            sum = 0; 
            sum2 = 0;
            for (i = r1; i <= r2; i++) {				// do for the range of rows in ROI
                for (j = c1; j <= c2; j++) {			// do for the range of columns in ROI
                    tmp   = image[i + Nr*j];			// get coresponding value in IMAGE
                    sum  += tmp ;				// take corresponding value and add to sum
                    sum2 += tmp*tmp;				// take square of corresponding value and add to sum2
                }
            }
            meanROI = sum / NeROI;				// gets mean (average) value of element in ROI
            varROI  = (sum2 / NeROI) - meanROI*meanROI;		// gets variance of ROI
            q0sqr   = varROI / (meanROI*meanROI);		// gets standard deviation of ROI
            
            
            #pragma omp parallel \
            shared(image, dN, dS, dW, dE, c, Nr, Nc, iN, iS, jW, jE, col_right_idx, col_left_idx, row_bottom_idx, row_top_idx) 
            {
                #pragma omp for \
                        private(i, j, k, Jc, G2, L, num, den, qsqr) \
                        schedule(static)

                for (j = col_left_idx; j < col_right_idx; j++) {    // do for the range of columns in IMAGE

                    for (i = row_top_idx; i < row_bottom_idx; i++) {// do for the range of rows in IMAGE 

                        // current index/pixel
                        k = i + Nr*j;				// get position of current element
                        Jc = image[k];				// get value of the current element

                        // directional derivates (every element of IMAGE)
                        dN[k] = image[iN[i] + Nr*j] - Jc;	// north direction derivative
                        dS[k] = image[iS[i] + Nr*j] - Jc;	// south direction derivative
                        dW[k] = image[i + Nr*jW[j]] - Jc;	// west direction derivative
                        dE[k] = image[i + Nr*jE[j]] - Jc;	// east direction derivative

                        // normalized discrete gradient mag squared (equ 52,53)
                        G2 = (dN[k]*dN[k] + dS[k]*dS[k]			// gradient (based on derivatives)
                            + dW[k]*dW[k] + dE[k]*dE[k]) / (Jc*Jc);

                        // normalized discrete laplacian (equ 54)
                        L = (dN[k] + dS[k] + dW[k] + dE[k]) / Jc;	// laplacian (based on derivatives)

                        // ICOV (equ 31/35) Instantaneous coeffcient of variation
                        num  = (0.5*G2) - ((1.0/16.0)*(L*L)) ;		// num (based on gradient and laplacian)
                        den  = 1 + (.25*L);				// den (based on laplacian)
                        qsqr = num/(den*den);				// qsqr (based on num and den)

                        // diffusion coefficent (equ 33) (every element of IMAGE)
                        den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)) ;	// den (based on qsqr and q0sqr)
                        c[k] = 1.0 / (1.0+den) ;			// diffusion coefficient (based on den)

                        // saturate diffusion coefficent to 0-1 range
                        if (c[k] < 0)					// if diffusion coefficient < 0
                                                {c[k] = 0;}												// ... set to 0
                        else if (c[k] > 1)				// if diffusion coefficient > 1
                                                {c[k] = 1;}												// ... set to 1
                    }
                }            
            }   
            
            lambda = lambdas[m];
            
            #pragma omp parallel \
            shared(image, c, Nr, Nc, lambda, previmg, col_right_idx, col_left_idx, row_bottom_idx, row_top_idx)
            {
		
                #pragma omp for \
                            private(i, j, k, D, cS, cN, cW, cE) \
                            schedule(static) \
                            reduction(+:delta)

                for (j=col_left_idx; j<col_right_idx; j++) {					// do for the range of columns in IMAGE

                    for (i=row_top_idx; i<row_bottom_idx; i++) {				// do for the range of rows in IMAGE

                        // current index
                        k = i + Nr*j;					// get position of current element

                        // diffusion coefficent
                        cN = c[k];					// north diffusion coefficient
                        cS = c[iS[i] + Nr*j];				// south diffusion coefficient
                        cW = c[k];					// west diffusion coefficient
                        cE = c[i + Nr*jE[j]];				// east diffusion coefficient

                        // divergence (equ 58)
                        D = cN*dN[k] + cS*dS[k] + cW*dW[k] + cE*dE[k];	// divergence

                        // image update (equ 61) (every element of IMAGE)
                        image[k] = image[k] + 0.25*lambda*D;		// updates image (based on input time step, "lambda" and divergence)

                        //If the pixel has changed more than the max percentage change permitted and stablished on the
                        //config file, increase delta
                        if (abs(previmg[k] - (int)log(image[k])*255) > thresholds_prct_changed[m] * previmg[i]) delta += 1;

                    }
                }
            }
            
        } while (delta <= thresholds_pix_changed[m]*Ne && loop++ < loops[m]);
        
        timing = omp_get_wtime() - timing;
        
        fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");
        fprintf(fileResults,"\t\t\t\t\t\tDelta = %d vs Threshold = %f\n", delta, thresholds_pix_changed[m]*Ne);
        fprintf(fileResults,"\t\t\t\t\t\tExec. Loops = %d vs Config. Loops = %d\n", loop, loops[m]);
        fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");        
        fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");
        fprintf(fileResults,"\t\t\t\t\t\tProcess Time: %f\n", timing);
        fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");

        time7 = get_time();
        
        for (i=0; i < Ne; i++) {		// do for the number of elements in IMAGE
            image[i] = log(image[i])*255;	// take logarithm of image, log compress
        }
        
        time8 = get_time();
        
        if ( m == 0){
            //If we're running the golden, assumed always be set as the first configuration line in 
            //the config file
            memcpy(golden, image, sizeof(fp) * Ne);
            
        } else {
            fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");
            fprintf(fileResults,"\t\t\t\t\t\tImage difference: %f\n", image_distance(Ne, golden, image));
            fprintf(fileResults,"…………………………………………………………………………………………………………………………………………………………………………………………………………\n");
        }
        
        //WRITE IMAGE AFTER PROCESSING
        sprintf(buffer,"../image_out_%d_%lld.pgm", m, time0);
        write_graphics(buffer, image, Nr, Nc, 1, 255);
        time9 = get_time();

    } 
    
    fclose(fileResults);
    
    //DEALLOCATE
    free(image_ori);
    free(image);
    free(previmg);
    free(iN); free(iS); free(jW); free(jE);			// deallocate surrounding pixel memory
    free(dN); free(dS); free(dW); free(dE);			// deallocate directional derivative memory
    free(c);							// deallocate diffusion coefficient memory

    time10 = get_time();

    //DISPLAY TIMING

    printf("Time spent in different stages of the application:\n");
    printf("%.12fs, %.12f : SETUP VARIABLES\n",(float) (time1-time0) / 1000000, (float) (time1-time0) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : READ COMMAND LINE PARAMETERS\n",(float) (time2-time1) / 1000000, (float) (time2-time1) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : READ IMAGE FROM FILE\n",(float) (time3-time2) / 1000000, (float) (time3-time2) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : RESIZE IMAGE\n", 	(float) (time4-time3) / 1000000, (float) (time4-time3) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : SETUP, MEMORY ALLOCATION\n",(float) (time5-time4) / 1000000, (float) (time5-time4) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : EXTRACT IMAGE\n",(float) (time6-time5) / 1000000, (float) (time6-time5) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : COMPUTE\n",(float) (time7-time6) / 1000000, (float) (time7-time6) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : COMPRESS IMAGE\n",(float) (time8-time7) / 1000000, (float) (time8-time7) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : SAVE IMAGE INTO FILE\n",(float) (time9-time8) / 1000000, (float) (time9-time8) / (float) (time10-time0) * 100);
    printf("%.12fs, %.12f : FREE MEMORY\n",(float) (time10-time9) / 1000000, (float) (time10-time9) / (float) (time10-time0) * 100);
    printf("Total time:\n");
    printf("%.12f s\n", 																					(float) (time10-time0) / 1000000);

//END OF FILE
}


