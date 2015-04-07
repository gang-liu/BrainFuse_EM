/*==============================================
 * Runs the EM algorithm for the label fusion
 *==============================================*/
/* $Revision: 1.10.6.6 $ */
#include <math.h>
#include "mex.h"


/* Constants */
#define PI  3.141592653589793238462643383279502884197169399375

/* Input Arguments */

#define	I 	   prhs[0]
#define	IREG   prhs[1]
#define	LREG   prhs[2]
#define	VAR    prhs[3]
#define	MAPS   prhs[4]
#define	INDS   prhs[5]


/* Output Arguments */

#define	W	   plhs[0]
#define	MUS    plhs[1]
#define	PRIORS plhs[2]
#define	POSTS  plhs[3]
#define	SEG    plhs[4]

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
        
{
    
    // IO, constants, etc
    
    int Nscans = mxGetN(IREG);
    const int *dimsMap = mxGetDimensions(MAPS);
    const int Nmaps = *(dimsMap+0);
    const int Nlabels = *(dimsMap+1);
    
    float *pI = (float *) mxGetData(I);
    float *pIREG = (float *) mxGetData(IREG);
    unsigned char *pLREG = (unsigned char *) mxGetData(LREG);
    float *pVAR = (float *) mxGetData(VAR);
    unsigned char *pMAPS = (unsigned char *) mxGetData(MAPS);
    unsigned char *pINDS = (unsigned char *) mxGetData(INDS);
    
    
    mwSize dim[2];
    dim[0]=Nlabels; dim[1]=Nscans;
    W = mxCreateNumericArray( 2, dim, mxSINGLE_CLASS, mxREAL);
    dim[0]=Nlabels; dim[1]=1;
    MUS = mxCreateNumericArray( 2, dim, mxSINGLE_CLASS, mxREAL);
    PRIORS = mxCreateNumericArray( 2, dim, mxSINGLE_CLASS, mxREAL);
    POSTS = mxCreateNumericArray( 2, dim, mxSINGLE_CLASS, mxREAL);
    dim[0]=1; dim[1]=1;
    SEG = mxCreateNumericArray( 2, dim, mxUINT8_CLASS, mxREAL);
    
    
    float *pW = (float *) mxGetData(W);
    float *pMUS = (float *) mxGetData(MUS);
    float *pPRIORS = (float *) mxGetData(PRIORS);
    float *pPOSTS = (float *) mxGetData(POSTS);
    unsigned char *pSEG = (unsigned char *) mxGetData(SEG);
    
    
    float gaussConstant = 1.0/sqrt(2*PI*(*pVAR));
    float prec = 1.0 / (*(pVAR));
    
    
    // Initialization
    
    double aux = 1.0/((double) Nlabels);
    for(int l=0; l < Nlabels; l++)
    {
        *(pPRIORS+l) = aux;
        *(pMUS+l) = 50.0;
    }
    
    
    // EM algorithm
    
    bool ready = false;
    float cost = 0;
    float normalizer = 0;
    float invNormalizer = 0;
    float dif = 0;
    float prob = 0;
    float weight;
    float costOld = 1e30;
    int its = 0;
    
    while(!ready)
    {
        its++;
        cost = 0;
        
        // E step
        for(int m=0; m<Nscans; m++)
        {
            
            normalizer = 1e-12;
            for(int l=0; l < Nlabels; l++)
            {
                if(*( pMAPS + (*(pINDS+m))-1 + (*(pLREG+m)-1)*Nmaps + l*Nmaps*Nlabels)>0)
                {
                    dif = *(pIREG+m) - *(pMUS+l);
                    prob = *(pPRIORS+l) * gaussConstant * exp(-0.5 * dif * dif * prec );
                    *(pW+l+m*Nlabels) = prob;
                    normalizer += prob; 
                }
                else
                {
                    *(pW+l+m*Nlabels) = 0;
                }
            }
            cost -= log(normalizer);
            invNormalizer = 1.0 / normalizer;
            for(int l=0; l < Nlabels; l++)   *(pW+l+m*Nlabels)=*(pW+l+m*Nlabels)*invNormalizer;
        }
        // printf("Iteration %d, cost %f\n",its,cost);
        
        
        // M step
        for(int l=0; l < Nlabels; l++)
        {
            *(pMUS+l) = 0;
            normalizer = 1e-12;
            for(int m=0; m<Nscans; m++)
            {
                weight = *(pW+l+m*Nlabels);
                *(pMUS+l) = *(pMUS+l) + weight * (*(pIREG+m));
                normalizer += weight;
            }
            *(pMUS+l) = *(pMUS+l) / normalizer;
            *(pPRIORS+l) = normalizer / ((float) Nscans);
        }
        
        if( its>200 || (costOld-cost)/cost<1e-3)
        {
            //  printf("Converged!\n",its,cost);
            ready=true;
        }
        else
        {
            costOld=cost;
        }
        
    }
    
    
    // computation of label posteriors
    // computation of label posteriors
    float *C = (float *) malloc( Nscans * sizeof(float));
    float *p_Im_l = (float *) malloc( Nscans * Nlabels * sizeof(float));
    for(int m=0; m<Nscans; m++)
    {
        for(int l=0; l < Nlabels; l++)
        {
            if(*( pMAPS + (*(pINDS+m))-1 + (*(pLREG+m)-1)*Nmaps +l*Nmaps*Nlabels)>0)
            {
                dif = *(pIREG+m) - *(pMUS+l);
                prob = gaussConstant * exp(-0.5 * dif * dif * prec );
                *(p_Im_l+l+m*Nlabels) = prob;
                *(C+m) = (*(C+m)) + prob * (*(pPRIORS+l));
            }
        }
    }


    normalizer = 1e-12;
    float max=-1e20;
    unsigned char labMax=-1;
    for(int l=0; l < Nlabels; l++)
    {
        *(pPOSTS+l)=0;
        for(int m=0; m<Nscans; m++)
        {
            if(*( pMAPS + (*(pINDS+m))-1 + (*(pLREG+m)-1)*Nmaps +l*Nmaps*Nlabels)>0)
            {
                 *(pPOSTS+l) = *(pPOSTS+l) + *(p_Im_l+l+m*Nlabels) /( 1e-12 + *(C+m) );  /* C_m \propto 1/Psi_m */
            }
        }
        if ( *(pPOSTS+l) > 0)
        {
            dif = *(pI) - *(pMUS+l);
            *(pPOSTS+l) = (*(pPOSTS+l)) * (*(pPRIORS+l)) * exp(-0.5 *dif * dif * prec );
            normalizer += (*(pPOSTS+l));
            if ((*(pPOSTS+l))>max)
            {
                labMax = l;
                max = (*(pPOSTS+l));
            }
        }
    }
    *pSEG = 1+labMax;
    invNormalizer = 1.0 / normalizer;
    for(int l=0; l < Nlabels; l++)
    {
        *(pPOSTS+l) = *(pPOSTS+l) * invNormalizer;


    }

    
    free(C);
    free(p_Im_l);
    
    
    return;

}


