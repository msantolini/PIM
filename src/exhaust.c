#define GETVALUE(vble)  ( (vble) = getvalue(#vble) )
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "mex.h"
/*#include "rangen.h"
#include "rangen.c"*/

#define WRITEPROBA

int N;

int q=4;

typedef struct Var{
  int n;
  double entropy;
  double ****J;  double ****C;  double **H;  double **he;
  char type[10];
} Var;

int sub(char s1, char *s2) {
  int i;
  for(i=0;s2[i]!='\0';i++)
    if(s1==s2[i]) return 1;
  return 0;
}

void free2(double **tab, int n1) {
  int i;
  for(i=0;i<n1;i++) free(tab[i]);
  free(tab);
  return;
}

void free3(double ***tab,int n1, int n2) {
  int i;
  for(i=0;i<n1;i++) free2(tab[i],n2);
  free(tab);
}

void free4(double ****tab,int n1, int n2, int n3) {
  int i;
  for(i=0;i<n1;i++) free3(tab[i],n2,n3);
  free(tab);
}

double * d1(int n) {
  double *tab=calloc(n,sizeof(double));
  return tab;
}

double ** d2(int n1,int n2) {
  int i;
  double **tab=calloc(n1,sizeof(double*));
  for(i=0;i<n1;i++) tab[i]=d1(n2);
  return tab;
}

double *** d3(int n1,int n2,int n3) {
  int i;
  double ***tab=calloc(n1,sizeof(double**));
  for(i=0;i<n1;i++) tab[i]=d2(n2,n3);
  return tab;
}

double **** d4(int n1,int n2,int n3,int n4) {
  int i;
  double ****tab=calloc(n1,sizeof(double***));
  for(i=0;i<n1;i++) tab[i]=d3(n2,n3,n4);
  return tab;
}



Var * fill_var(int n,char *s) {
  int i;
  Var *m=malloc(sizeof(Var));
  m->n=n;
  strcpy(m->type,s);
  if(sub('H',s)) m->H=d2(n,q);      if(sub('e',s))  m->he=d2(n,q);   
  if(sub('C',s)) m->C=d4(n,n,q,q);    if(sub('J',s))  m->J=d4(n,n,q,q); 
  return m;
}

void free_var(Var *m) {
  int i,n=m->n;
  if(sub('H',m->type)) free2(m->H,n);      if(sub('e',m->type))  free2(m->he,n);   
  if(sub('C',m->type)) free4(m->C,n,n,q);    if(sub('J',m->type))  free4(m->J,n,n,q);
  free(m);
  return;
}

void copy_bits(char *bits1,char *bits2,int n) 
{  int i;  for(i=0;i<n;i++) bits1[i]=bits2[i];  return;}


int conv(int c,int n,char *bits) {
  int i;
  int nchg=0,ichg=0;
  for(i=0;i<n;i++) {
      if(bits[i]!=c%q) {
        nchg++;
        ichg=i;
        bits[i]=c%q;
      }
    c=c/q;
  }
  if(nchg>1) return n; else return n; // ???
}

double energie(char *bits,struct Var *var) {
  double sum=0;
  int i,j,n=var->n;
  for(i=0;i<n;i++) {
    sum+=var->he[i][bits[i]];
    for(j=0;j<i;j++)
      sum+=var->J[i][j][bits[i]][bits[j]];
  }
  return(sum);
}

void exhaust(Var *var,int nn, double *betas,double *me,double *me2, double *entropy, double *partfunc,double *offset, char *rawprobafile) {
  int a,c=0,i,j,n=var->n,g,h;
  /*for(i=0;i<n;i++) {
      for(g=0;g<q;g++)
          printf("%lg ",var->he[i][g]);
      printf("\n");
  }*/
  double *s1=calloc(nn,sizeof(double));
  double *entr=calloc(nn,sizeof(double));
  double s10=0,entr0=0;
  for(a=0;a<nn;a++) {
      s1[a]=0;
      entr[a]=0;
  }
  double w;
  for(i=0;i<n;i++) for(g=0;g<q;g++) var->H[i][g]=0;
  for(i=0;i<n;i++) for(j=0;j<n;j++) for(g=0;g<q;g++) for(h=0;h<q;h++) var->C[i][j][g][h]=0;  
  char *bits=calloc(n,sizeof(char));
  int nconfig=(int)pow(q,n);
#ifdef WRITEPROBA
  float *Proba;
  if(rawprobafile)
      Proba=calloc(nconfig, sizeof(float));
#endif
  *me=0;
  *me2=0;
  double ene;
  int ichg;
  for(c=0;c<nconfig;c++) {      
    /*printf("c = %d\n",c);fflush(stdout);*/
    ichg=conv(c,n,bits);
    ene=energie(bits, var);
    
    w=exp(ene);
    
    #ifdef WRITEPROBA
    if(rawprobafile)
        Proba[c]=ene;
    #endif
    s10+=w;
    entr0-=w*ene;
    for(a=0;a<nn;a++) {
      w=exp(betas[a]*(ene-offset[a]));
      entr[a]-=w*ene;
      s1[a]+=w;
      me2[a]+=w*ene*ene;
    }
    
    
    
    for(i=0;i<n;i++) var->H[i][bits[i]]+=w;
    for(i=0;i<n;i++) for(j=0;j<i;j++) var->C[i][j][bits[i]][bits[j]]+=w;

  }
#ifdef WRITEPROBA
  FILE *out;
  if(rawprobafile) {
    out=fopen(rawprobafile, "w");
    for(c=0;c<nconfig;c++) fprintf(out, "%lg\n", Proba[c]);
    fclose(out);
    free(Proba);
  }
#endif
  for(a=0;a<nn;a++) { 
    me[a]=-entr[a]/s1[a];
    me2[a]=me2[a]/s1[a]-me[a]*me[a];
    entropy[a]=log(s1[a])+betas[a]*offset[a]-me[a];
    partfunc[a]=s1[a];
  }
  /*var->entropy[0]=(entr/s1+log(s1))/log(2.);
    printf("Entropy: %lg\n",var->entropy[0]);*/
  for(i=0;i<n;i++) for(g=0;g<q;g++) {
        var->H[i][g]=var->H[i][g]/s10;
        for(h=0;h<q;h++) var->C[i][i][g][h]=var->C[i][i][h][g]=(h==g)?var->H[i][g]:0;
  }
  for(i=0;i<n;i++) for(j=0;j<i;j++) for(g=0;g<q;g++) for(h=0;h<q;h++) var->C[i][j][g][h]=var->C[j][i][h][g]=var->C[i][j][g][h]/s10;
  free(bits);
  free(s1);
  free(entr);
  return;
}




/* function adapted from Anton Krukowski's code */
double getvalue(char *vble)

{
    const mxArray *pscalar;
    double val;
    
    if ( (pscalar = mexGetVariablePtr("caller",vble)) == NULL) {
        mexPrintf("Couldn't find MATLAB matrix %s.\n",vble);
    }
    else {
        val = mxGetScalar(pscalar);
        return(val);
    }
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {

    double *Hs;
    double *Cs;
    double *Js;
    double *hs;
    double *me;
    
    /*GETVALUE(N);*/

    
        
    if (nrhs != 2 && nrhs != 3 && nrhs !=4 && nrhs !=5)
        mxErrMsgTxt("Two, three, four of five input variables required.\nh, J, betas, offsets, rawprobafileout.\n");
    if (nlhs != 3)
        mxErrMsgTxt("Three output variables required.");
    N=mxGetDimensions(prhs[0])[1];
    q=mxGetDimensions(prhs[0])[0];
    hs=mxGetPr(prhs[0]);
    Js=mxGetPr(prhs[1]);
    int i,j,g,h;
    
    Var *var=fill_var(N,"CHeJ");
    
    for(i=0;i<N;i++) for(g=0;g<q;g++) var->he[i][g]=hs[i*q+g];
    for(i=1;i<N;i++) for(j=0;j<i;j++) for(g=0;g<q;g++) for(h=0;h<q;h++)
        var->J[i][j][g][h]=var->J[j][i][h][g]=Js[i*N*q*q+j*q*q+g*q+h];
    
    
    
    int nn,a;
    double *betas, *offset;
    int noffset=0,nobeta=0;
    
    if(nrhs>2) {
      nn=mxGetDimensions(prhs[2])[0];
      if(nn==0) {
          nn=1; nobeta=1;
          betas=malloc(sizeof(double));
          betas[0]=1;
      } else {
          betas=mxGetPr(prhs[2]);
      }
    }
    else {
      betas=malloc(sizeof(double));
      betas[0]=1;
      nn=1; nobeta=1;
    }

    if(nrhs<=3) {
        noffset=1;
        offset=calloc(nn,sizeof(double));
        for(a=0;a<nn;a++) offset[a]=0;
    } else {
        if(mxGetDimensions(prhs[3])[0]==0) {
            noffset=1;
            offset=calloc(nn,sizeof(double));
            for(a=0;a<nn;a++) offset[a]=0;
        } else
            offset=mxGetPr(prhs[3]);
    }
        
    char *rawprobafile;
    rawprobafile=0;
    if(nrhs>4)
        rawprobafile=mxArrayToString(prhs[4]);
    
    /*printf("%d\n",nn);*/
    
    plhs[0]=mxCreateDoubleMatrix(N*q,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(N*N*q*q,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(nn,4,mxREAL);
    

    Hs=mxGetPr(plhs[0]);
    Cs=mxGetPr(plhs[1]);
    me=mxGetPr(plhs[2]);

    exhaust(var,nn,betas,me,me+nn,me+2*nn,me+3*nn,offset,rawprobafile);

    for(i=0;i<N;i++) for(g=0;g<q;g++) Hs[i*q+g]=var->H[i][g];
    for(i=0;i<N;i++) for(j=0;j<N;j++) for(g=0;g<q;g++) for(h=0;h<q;h++) Cs[i*N*q*q+j*q*q+g*q+h]=var->C[i][j][g][h];

    free_var(var);
   

    if(nobeta) free(betas);
    if(noffset) free(offset);

    return;    
}



