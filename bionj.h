#ifndef BIONJ_H
#define BIONJ_H
/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         BIONJ program                                     ;
;       was obtained from http://www.lirmm.fr/~w3ifa/MAAS/BIONJ/BIONJ.html  ;                                    ;
;                                                                           ;
;                         Olivier Gascuel                                   ;
;                                                                           ;
;                         GERAD - Montreal- Canada                          ;
;                         olivierg@crt.umontreal.ca                         ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         gascuel@lirmm.fr                                  ;
;                                                                           ;
;                         UNIX version, written in C                        ;
;                         by Hoa Sien Cuong (Univ. Montreal)                ; 
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

#include <bits/stdc++.h>
using namespace std;

#define PREC 8                             /* precision of branch-lengths  */
#define PRC  100
#define LEN  1000                            /* length of taxon names        */

class BioNj { // stupid code
typedef struct word
{
  char name[LEN];
  struct word *suiv;
}WORD;

typedef struct pointers
{
  WORD *head;
  WORD *tail;
}POINTERS;

/*
void   Initialize(float **delta, FILE *input, int n, POINTERS *trees);

void   Compute_sums_Sx(float **delta, int n);

void   Best_pair(float **delta, int r, int *a, int *b, int n);

void   Finish(float **delta, int n, POINTERS *trees, FILE *output);

void   Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post);

void   Print_output(int i, POINTERS *trees, FILE *output);

float Distance(int i, int j, float **delta);

float Variance(int i, int j, float **delta);

float Sum_S(int i, float **delta);

float Agglomerative_criterion(int i, int j, float **delta, int r);

float Branch_length(int a, int b, float **delta, int r);

float Reduction4(int a, float la, int b, float lb, int i, float lamda,
		 float **delta);

float Reduction10(int a, int b, int i, float lamda, float vab, float
		  **delta);
float Lamda(int a, int b, float vab, float **delta, int n, int r);

float Finish_branch_length(int i, int j, int k, float **delta);

int    Emptied(int i, float **delta);

int    Symmetrize(float **delta, int n);

*/
/*;;;;;;;;;;;  INPUT, OUTPUT, INITIALIZATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                                                                           ;
;              The delta matrix is read from the input-file.                ;
;              It is recommended to put it and the executable in            ;
;              a special directory. The input-file and output-file          ;
;              can be given as arguments to the executable by               ;
;              typing them after the executable (Bionj input-file           ;
;              output-file) or by typing them when asked by the             ;
;              program. The input-file has to be formated according         ;
;              the PHYLIP standard. The output file is formated             ;
;              according to the NEWSWICK standard.                          ;
;                                                                           ;
;              The lower-half of the delta matrix is occupied by            ;
;              dissimilarities. The upper-half of the matrix is             ;
;              occupied by variances. The first column                      ;
;              is initialized as 0; during the algorithm some               ;
;              indices are no more used, and the corresponding              ;
;              positions in the first column are set to 1.                  ;
;                                                                           ;
;              This delta matix is made symmetrical using the rule:         ;
;              Dij = Dji <- (Dij + Dji)/2. The diagonal is set to 0;        ;
;              during the further steps of the algorithm, it is used        ;
;              to store the sums Sx.                                        ;
;                                                                           ;
;              A second array, trees, is used to store taxon names.         ;
;              During the further steps of the algoritm, some               ;
;              positions in this array are emptied while the others         ;
;              are used to store subtrees.                                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


/*;;;;;;;;;;;;;;;;;;;;;;;;;; Initialize        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function reads an input file and return the            ;
;               delta matrix and trees: the list of taxa.                   ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              FILE *input    : pointer to input file                       ;
;              int n          : number of taxa                              ;
;              char **trees   : list of taxa                                ;
;                                                                           ;
; return value:                                                             ;
;              float **delta : delta matrix                                 ;
;              char *trees    : list of taxa                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Initialize(float **delta, FILE *input, int n, POINTERS *trees)
{
  int lig;                                          /* matrix line       */
  int col;                                          /* matrix column     */
  float distance;
  char name_taxon[LEN];                             /* taxon�s name      */
  WORD *name;

  for(lig=1; lig <= n; lig++)
    {
      fscanf(input,"%s",name_taxon);                  /* read taxon�s name */
      name=(WORD *)calloc(1,sizeof(WORD));            /* taxon�s name is   */
      if(name == NULL)                                /* put in trees      */
	{
	  printf("Out of memories !!");
	  exit(0);
	}
      else
	{
	  strcpy(name->name,name_taxon);
	  name->suiv=NULL;
	  trees[lig].head=name;
	  trees[lig].tail=name;
	  for(col= 1; col <= n; col++)
	    {
	      fscanf(input,"%f",&distance);             /* read the distance  */
	      delta[lig][col]=distance;
	    }
	}
    }
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Print_output;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function prints out the tree in the output file.       ;
;                                                                           ;
; input       :                                                             ;
;              POINTERS *trees : pointer to the subtrees.                   ;
;              int i          : indicate the subtree i to be printed.       ;
:              FILE *output   : pointer to the output file.                 ;
;                                                                           ;
; return value: The phylogenetic tree in the output file.                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


void Print_output(int i, POINTERS *trees, FILE *output)
{
  WORD *parcour;
  parcour=trees[i].head;
  while(parcour != NULL)
    {
      fprintf(output,"%s",parcour->name);
      parcour=parcour->suiv;
    }

}



/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                             Utilities                                     ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/



/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Symmetrize  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifies if the delta matrix is symmetric;    ;
;               if not the matrix is made symmetric.                        ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              int n          : number of taxa                              ;
;                                                                           ;
; return value:                                                             ;
;              int symmetric  : indicate if the matrix has been made        ;
;                               symmetric or not                            ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Symmetrize(float **delta, int n)
{
  int lig;                                         /* matrix line        */
  int col;                                         /* matrix column      */
  float value;                                     /* symmetrized value  */
  int symmetric;

  symmetric=1;
  for(lig=1; lig  <=  n; lig++)
    {
      for(col=1; col< lig; col++)
	{
	  if(delta[lig][col] != delta[col][lig])
	    {
	      value= (delta[lig][col]+delta[col][lig])/2;
	      delta[lig][col]=value;
	      delta[col][lig]=value;
	      symmetric=0;
	    }
        }
    }
  if(!symmetric)
    printf("The matrix is not symmetric");
  return(symmetric);
}




/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Concatenate ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function concatenates a string to another.             ;
;                                                                           ;
; input       :                                                             ;
;      char *chain1    : the string to be concatenated.                     ;
;      int ind         : indicate the subtree to which concatenate the      ;
;                        string                                             ;
;      POINTERS *trees  : pointer to subtrees.                              ;
;      int post        : position to which concatenate (front (0) or        ;
;                        end (1))                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post)
{
  WORD *bran;

  bran=(WORD *)calloc(1,sizeof(WORD));
  if(bran == NULL)
    {
      printf("Out of memories");
      exit(0);
    }
  else
    {
      strcpy(bran->name,chain1);
      bran->suiv=NULL;
    }
  if(post == 0)
    {
      bran->suiv=trees[ind].head;
      trees[ind].head=bran;
    }
  else
    {
      trees[ind].tail->suiv=bran;
      trees[ind].tail=trees[ind].tail->suiv;
    }
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Distance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve ant return de distance between taxa  ;
;               i and j from the delta matrix.                              ;
;                                                                           ;
; input       :                                                             ;
;               int i          : taxon i                                    ;
;               int j          : taxon j                                    ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               float distance : dissimilarity between the two taxa         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Distance(int i, int j, float **delta)
{
  if(i > j)
    return(delta[i][j]);
  else
    return(delta[j][i]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Variance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve and return the variance of the       ;
;               distance between i and j, from the delta matrix.            ;
;                                                                           ;
; input       :                                                             ;
;               int i           : taxon i                                   ;
;               int j           : taxon j                                   ;
;               float **delta  : the delta matrix                           ;
;                                                                           ;
; return value:                                                             ;
;               float distance : the variance of  Dij                       ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Variance(int i, int j, float **delta)
{
  if(i > j)
    return(delta[j][i]);
  else
    return(delta[i][j]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Emptied ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifie if a line is emptied or not.          ;
;                                                                           ;
; input       :                                                             ;
;               int i          : subtree (or line) i                        ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               0              : if not emptied.                            ;
;               1              : if emptied.                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Emptied(int i, float **delta)      /* test if the ith line is emptied */
{
  return((int)delta[i][0]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Sum_S;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function retrieves the sum Sx from the diagonal       ;
;                of the delta matrix.                                       ;
;                                                                           ;
;  input       :                                                            ;
;               int i          : subtree i                                  ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
;  return value:                                                            ;
;                float delta[i][i] : sum Si                                 ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Sum_S(int i, float **delta)          /* get sum Si form the diagonal */
{
  return(delta[i][i]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;Compute_sums_Sx;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function computes the sums Sx and store them in the    ;
;               diagonal the delta matrix.                                  ;
;                                                                           ;
; input       :                                                             ;
;     	         float **delta : the delta matrix.                      ;
;     	         int n          : the number of taxa                    ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Compute_sums_Sx(float **delta, int n)
{
  float sum;
  sum = 0.0;
  int i;
  int j;

  for(i= 1; i <= n ; i++)
    {
      if(!Emptied(i,delta))
	{
	  sum=0;
	  for(j=1; j <=n; j++)
	    {
	      if(i != j && !Emptied(j,delta))           /* compute the sum Si */
		sum=sum + Distance(i,j,delta);
	    }
	}
      delta[i][i]=sum;                           /* store the sum Si in */
    }                                               /* delta�s diagonal    */
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Best_pair;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function finds the best pair to be agglomerated by    ;
;                minimizing the agglomerative criterion (1).                ;
;                                                                           ;
;  input       :                                                            ;
;                float **delta : the delta matrix                           ;
;                int r          : number of subtrees                        ;
;                int *a         : contain the first taxon of the pair       ;
;                int *b         : contain the second taxon of the pair      ;
;                int n          : number of taxa                            ;
;                                                                           ;
;  return value:                                                            ;
;                int *a         : the first taxon of the pair               ;
;                int *b         : the second taxon of the pair              ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Best_pair(float **delta, int r, int *a, int *b, int n)
{
  float Qxy;                         /* value of the criterion calculated*/
  int x,y;                           /* the pair which is tested         */
  float Qmin;                        /* current minimun of the criterion */

  Qmin=1.0e300;
  for(x=1; x <= n; x++)
    {
      if(!Emptied(x,delta))
        {
	  for(y=1; y < x; y++)
	    {
	      if(!Emptied(y,delta))
		{
		  Qxy=Agglomerative_criterion(x,y,delta,r);
		  if(Qxy < Qmin-0.000001)
		    {
		      Qmin=Qxy;
		      *a=x;
		      *b=y;
        }
		}
	    }
        }
    }
    //cout << Qmin << '\n';
}


/*;;;;;;;;;;;;;;;;;;;;;;Finish_branch_length;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description :  Compute the length of the branch attached                 ;
;                 to the subtree i, during the final step                   ;
;                                                                           ;
;  input       :                                                            ;
;                int i          : position of subtree i                     ;
;                int j          : position of subtree j                     ;
;                int k          : position of subtree k                     ;
;                float **delta :                                            ;
;                                                                           ;
;  return value:                                                            ;
;                float length  : The length of the branch                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Finish_branch_length(int i, int j, int k, float **delta)
{
  float length;
  length=0.5*(Distance(i,j,delta) + Distance(i,k,delta)
	      -Distance(j,k,delta));
  return(length);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Finish;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function compute the length of the lasts three        ;
;                subtrees and write the tree in the output file.            ;
;                                                                           ;
;  input       :                                                            ;
;                float **delta  : the delta matrix                          ;
;                int n           : the number of taxa                       ;
;                WORD *trees   : list of subtrees                           ;
;                                                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Finish(float **delta, int n, POINTERS *trees, FILE *output)
{
  int l=1;
  int i=0;
  float length;
  char *str;
  WORD *bidon;
  WORD *ele;
  int last[3];                            /* the last three subtrees     */

  str=(char *)calloc(LEN,sizeof(char));

  if(str == NULL)
    {
      printf("Out of memories !!");
      exit(0);
    }
  while(l <= n)
    {                                       /* find the last tree subtree  */
      if(!Emptied(l, delta))
	{
	  last[i]=l;
	  i++;
	}
      l++;
    }

  length=Finish_branch_length(last[0],last[1],last[2],delta);
  fprintf(output,"(");
  Print_output(last[0],trees,output);
  fprintf(output,":");
/*   gcvt(length,PREC, str); */
/*   fprintf(output,"%s,",str); */
  fprintf(output,"%f,",length);

  length=Finish_branch_length(last[1],last[0],last[2],delta);
  Print_output(last[1],trees,output);
  fprintf(output,":");
/*   gcvt(length,PREC, str); */
/*   fprintf(output,"%s,",str); */
  fprintf(output,"%f,",length);

  length=Finish_branch_length(last[2],last[1],last[0],delta);
  Print_output(last[2],trees,output);
  fprintf(output,":");
/*   gcvt(length,PREC,str); */
/*   fprintf(output,"%s",str); */
  fprintf(output,"%f",length);
  fprintf(output,");");
  fprintf(output,"\n");

  for(i=0; i < 3; i++)
    {
      bidon=trees[last[i]].head;
      ele=bidon;
      while(bidon!=NULL)
	{
	  ele=ele->suiv;
	  free(bidon);
	  bidon=ele;
	}
    }
  free(str);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                          Formulae                                         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


float Agglomerative_criterion(int i, int j, float **delta, int r)
{
  float Qij;
  Qij=(r-2)*Distance(i,j,delta)                           /* Formula (1) */
    -Sum_S(i,delta)
    -Sum_S(j,delta);

  return(Qij);
}


float Branch_length(int a, int b, float **delta, int r)
{
  float length;
  length=0.5*(Distance(a,b,delta)                         /* Formula (2) */
	      +(Sum_S(a,delta)
		-Sum_S(b,delta))/(r-2));
  return(length);
}


float Reduction4(int a, float la, int b, float lb, int i, float lamda,
		 float **delta)
{
  float Dui;
  Dui=lamda*(Distance(a,i,delta)-la)
    +(1-lamda)*(Distance(b,i,delta)-lb);                /* Formula (4) */
  return(Dui);
}


float Reduction10(int a, int b, int i, float lamda, float vab,
		  float **delta)
{
  float Vci;
  Vci=lamda*Variance(a,i,delta)+(1-lamda)*Variance(b,i,delta)
    -lamda*(1-lamda)*vab;                              /*Formula (10)  */
  return(Vci);
}

float Lamda(int a, int b, float vab, float **delta, int n, int r)
{
  float lamda=0.0;
  int i;

  if(vab==0.0)
    lamda=0.5;
  else
    {
      for(i=1; i <= n ; i++)
	{
          if(a != i && b != i && !Emptied(i,delta))
            lamda=lamda + (Variance(b,i,delta) - Variance(a,i,delta));
	}
      lamda=0.5 + lamda/(2*(r-2)*vab);
    }                                              /* Formula (9) and the  */
  if(lamda > 1.0)                                /* constraint that lamda*/
    lamda = 1.0;                             /* belongs to [0,1]     */
  if(lamda < 0.0)
    lamda=0.0;
  return(lamda);
}
/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         Main program                                      ;
;                                                                           ;
;                         argc is the number of arguments                   ;
;                         **argv contains the arguments:                    ;
;                         the first argument has to be BIONJ;               ;
;                         the second is the inptu-file;                     ;
;                         the third is the output-file.                     ;
;                         When the input and output files are               ;
;                         not given, the user is asked for them.            ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

public :
int create(const char *inputFile, const char *outputFile) {

  FILE *input;                            /* pointer to input file       */
  FILE *output;                           /* pointer to output file      */
  POINTERS *trees;                        /* list of subtrees            */
  char *Name_fich1;                       /* name of the input file      */
  char *Name_fich2;                       /* name of the output file     */
  char *chain1;                           /* stringized branch-length    */
  char *chain2;                           /* idem                        */
  int *a, *b;                             /* pair to be agglomerated     */
  float **delta;                          /* delta matrix                */
  float la;                               /* first taxon�s branch-length */
  float lb;                               /* second taxon�s branch-length*/
  float vab;                              /* variance of Dab             */
  float lamda;
  int i;
  int ok;
  int r;                                  /* number of subtrees          */
  int n;                                  /* number of taxa              */
  int x, y;
  //float t;


  /*   Allocation of memories    */

  Name_fich1=(char*)calloc(LEN,sizeof(char));
  Name_fich2=(char*)calloc(LEN,sizeof(char));
  a=(int*)calloc(1,sizeof(int));
  b=(int*)calloc(1,sizeof(int));
  chain1=(char *)calloc(LEN,sizeof(char));
  chain2=(char *)calloc(LEN,sizeof(char));

  input= fopen(inputFile,"r");
  fscanf(input,"%d",&n);

  output= fopen(outputFile,"w");
  /*      Create the delta matrix     */

  delta=(float **)calloc(n+1,sizeof(float*));
  for(i=1; i<= n; i++)
    {
      delta[i]=(float *)calloc(n+1, sizeof(float));
      if(delta[i] == NULL)
	{
	  printf("Out of memories!!");
	  exit(0);
	}
    }
  trees=(POINTERS *)calloc(n+1,sizeof(POINTERS));
  if(trees == NULL)
    {
      printf("Out of memories!!");
      exit(0);
    }
  /*   initialise and symmetrize the running delta matrix    */

  rewind(input);
  while(fscanf(input,"%d",&n) != EOF )
    {
      r=n;
      *a=0;
      *b=0;
      Initialize(delta, input, n, trees);
      ok=Symmetrize(delta, n);
      if(!ok)
	printf("\n The matrix  is not symmetric.\n ");
      while (r > 3)                             /* until r=3                 */
	{
	  Compute_sums_Sx(delta, n);             /* compute the sum Sx       */
	  Best_pair(delta, r, a, b, n);          /* find the best pair by    */
	  vab=Variance(*a, *b, delta);           /* minimizing (1)           */
	  la=Branch_length(*a, *b, delta, r);    /* compute branch-lengths   */
	  lb=Branch_length(*b, *a, delta, r);    /* using formula (2)        */
	  lamda=Lamda(*a, *b, vab, delta, n, r); /* compute lambda* using (9)*/
	  for(i=1; i <= n; i++)
	    {
	      if(!Emptied(i,delta) && (i != *a) && (i != *b))
		{
		  if(*a > i)
		    {
		      x=*a;
		      y=i;
		    }
		  else
		    {
		      x=i;
		      y=*a;                           /* apply reduction formulae */
		    }                                  /* 4 and 10 to delta        */
		  delta[x][y]=Reduction4(*a, la, *b, lb, i, lamda, delta);
		  delta[y][x]=Reduction10(*a, *b, i, lamda, vab, delta);
		}
	    }
	  strcpy(chain1,"");                     /* agglomerate the subtrees */
	  strcat(chain1,"(");                    /* a and b together with the*/
	  Concatenate(chain1, *a, trees, 0);     /* branch-lengths according */
	  strcpy(chain1,"");                     /* to the NEWSWICK format   */
	  strcat(chain1,":");

	  sprintf(chain1+strlen(chain1),"%f",la);
/* 	  gcvt(la,PREC, chain2); */
/* 	  strcat(chain1, chain2); */

	  strcat(chain1,",");
	  Concatenate(chain1,*a, trees, 1);
	  trees[*a].tail->suiv=trees[*b].head;
	  trees[*a].tail=trees[*b].tail;
	  strcpy(chain1,"");
	  strcat(chain1,":");

	  sprintf(chain1+strlen(chain1),"%f",lb);
/* 	  gcvt(lb, PREC, chain2); */
/* 	  strcat(chain1, chain2); */
	  strcat(chain1,")");
	  Concatenate(chain1, *a, trees, 1);
	  delta[*b][0]=1.0;                     /* make the b line empty     */
	  trees[*b].head=NULL;
	  trees[*b].tail=NULL;
	  r=r-1;                                /* decrease r                */
	}
      Finish(delta, n, trees, output);       /* compute the branch-lengths*/
      for(i=1; i<=n; i++)       	          /* of the last three subtrees*/
	{				                /* and print the tree in the */
	  delta[i][0]=0.0;		          /* output-file               */
	  trees[i].head=NULL;
	  trees[i].tail=NULL;
	}
    }
  free(trees);
  for(i=n; i>=1; i--)
    {
      free(delta[i]);
    }
  free(delta);
  /* Minh free memory-leak */
  free(chain2);
  free(chain1);
  free(b);
  free(a);
  free(Name_fich2);
  free(Name_fich1);
  /* Minh done */
  fclose(input);
  fclose(output);

  return 0;
}
};

class RapidNJ {
  FILE *input, *output;  
  struct TREE {
    vector<vector<pair<int, float>>> child;
    vector<string> nameTax;
    int cur = 0;
    FILE *output;
    TREE() {}
    void initialize(vector<string> nameTax, FILE *output) {
      int n = nameTax.size();
      cur = n;
      child.resize(n * 2 + 2);
      this->nameTax = nameTax; 
      this->output = output;
    }
    int createNode() {
      return cur++;
    }
    void addEdge(int u, int v, float w) { // u is v's parent
      child[u].push_back({v, w});
    }
    vector<pair<int, float>> calcTotalDistanceToLeaves(int top, float branchLen) {
      vector<pair<int, float>> ans;
      queue<pair<int, float>> q;
      q.push({top, branchLen});
      while (q.size()) {
        int u = q.front().first; 
        float d = q.front().second;
        q.pop();
        if (child[u].empty()) {
          ans.push_back({u, d});
          continue;
        }
        for (auto [v, w]: child[u]) {
          q.push({v, d + w});
        }
      }
      return ans;
    }
    void print(int u, float l) {
      if (child[u].empty()) {
        fprintf(output, "%s:%f", nameTax[u].c_str(), l);
      } else {
        fprintf(output, "(");
        for (int i = 0; i < child[u].size(); ++i) {
          int v = child[u][i].first;
          float w = child[u][i].second;
          print(v, w);
          if (i != child[u].size() - 1) fprintf(output, ",");
        }
        fprintf(output, "):%f", l);
      }
    }
    void print() {
      int u = cur - 1;
      fprintf(output, "(");
      for (int i = 0; i < child[u].size(); ++i) {
        int v = child[u][i].first;
        float w = child[u][i].second;
        print(v, w);
        if (i != child[u].size() - 1) fprintf(output, ",");
      }
      fprintf(output, ");");
    }
  } tree;
  
                        
  float **delta, **dist;                        
  int r, n, time = 0;
  bool *garbage;
  pair<float, int> **sorted;
  int *len, *rank, *realNode;
  
  void Initialize() { 
    float distance;
    char cur[LEN];
    vector<string> name(n);    
    for (int i = 0; i < n; ++i) {
      fscanf(input, "%s", cur);
      name[i] = cur;
      for (int j = 0; j < n; ++j) {
        fscanf(input, "%f", &distance);
        delta[i][j] = dist[i][j] = distance;
      }
    }
    tree.initialize(name, output);
  }

  int Symmetrize() {
    int lig;                                  
    int col;                                 
    float value;                              
    int symmetric = 1;
    for (lig = 0; lig < n; lig++) {
      for (col = 0; col < lig; col++) {
        if (delta[lig][col] != delta[col][lig]) {
          value = (delta[lig][col] + delta[col][lig]) / 2;
          delta[lig][col] = value;
          delta[col][lig] = value;
          symmetric = 0;
        }
      }
    }
    if(!symmetric)
        printf("The matrix is not symmetric");
    return(symmetric);
  }
  float Distance(int i, int j) {
    if (i > j) return (delta[i][j]);
    else return(delta[j][i]);
  }

  float Variance(int i, int j) {
    if (i > j) return (delta[j][i]);
    else return (delta[i][j]);
  }

  int Emptied(int i) {
    return (realNode[i] == -1);
  }

  float Sum_S(int i) {
    return(delta[i][i]);
  }

  void Compute_sums_Sx() {
    float sum = 0;
    for(int i = 0; i < n; i++) {
      if(!Emptied(i)) {
        sum = 0;
        for(int j = 0; j < n; j++) {
          if(i != j && !Emptied(j))          
          sum += Distance(i, j);
        }
      }
      delta[i][i] = sum;                          
    }                                              
  }
  int last = 1;
  void Best_pair(int *a, int *b) {
    float Qmin = 1e300;
    float Umax = -1e300;
    for (int i = 0; i < n; ++i) 
        if (!Emptied(i))
            Umax = max(Umax, Sum_S(i));  
    auto findMinRow = [&](int x) { 
      int dead = 0;
      for (int i = 0; i < len[x]; ++i) {
        int y = sorted[x][i].second;
        if (Emptied(y) || rank[x] < rank[y]) {
          ++dead;
        } else {
          if (sorted[x][i].first * (r - 2) - Sum_S(x) - Umax > Qmin) break;
          float Qxy = Agglomerative_criterion(x, y);
          if (Qxy < Qmin - 0.000001) {
              Qmin = Qxy;
              *a = x;
              *b = y;
          } 
        }
      }
      return dead;
    };
    auto findMinRowGarbage = [&](int x) { 
      int cur = 0;
      for (int i = 0; i < len[x]; ++i) {
        int y = sorted[x][i].second;
        if (Emptied(y) || rank[x] < rank[y]) continue;
        float Qxy = Agglomerative_criterion(x, y);
        if (Qxy < Qmin - 0.000001) {
          Qmin = Qxy;
          *a = x;
          *b = y;
        }
        sorted[x][cur++] = sorted[x][i]; 
      }
      len[x] = cur;
    };
    int x = last;
    do {
      if(!Emptied(x)) {
        if (garbage[x]) {
          findMinRowGarbage(x);
          garbage[x] = 0;
        } else {
          int dead = findMinRow(x);
          if (dead > r / 2) garbage[x] = 1;
        }
      }
      x = (x == n - 1 ? 0 : x + 1);
    } while (x != last);
  }

  float Finish_branch_length(int i, int j, int k) {
    return 0.5 * (Distance(i, j) + Distance(i, k) - Distance(j, k));
  }

  void Finish() {
    int l = 0;
    int i = 0;
    float length;
    int last[3];                  
    while (l < n) {                            
        if (!Emptied(l)) {
            last[i] = l;
            i++;
        }
        l++;
    }
    int root = tree.createNode();
    tree.addEdge(root, realNode[last[0]], Finish_branch_length(last[0], last[1], last[2]));
    tree.addEdge(root, realNode[last[1]], Finish_branch_length(last[1], last[0], last[2]));
    tree.addEdge(root, realNode[last[2]], Finish_branch_length(last[2], last[1], last[0]));
    tree.print();
  }

  float Agglomerative_criterion(int i, int j) {
    return (r - 2) * Distance(i, j) - Sum_S(i) - Sum_S(j);
  }

  float Branch_length(int a, int b) {
    return 0.5 * (Distance(a, b) + (Sum_S(a) - Sum_S(b)) / (r - 2));
  }

  float Reduction4(int a, float la, int b, float lb, int i, float lamda) {
    return lamda * (Distance(a, i) - la) + (1 - lamda) * (Distance(b, i) - lb);              
  }

  float Reduction10(int a, int b, int i, float lamda, float vab) {
    return lamda * Variance(a, i) + (1 - lamda) * Variance(b, i) - lamda * (1 - lamda) * vab; 
  }

  float Lamda(int a, int b, float vab) {
    float lamda = 0.0;
    if (vab == 0.0) lamda = 0.5;
    else {
      for (int i = 0; i < n; i++) {
          if (a != i && b != i && !Emptied(i))
              lamda += (Variance(b, i) - Variance(a, i));
      }
      lamda = 0.5 + lamda / (2 * (r - 2) * vab);
    }                                             
    if(lamda > 1.0) lamda = 1.0;                            
    if(lamda < 0.0) lamda = 0.0;
    return(lamda);
  }

public :
  void readInp(const char *inputFile, const char *outputFile) {
    input = fopen(inputFile, "r");
    output = fopen(outputFile, "w");

    fscanf(input, "%d", &n);

    garbage = (bool*) calloc(n, sizeof(bool));
    rank = (int*) calloc(n, sizeof(int));
    len = (int*) calloc(n, sizeof(int));
    delta = (float**) calloc(n, sizeof(float*));
    dist = (float**) calloc(n, sizeof(float*));
    realNode = (int*) calloc(n, sizeof(int));
    for (int i = 0; i < n; i++) {
        delta[i] = (float*) calloc(n, sizeof(float));
        dist[i] = (float*) calloc(n, sizeof(float));
        if(delta[i] == NULL) {
            printf("Out of memories!!");
            exit(0);
        }
    }
    r = n;
    Initialize();
    if(!Symmetrize()) printf("\n The matrix  is not symmetric.\n ");
  }
  int create(const char *inputFile, const char *outputFile) {
    int* a = (int*)calloc(1, sizeof(int));
    int* b = (int*)calloc(1, sizeof(int));

    *a = 0;
    *b = 0;

    for (int i = 0; i < n; ++i) rank[i] = 0, realNode[i] = i;
    sorted = new pair<float, int>* [n];
    for (int i = 0; i < n; ++i) {
        sorted[i] = new pair<float, int> [n];
        for (int j = 0; j < i; ++j) {
            sorted[i][j].first = Distance(i, j);
            sorted[i][j].second = j; 
        }
        len[i] = i;
        sort(sorted[i], sorted[i] + i);
    }
    float *tmp = new float[n];
    Compute_sums_Sx();
    while (r > 3) {
      Best_pair(a, b);       
      float vab = Variance(*a, *b);    
      float la = Branch_length(*a, *b); 
      float lb = Branch_length(*b, *a);  
      float lamda = Lamda(*a, *b, vab); 
      
      rank[*a] = ++time;
      len[*a] = 0;
      float sum = 0;
      int x, y;
      for (int i = 0; i < n; ++i) tmp[i] = 0;
      for (int i = 0; i < n; i++) {
        if (!Emptied(i) && (i != *a) && (i != *b)) {
          if (*a > i) {
              x = *a;
              y = i;
          } else {
              x = i;
              y = *a;                     
          }                                
          tmp[i] -= Distance(*a, i) + Distance(*b, i);                        
          delta[x][y] = Reduction4(*a, la, *b, lb, i, lamda);
          tmp[i] += Distance(*a, i);

          sum += Distance(*a, i);
          delta[y][x] = Reduction10(*a, *b, i, lamda, vab);
          sorted[*a][len[*a]++] = {delta[x][y], i};
        }
      }
      for (int i = 0; i < n; ++i) delta[i][i] += tmp[i];
      sort(sorted[*a], sorted[*a] + len[*a]);
      delta[*a][*a] = sum;          
      last = *a;

      int par = tree.createNode();
      tree.addEdge(par, realNode[*a], la);
      tree.addEdge(par, realNode[*b], lb);

      realNode[*b] = -1;
      realNode[*a] = par;
      --r;
    }
    Finish();   
    for(int i = n - 1; i >= 0; i--) free(delta[i]);
    free(realNode);
    free(delta);
    free(b);
    free(a);
    fclose(input);
    fclose(output);
    return 0;
  }
  void calcError() {
    float sumSquares = 0;
    for (int u = n; u < tree.cur; ++u) {
      vector<vector<pair<int, float>>> totalDist;
      for (auto [v, w]: tree.child[u]) 
        totalDist.push_back(tree.calcTotalDistanceToLeaves(v, w));
      for (int i = 0; i < totalDist.size(); ++i)
        for (int j = i + 1; j < totalDist.size(); ++j) {
          for (auto [c1, d1]: totalDist[i])
            for (auto [c2, d2]: totalDist[j]) {
              float diff = (d1 + d2 - dist[c1][c2]);
              sumSquares += diff * diff;
            }
        }
    }
    float rms = sqrt(sumSquares * 2.0 / n / (n - 1.0));
    cout << "Root Mean Square Error: " << rms << '\n';
    free(dist);
  }
};

#endif
