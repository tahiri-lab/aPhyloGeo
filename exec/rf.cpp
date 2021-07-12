/*
	int a = 0;
* --------------------------------------------*
* Computation of the Robinson and Foulds      *
* topological distance between two (or more)  *
* trees given their distance matrices         *
*                                             *
* Input data file format:                     *
* 1. Number of objects (matrix size) n        *
* 2. A sequence of m squared additive         *
*    distance matrices of size n by n         *
*    assosiated with trees                    *
*                                             *
* Output data:                                *
* 1. A sequence of m-1 integer numbers        *
*    containing the Robinson and Foulds       *
*    distances between the tree associated    *
*    with the first matrix and all the other  *
*    trees associated with the following ones *
*                                             *
* --------------------------------------------*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#pragma warning(disable:4996)

#define INFINI 999999.99
// 1024*1024*500 element
#define SEUIL 104857600

//unsigned int n;
FILE* Output1;
int a;

/* Function's Propotypes */

void odp(double**, int*, int*, int*, int n);
int Bipartition_Table(double**, int**, int*, int n);
int Table_Comparaison(int**, int**, int*, int*, int, int, int n);


int nbSpeciesPhylip(const char* fichier) {
	FILE* in;
	int nbItem;

	in = fopen(fichier, "r+");
	fscanf(in, "%d", &nbItem);

	fclose(in);
	a = 0;
	if ((nbItem > 10000) || (nbItem < 3))
		exit(8);
	return nbItem;
}

int fileInNewickFormat(char* fichier) {

	char car;
	//	int retour;
	int i = 0;
	FILE* in = fopen(fichier, "r+");
	//= rechercher le premier caractere
	do {
		fscanf(in, "%c", &car);
		i++;
		if (i > SEUIL) exit(-1);
	} while ((car == ' ') || (car == '\t') || (car == '\n') || (car == '\r'));


	fclose(in);
	a = 0;
	if (car == '(') return 1;
	else return 0;

}
void filtrerMatrice2(double** dissSpecies, double** dissGene, char** nomsSpecies, char** nomsGene, int nbSpecies, int nbGene, char* fichier) {

	int i, j, temoin;
	/*FILE * out;

	if((out = fopen(fichier,"a+")) == NULL){
		printf("Can't open %s",fichier);
		exit(-1);
	}

	fprintf(out,"\n=>Filtre");

	*/
	for (i = 1; i <= nbSpecies; i++) {
		//printf("\n%s -",nomsSpecies[i]);
		temoin = 0;
		for (j = 1; j <= nbGene; j++) {
			//printf(" %s",nomsGene[j]);
			if (strcmp(nomsSpecies[i], nomsGene[j]) == 0)
				temoin = 1;
		}
		if (temoin == 0) {
			//fprintf(out,"\nSpecies:%s",nomsSpecies[i]);
			for (j = 1; j <= nbSpecies; j++) {
				dissSpecies[i][j] = dissSpecies[j][i] = -1;
			}
			strcpy(nomsSpecies[i], "");
		}
	}
	//printf("==============");
	for (i = 1; i <= nbGene; i++) {
		//printf("\n%s -",nomsGene[i]);
		temoin = 0;
		for (j = 1; j <= nbSpecies; j++) {
			//	printf(" %s",nomsSpecies[j]); 	
				//if(strlen(nomsSpecies[j]) > 1) 
			if (strcmp(nomsSpecies[j], nomsGene[i]) == 0)
				temoin = 1;
		}
		if (temoin == 0) {
			//fprintf(out,"\nGene:%s",nomsGene[i]);
			for (j = 1; j <= nbGene; j++) {
				dissGene[i][j] = dissGene[j][i] = -1;
			}
			strcpy(nomsGene[i], "");
		}
	}

	//fclose(out);
}

void filtrerMatrice(double** dissSpecies, double** dissGene, char** nomsSpecies, char** nomsGene, int nbSpecies, int nbGene) {

	int i, j, temoin;

	//printf("\nn=%d",nbSpecies);
	for (i = 1; i <= nbSpecies; i++) {
		//	printf("\n%s -",nomsSpecies[i-1]);
		temoin = 0;
		for (j = 1; j <= nbGene; j++) {
			//printf(" %s",nomsGene[j-1]);
			if (strcmp(nomsSpecies[i - 1], nomsGene[j - 1]) == 0)
				temoin = 1;
		}
		if (temoin == 0) {
			//		printf("*");
			for (j = 1; j <= nbSpecies; j++) {
				//dissSpecies[i+1][j+1] = dissSpecies[j+1][i+1] = -1;
				dissSpecies[i][j] = dissSpecies[j][i] = -1;
			}
			strcpy(nomsSpecies[i - 1], "");
		}
	}

	//	printf("\n\nn=%d",nbGene);
	for (i = 1; i <= nbGene; i++) {
		//		printf("\n%s -",nomsGene[i-1]);
		temoin = 0;
		for (j = 1; j <= nbSpecies; j++) {
			//	printf(" %s",nomsSpecies[j]);
			//	if(strlen(nomsSpecies[j]) > 1) 
			if (strcmp(nomsSpecies[j - 1], nomsGene[i - 1]) == 0)
				temoin = 1;
		}
		if (temoin == 0) {
			//		printf("*");
			for (j = 1; j <= nbGene; j++) {
				//dissGene[i+1][j+1] = dissGene[j+1][i+1] = -1;
				dissGene[i][j] = dissGene[j][i] = -1;
			}
			strcpy(nomsGene[i - 1], "");
		}
	}
}

int ecrireMatrice2(double** source, double** dest, int taille, char** nomsSource, char** nomsDest) {
	int i, ii, j, jj, finalTaille;
	//	FILE *out;

	finalTaille = 0;

	for (i = 0; i < taille; i++)
		if (strcmp(nomsSource[i], "") != 0)
			finalTaille = finalTaille + 1;
	if (finalTaille < 3) {
		printf("There are less than 3 same species !");
		exit(5);
	}
	ii = 1;
	for (i = 0; i < taille; i++) {
		if (strcmp(nomsSource[i], "") != 0) {//if(strlen(noms[i]) > 1){
			strcpy(nomsDest[ii - 1], nomsSource[i]);
			jj = 1;
			for (j = 0; j < taille; j++)
				if (source[i + 1][j + 1] != -1) {
					dest[ii][jj] = source[i + 1][j + 1];
					jj++;
				}
			ii++;
		}
	}
	return finalTaille;
}

//=========================================================
//
//=========================================================
double** createDoubleMatrix(int size) {
	double** matrice;
	int i;

	matrice = (double**)malloc((size) * sizeof(double*));

	for (i = 0; i < size; i++)
		matrice[i] = (double*)malloc(size * sizeof(double));
	return matrice;
}

//==========================================================
//
//==========================================================
void freeDoubleMatrix(double** matrice, int size) {
	int i;
	for (i = 0; i < size; i++)
		free(matrice[i]);
	//free(matrice);	
}


char** createStringMatrix(int nbString, int stringSize) {
	char** noms;
	int i;

	noms = (char**)malloc(nbString * sizeof(char*));

	for (i = 0; i < nbString; i++)
		noms[i] = (char*)malloc(stringSize);

	return noms;
}


int nbSpeciesNewick(char* fichier) {
	FILE* data;
	int i = 0;
	int n = 0;
	char symbol;
	char symbolOld = ' ';
	int temoin = 0;
	if ((data = fopen(fichier, "r")) == 0) { printf("\n%s:Open Failed....", fichier); return -1; }

	//Correctness of the Newick format verification

	do
	{
		symbol = getc(data);
		if (symbol == ':') i++;
		if (symbol == ':' && symbolOld != ')' && temoin != 1) n++;
		if (symbol >= '0' && symbol <= '9' && symbolOld == ')') temoin = 1;
		if (symbol == ':' && temoin == 1) temoin = 0;
		symbolOld = symbol;
	} while (symbol != ';');

	/*
	while ((symbol=getc(data))!=';')
	{
		if (symbol==':') i++;
		if (symbol == ':' && symbolOld != ')') n++;
		symbolOld = symbol;
	}
	*/
	fclose(data);
	return n;
}


void TreeMatrix(double** DISS, long int* ARETE, double* LONGUEUR, int n, int kt)
{
	int i, j, k, A, B;
	double** Adjacence;
	double longueur;

	Adjacence = (double**)malloc((2 * n - 1) * sizeof(double*));

	for (i = 0; i <= 2 * n - 2; i++)
	{
		Adjacence[i] = (double*)malloc((2 * n - 1) * sizeof(double));
	}

	for (i = 1; i <= 2 * n - 2; i++)
		for (j = 1; j <= 2 * n - 2; j++)
			Adjacence[i][j] = INFINI;

	for (i = 1; i <= 2 * n - 3 - kt; i++)
	{
		A = ARETE[2 * i - 2];
		B = ARETE[2 * i - 1];
		longueur = LONGUEUR[i - 1];
		Adjacence[A][B] = Adjacence[B][A] = longueur;
	}

	for (i = 1; i <= 2 * n - 2; i++)
		for (j = i + 1; j <= 2 * n - 2; j++)
		{
			DISS[i][j] = DISS[j][i] = Adjacence[i][j];
		}

	for (i = 1; i <= 2 * n - 2; i++)
		DISS[i][i] = 0.0;
	//int bidon;
	for (i = 1; i <= 2 * n - 2; i++)
		for (j = 1; j <= 2 * n - 2; j++)
			for (k = 1; k <= 2 * n - 2; k++)
			{
				if ((DISS[j][i] + DISS[i][k]) < DISS[j][k])
					DISS[j][k] = DISS[k][j] = DISS[j][i] + DISS[i][k];
			}

	//	afficherMatriceDouble(DISS,n);
//			free(Adjacence);
/*	for(i=1;i<=2*n-2;i++)
		for(j=1;j<=2*n-2;j++)
			DISS[i-1][j-1] = DISS[i][j];*/
}

static void xtoa(unsigned long val, char* buf, unsigned radix, int is_neg) {
	char* p;                /* pointer to traverse string */
	char* firstdig;         /* pointer to first digit */
	char temp;              /* temp char */
	unsigned digval;        /* value of digit */

	p = buf;

	if (is_neg) {
		/* negative, so output '-' and negate */
		*p++ = '-';
		val = (unsigned long)(-(long)val);
	}

	firstdig = p;           /* save pointer to first digit */

	do {
		digval = (unsigned)(val % radix);
		val /= radix;       /* get next digit */

		/* convert to ascii and store */
		if (digval > 9)
			*p++ = (char)(digval - 10 + 'a');  /* a letter */
		else
			*p++ = (char)(digval + '0');       /* a digit */
	} while (val > 0);

	/* We now have the digit of the number in the buffer, but in reverse
	order.  Thus we reverse them now. */

	*p-- = '\0';            /* terminate string; p points to last digit */

	do {
		temp = *p;
		*p = *firstdig;
		*firstdig = temp;   /* swap *p and *firstdig */
		--p;
		++firstdig;         /* advance to next two digits */
	} while (firstdig < p); /* repeat until halfway */
}

/* Actual functions just call conversion helper with neg flag set correctly,
and return pointer to buffer. */

char* itoa(int val, char* buf, int radix) {
	if (radix == 10 && val < 0)
		xtoa((unsigned long)val, buf, radix, 1);
	else
		xtoa((unsigned long)(unsigned int)val, buf, radix, 0);
	return buf;
}

//===============================================================================================
//
//===============================================================================================
int lectureNewick(const char* fichier, long int* ARETE, double* LONGUEUR, char** lesNoms, int* kt)
{
	int n;
	int cpt_x;
	// Ce sous programme permet de lire un arbre au format newick et de le transcrire dans
	// des vecteurs arete-longueur commencant à 1
	// ATTENTION: les noms commencent à 0
	// 7 octobre 2004
	// Elmaestro
	//printf("\nlecture Newick : ");
	// TODO: Add your command handler code here
	int FAIL = -1;
	int i, j, j1, k, a, a1, a2, a3, VertexNumber, numero;
	char symbol, * string, * string1, * string2/* *string3,c ,**Et*/;
	int taxaPos; // le nombre de taxas recupéré
	int aretePos; // le nombre d'aretes recupéré
	char symbolOld = ' ';
	int zz, xx, jj, ii;
	double longueur;
	char* tempString;
	int cpt = 0;
	int na;
	char* string4 = (char*)malloc(1000);
	int temoin = 0;
	char* newick = (char*)malloc(10000);

	//long int *ARETE;
	//double *LONGUEUR; 	

	FILE* data;
	if ((data = fopen(fichier, "r")) == 0) { printf("\n%s:Open Failed....", fichier); return FAIL; }
	i = 0;
	while ((symbol = getc(data)) != EOF)
	{
		newick[i] = symbol;
		i++;
	}
	newick[i] = '\0';
	fclose(data);
	//printf(" fclose(data)");
	a = 0;

	//Correctness of the Newick format verification
	i = 0;
	n = 0;

	do
	{
		symbol = newick[cpt++];
		if (symbol == ':') i++;
		if (symbol == ':' && symbolOld != ')' && temoin != 1) n++;
		if (symbol >= '0' && symbol <= '9' && symbolOld == ')') temoin = 1;
		if (symbol == ':' && temoin == 1) temoin = 0;
		symbolOld = symbol;
	} while (symbol != '\0');

	cpt = 0;
	na = i;
	if (i <= 2 * n - 3)(*kt) = i;
	else (*kt) = 2 * n - 3;

	if (i == 0) { printf("Incorrect Newick file format. Edge lengths must be indicated after a ':' characters."); return FAIL; }

	if ((i > 2 * n - 3) || (i < n)) {
		//printf("Unable to read your input data, please check it and try again...");
		//exit(-1);
	}

	j = 0;
	do {
		symbol = newick[cpt++];
		if (symbol == '(') j++;
	} while (symbol != '\0');

	cpt = 0;
	j1 = 0;

	do {
		symbol = newick[cpt++];
		if (symbol == ')') j1++;
	} while (symbol != '\0');

	cpt = 0;

	// verification des arités de l'arbre
	if (j1 != j) { printf("Incorrect Newick file format. Number of right parentheses must be equal to number of left parentheses."); return FAIL; }
	//else if (j!=n-2) { printf("Incorrect Newick file format. Only trees with vertices of degree 1 and 3 are allowed by T-REX."); fclose (data); return FAIL;}

	k = 0;

	do {
		symbol = newick[cpt++];
		if (symbol == ',') k++;
	} while (symbol != '\0');

	cpt = 0;
	//if (k!=(n-1)) { printf("Incorrect Newick file format. Number of objects must be equal to number of commas plus 1."); fclose (data); return FAIL;}

	a = 0;

	do {
		symbol = newick[cpt++];
		if (symbol == ';') a++;
	} while (symbol != '\0');
	cpt = 0;

	if (a == 0) { printf("Incorrect Newick file format. Newick string must be followed by a ';' character."); return FAIL; }
	else if (a > 1) { printf("Incorrect Newick file format. Newick string must contain (in the end) only one ';' character."); return FAIL; }

	a = 0;
	do {
		symbol = newick[cpt++]; a++;
	} while (symbol == ' ');

	cpt = 0;

	if (symbol != '(') { printf("Incorrect Newick file format. Newick string must begin with a '(' character."); return FAIL; }

	a = 0;

	do {
		symbol = newick[cpt++];
		if (symbol == '%') a++;
	} while (symbol != '\0');

	cpt = 0;

	if (a > 0) { printf("Incorrect Newick file format. Newick string cannot contain '%%' character."); return FAIL; }

	do
	{
		symbol = newick[cpt++];
		if ((symbol == '(') || (symbol == ','))
		{
			symbol = newick[cpt++]; a = 0;
			if ((symbol != '(') && (symbol != ',') && (symbol != ';') && (symbol != ':'))
			{
				cpt--;
				do {
					symbol = newick[cpt++]; a++;
				} while ((symbol != ':') && (symbol != '\0'));
			}
			else cpt--;
			if (a > 50) { printf("Incorrect Newick file format. Names of objects must not exceed 50 characters.");  return FAIL; }
		}
	} while (symbol != '\0');
	cpt = 0;

	string = (char*)malloc((100000 * n) * sizeof(char));
	string2 = (char*)malloc((100000 * n) * sizeof(char));
	string1 = (char*)malloc((200000) * sizeof(char));

	if ((string == NULL) || (string1 == NULL) || (string2 == NULL)/*||(string3 == NULL)*/)
	{
		printf("Input data are too large or not a correct Newick file chosen"); return FAIL;
	}

	a = 0;

	do {
		symbol = newick[cpt++];
		if ((symbol != ' ') && (symbol != '\n') && (symbol != '\t')) { string[a++] = symbol; }
	} while (symbol != '\0');

	int boot_value;
	int temoin3 = 0;
	k = 0; VertexNumber = n;
	//a1 = 0;
	//a2 = 0;
	taxaPos = 1;    // nous allons commencer à mettre les taxas à la position 1
	aretePos = 1;
	boot_value = 0;

	while (string[0] == '(')   // traiter toute la chaine
	{
		a1 = 0;
		a2 = 0;
		while (string[a2] != ')')  // traiter la paire () la plus profonde
		{
			if (string[a2] == '(') a1 = a2;  // retrouver ;a parenthèse ouvrante
			a2++;
		}


		// a   => contient la longueur de la chaine
		// a1  => contient le debut d'un noeud à traiter
		// a2  => contient la fin d'un noeud à traiter
		// a3  => délimite le noeud et sa longueur

		zz = a1 + 1;
		VertexNumber++;  // augmenter le nombre de noeuds
		boot_value = 0;
		for (ii = a1 + 1; ii <= a2; ii++)
		{// decortiquer cette chaine

			if (string[ii] == ':')
			{
				xx = 0;
				a3 = ii + 1;

				if (string[zz] == '%')
				{ // cela veut dire que c'est un  noeud que l'on traite

					for (jj = zz + 1; jj < ii; jj++)
					{
						if (string[jj] == '|')
							break;
						string1[xx++] = string[jj];
					}
					temoin3 = 1;
					string1[xx++] = '\0';
					numero = atoi(string1);

					if (string[jj] == '|') {
						boot_value = 1;
						cpt_x = 0;
						jj++;
						while (string[jj] != ':')
							string4[cpt_x++] = string[jj++];
						string4[cpt_x] = '\0';
					}

				}
				else
				{
					// on recupère le nom du taxa

					for (jj = zz; jj < ii; jj++)
					{
						lesNoms[taxaPos - 1][xx++] = string[jj];
					}
					numero = taxaPos;
					lesNoms[taxaPos - 1][xx] = '\0';  // mettre la fin de chaine
					taxaPos++;  // augmenter le nombre de taxas pris
				}

			}
			else if (string[ii] == ',' || string[ii] == ')')
			{
				xx = 0;
				zz = ii + 1;   // faire pointer sur le prochain noeud
				for (jj = a3; jj < ii; jj++)
				{
					string1[xx++] = string[jj];
				}
				string1[xx++] = '\0';
				longueur = atof(string1);
				ARETE[aretePos++] = VertexNumber;
				ARETE[aretePos++] = numero;
				LONGUEUR[(aretePos - 1) / 2] = longueur;

				if (boot_value == 1) {
					//printf("\n%d--%d : %lf (%s)",VertexNumber,numero,longueur,string4);
					boot_value = 0;
				}
			}

		}

		// fin for pour traiter noeud
		//transcrire la nouvelle chaine
		xx = 0;
		for (jj = 0; jj < (int)a1; jj++)
		{
			string2[xx++] = string[jj];
		}

		// ecrire le vertex
		//	char buffer[50];
		itoa(VertexNumber, string1, 10);
		string2[xx++] = '%';   // indiquer que c'est un noeud
		for (jj = 0; jj < (int)strlen(string1); jj++)
		{
			string2[xx++] = string1[jj];
		}

		int temoin = 0;

		// transcrire la fin
		for (jj = a2 + 1; jj <= a; jj++)  // il faut voir si c'est  <= a ou c < a
		{
			if (string[jj] != ':' && temoin == 0) {
				string2[xx++] = '|';
			}
			temoin = 1;
			//if(temoin==1)
			string2[xx++] = string[jj];

		}

		// remplacer string par string2 en remplacant les pointeurs

		tempString = string;
		string = string2;
		string2 = tempString;
		tempString = 0;
		a = xx;  // mettre la longueur à jour 

	} // fin du while pour traiter toute la string

	/*for( jj=n;jj>0;jj--)
		strcpy(lesNoms[jj],lesNoms[jj-1]);*/

	int root_existance = -1;
	//for( jj=n;jj>0;jj--){
	for (jj = 0; jj < n; jj++) {
		//strcpy(lesNoms[jj],lesNoms[jj+1]);
		if (strcmp(lesNoms[jj], "Root") == 0)
			root_existance = jj;
	}
	/*printf("\nRoot_existance = %d",root_existance);
	for( jj=0;jj<n;jj++){
		printf("\n(%d) - %s",jj,lesNoms[jj]);
	}
	*/
	/*
		printf("\n\naretePos=%d\n",aretePos);
		for(i=1;i<=2*(2*n-3);i++){
			printf("\n%d",ARETE[i]);
		}

		*/
	ARETE[aretePos++] = 0;
	ARETE[aretePos++] = 0;

	/*for(i=1;i<=2*(na);i++){
		ARETE[i-1] = ARETE[i];
	}*/

	i = 0;
	int cpt_branches = 0;
	do {
		i++; cpt_branches++;
		//	printf("\n%d : %d",cpt_branches,ARETE[i]);
		ARETE[i - 1] = ARETE[i];
		i++;
		//	printf("--%d",ARETE[i]);
		ARETE[i - 1] = ARETE[i];
	} while (ARETE[i] != 0);

	for (i = 1; i <= na; i++) {
		LONGUEUR[i - 1] = LONGUEUR[i];
	}

	//printf("\nNombre de branches dans l'arbre : %d",na);
	/*
	printf("\n\naretePos=%d,n=%d\n",aretePos,n);
	for(i=1;i<=na;i++){
		printf("\n%d : %d-%d",i,ARETE[2*i-2],ARETE[2*i-1]);
	}*/
	//== recherche des branches connexes a celle de la racine;
	if (root_existance > 0) {
		int noeud_interne = -1;
		for (i = 1; i <= na; i++) {
			if (ARETE[2 * i - 1] == root_existance)
				noeud_interne = ARETE[2 * i - 2];
			if (ARETE[2 * i - 2] == root_existance)
				noeud_interne = ARETE[2 * i - 1];
		}
		//printf("\n[%d,%d]",noeud_interne,root_existance);
		for (i = 1; i <= na; i++) {
			if ((ARETE[2 * i - 1] != root_existance) && (noeud_interne == ARETE[2 * i - 2])) {
				LONGUEUR[i - 1] = 50;
				//		printf("\n[%d,%d]",noeud_interne,ARETE[2*i-1]);
			}
			if ((ARETE[2 * i - 2] != root_existance) && (noeud_interne == ARETE[2 * i - 1])) {
				LONGUEUR[i - 1] = 50;
				//		printf("\n[%d,%d]",noeud_interne,ARETE[2*i-2]);
			}
		}


	}


	//=== on teste si il y a un noeud de degre 2 et un noeud de degre 1

/*	printf("\nn=%d",n);
	for(i=1;i<=na;i++){
		printf("\n%d : %d-%d --> %lf",i,ARETE[2*i-1],ARETE[2*i-2],LONGUEUR[i-1]);
	}
	printf("\n");
  */
	int* tableau = (int*)malloc((2 * n) * sizeof(int));
	int deg2 = -1, deg1 = -1;
	for (i = 1; i <= 2 * n; i++)
		tableau[i - 1] = 0;
	for (i = 1; i <= na; i++) {
		tableau[ARETE[2 * i - 1]]++;
		tableau[ARETE[2 * i - 2]]++;
	}
	//	printf("\n");
	i = n + 1;
	while (tableau[i] > 0) {
		//printf("\n%d-->%d",i,tableau[i]);
		if (tableau[i] == 2) deg2 = i;
		if (tableau[i] == 1) deg1 = i;
		i++;
	}
	//printf("\ndeg2=%d , deg1=%d",deg2,deg1);
	int pos_racine = -1;

	for (i = 1; i <= na; i++) {
		if (ARETE[2 * i - 1] == deg1) { ARETE[2 * i - 1] = deg2; pos_racine = i;  break; }
		if (ARETE[2 * i - 2] == deg1) { ARETE[2 * i - 2] = deg2; pos_racine = i;  break; }
	}
	if (pos_racine != -1) {
		LONGUEUR[pos_racine - 1] = 100;
	}
	else if (deg2 != -1) {
		int pos1 = -1, pos2 = -1;
		for (i = 1; i <= na; i++) {
			if ((ARETE[2 * i - 1] == deg2 || ARETE[2 * i - 2] == deg2) && pos1 == -1) { pos1 = i; }
			else if ((ARETE[2 * i - 1] == deg2 || ARETE[2 * i - 2] == deg2) && pos1 != -1) { pos2 = i; }
		}
		//printf("\npos1=%d, pos2=%d",pos1,pos2);

		//== modification de la branche pos1
		if (ARETE[2 * pos1 - 1] == deg2)
			ARETE[2 * pos1 - 1] = (ARETE[2 * pos2 - 1] == deg2) ? ARETE[2 * pos2 - 2] : ARETE[2 * pos2 - 1];
		else
			ARETE[2 * pos1 - 2] = (ARETE[2 * pos2 - 1] == deg2) ? ARETE[2 * pos2 - 2] : ARETE[2 * pos2 - 1];

		//== suppression de la branche pos2
		for (i = pos2; i <= na; i++) {
			LONGUEUR[i - 1] = LONGUEUR[i];
			ARETE[2 * i - 1] = ARETE[2 * (i + 1) - 1];
			ARETE[2 * i - 2] = ARETE[2 * (i + 1) - 2];
		}
		na--;

	}

	/*printf("\nn=%d",n);
	for(i=1;i<=na;i++){
		printf("\n%d : %d-%d --> %lf",i,ARETE[2*i-1],ARETE[2*i-2],LONGUEUR[i-1]);
	}
	printf("\n");*/
	(*kt) = 2 * n - 3 - (*kt);

	//printf ("\nkt=%d,na=%d,n=%d",*kt,na,n);






	return n;
	//=== on teste si il y a un noeud de degre 2 et un noeud de degre 1

	//int * tableau = (int*)malloc((2*n)*sizeof(int));
	//int deg2=0,deg1=0;
	for (i = n + 1; i <= 2 * n; i++)
		tableau[i - 1] = 0;
	for (i = 1; i <= 2 * n - 3 - (*kt); i++) {
		tableau[ARETE[2 * i - 1]]++;
		tableau[ARETE[2 * i - 2]]++;
	}
	//	printf("\n");
	i = n + 1;
	while (tableau[i] > 0) {
		printf("\n%d-->%d", i, tableau[i]);
		if (tableau[i] == 2) deg2 = i;
		if (tableau[i] == 1) deg1 = i;
		i++;
	}

	for (i = 1; i <= 2 * n - 3 - (*kt); i++) {
		if (ARETE[2 * i - 1] == deg1) ARETE[2 * i - 1] = deg2;
		if (ARETE[2 * i - 2] == deg1) ARETE[2 * i - 2] = deg2;
	}

	printf("\nn=%d", n);
	for (i = 1; i <= na; i++) {
		printf("\n%d : %d-%d --> %lf", i, ARETE[2 * i - 1], ARETE[2 * i - 2], LONGUEUR[i - 1]);
	}



	free(string);
	free(string1);
	free(string2);
	free(newick);
	//free(string3);



	return n;

}

//===============================================================================================
//
//===============================================================================================
int lectureNewick_(char* fichier, long int* ARETE, double* LONGUEUR, char** lesNoms, int* kt)
{
	int n;
	printf("\nlecture Newick : ");
	// Ce sous programme permet de lire un arbre au format newick et de le transcrire dans
	// des vecteurs arete-longueur commencant à 1
	// ATTENTION: les noms commencent à 0
	// 7 octobre 2004
	// Elmaestro

	// TODO: Add your command handler code here
	int FAIL = -1;
	int i, j, j1, k, NewickValid = 1, a, a1, a2, a3, VertexNumber, numero;
	char symbol, * string, * string1, * string2/* *string3,c ,**Et*/;
	int taxaPos; // le nombre de taxas recupéré
	int aretePos; // le nombre d'aretes recupéré
	char symbolOld = ' ';
	int zz, xx, jj, ii;
	double longueur;
	char* tempString;
	//long int *ARETE;
	//double *LONGUEUR;


	FILE* data;
	if ((data = fopen(fichier, "r")) == 0) { printf("\n%s:Open Failed....", fichier); return FAIL; }

	//Correctness of the Newick format verification
	i = 0;
	n = 0;

	while ((symbol = getc(data)) != EOF)
	{
		if (symbol == ':') i++;
		if (symbol == ':' && symbolOld != ')') n++;
		symbolOld = symbol;
	}
	fseek(data, 0, 0);

	(*kt) = i;

	if (i == 0) { printf("Incorrect Newick file format. Edge lengths must be indicated after a ':' characters."); fclose(data); return FAIL; }

	j = 0;
	while ((symbol = getc(data)) != EOF)
	{
		if (symbol == '(') j++;
	}
	fseek(data, 0, 0);
	j1 = 0;
	while ((symbol = getc(data)) != EOF)
	{
		if (symbol == ')') j1++;
	}
	fseek(data, 0, 0);

	// verification des arités de l'arbre
	if (j1 != j) { printf("Incorrect Newick file format. Number of right parentheses must be equal to number of left parentheses."); fclose(data); return FAIL; }
	//else if (j!=n-2) { printf("Incorrect Newick file format. Only trees with vertices of degree 1 and 3 are allowed by T-REX."); fclose (data); return FAIL;}

	k = 0;
	while ((symbol = getc(data)) != EOF)
	{
		if (symbol == ',') k++;
	}
	fseek(data, 0, 0);
	//if (k!=(n-1)) { printf("Incorrect Newick file format. Number of objects must be equal to number of commas plus 1."); fclose (data); return FAIL;}

	a = 0;
	while ((symbol = getc(data)) != EOF)
	{
		if (symbol == ';') a++;
	}
	fseek(data, 0, 0);

	if (a == 0) { printf("Incorrect Newick file format. Newick string must be followed by a ';' character."); fclose(data); return FAIL; }
	else if (a > 1) { printf("Incorrect Newick file format. Newick string must contain (in the end) only one ';' character."); fclose(data); return FAIL; }

	a = 0;
	while ((symbol = getc(data)) == ' ') a++;
	fseek(data, 0, 0);
	if (symbol != '(') { printf("Incorrect Newick file format. Newick string must begin with a '(' character."); fclose(data); return FAIL; }

	a = 0;
	while ((symbol = getc(data)) != EOF)
	{
		if (symbol == '%') a++;
	}
	fseek(data, 0, 0);
	if (a > 0) { printf("Incorrect Newick file format. Newick string cannot contain '%%' character."); fclose(data); return FAIL; }

	while ((symbol = getc(data)) != EOF)
	{
		if ((symbol == '(') || (symbol == ','))
		{
			symbol = getc(data); a = 0;
			if ((symbol != '(') && (symbol != ',') && (symbol != ';') && (symbol != ':'))
			{
				fseek(data, -1, 1); while (((symbol = getc(data)) != ':') && (symbol != EOF)) a++;
			}
			else fseek(data, -1, 1);
			if (a > 50) { printf("Incorrect Newick file format. Names of objects must not exceed 50 characters."); fclose(data); return FAIL; }
		}
	}
	fseek(data, 0, 0);

	string = (char*)malloc((1000 * n) * sizeof(char));
	string2 = (char*)malloc((1000 * n) * sizeof(char));
	//string3 = (char*)malloc((1000*n) * sizeof(char));
	string1 = (char*)malloc((2000) * sizeof(char));

	if ((string == NULL) || (string1 == NULL) || (string2 == NULL)/*||(string3 == NULL)*/)
	{
		printf("Input data are too large or not a correct Newick file chosen"); fclose(data); return FAIL;
	}

	a = 0;
	while ((symbol = getc(data)) != EOF)
	{
		if ((symbol != ' ') && (symbol != '\n') && (symbol != '\t')) { string[a++] = symbol; }
	}
	fclose(data);
	printf(" fclose(data)");

	k = 0; VertexNumber = n;
	//a1 = 0;
	//a2 = 0;
	taxaPos = 1;    // nous allons commencer à mettre les taxas à la position 1
	aretePos = 1;
	while (string[0] == '(')   // traiter toute la chaine
	{
		a1 = 0;
		a2 = 0;
		while (string[a2] != ')')  // traiter la paire () la plus profonde
		{
			if (string[a2] == '(') a1 = a2;  // retrouver ;a parenthèse ouvrante
			a2++;
		}


		// a   => contient la longueur de la chaine
		// a1  => contient le debut d'un noeud à traiter
		// a2  => contient la fin d'un noeud à traiter
		// a3  => délimite le noeud et sa longueur

		zz = a1 + 1;
		VertexNumber++;  // augmenter le nombre de noeuds
		for (ii = a1 + 1; ii <= a2; ii++)
		{// decortiquer cette chaine

			if (string[ii] == ':')
			{
				xx = 0;
				a3 = ii + 1;

				if (string[zz] == '%')
				{ // cela veut dire que c'est un  noeud que l'on traite

					for (jj = zz + 1; jj < ii; jj++)
					{
						string1[xx++] = string[jj];
					}
					string1[xx++] = '\0';
					numero = atoi(string1);
				}
				else
				{
					// on recupère le nom du taxa

					for (jj = zz; jj < ii; jj++)
					{
						lesNoms[taxaPos - 1][xx++] = string[jj];
					}
					numero = taxaPos;
					lesNoms[taxaPos - 1][xx] = '\0';  // mettre la fin de chaine
					taxaPos++;  // augmenter le nombre de taxas pris
				}

			}
			else if (string[ii] == ',' || string[ii] == ')')
			{
				xx = 0;
				zz = ii + 1;   // faire pointer sur le prochain noeud
				for (jj = a3; jj < ii; jj++)
				{
					string1[xx++] = string[jj];
				}
				string1[xx++] = '\0';
				longueur = atof(string1);
				ARETE[aretePos++] = VertexNumber;
				ARETE[aretePos++] = numero;
				LONGUEUR[(aretePos - 1) / 2] = longueur;
			}

		}

		// fin for pour traiter noeud
		//transcrire la nouvelle chaine
		xx = 0;
		for (jj = 0; jj < (int)a1; jj++)
		{
			string2[xx++] = string[jj];
		}

		// ecrire le vertex
		//	char buffer[50];
		itoa(VertexNumber, string1, 10);
		string2[xx++] = '%';   // indiquer que c'est un noeud
		for (jj = 0; jj < (int)strlen(string1); jj++)
		{
			string2[xx++] = string1[jj];
		}

		// transcrire la fin
		for (jj = a2 + 1; jj <= a; jj++)  // il faut voir si c'est  <= a ou c < a
		{
			string2[xx++] = string[jj];
		}

		// remplacer string par string2 en remplacant les pointeurs

		tempString = string;
		string = string2;
		string2 = tempString;
		tempString = 0;
		a = xx;  // mettre la longueur à jour 



	} // fin du while pour traiter toute la string


	ARETE[aretePos++] = 0;
	ARETE[aretePos++] = 0;

	//fclose (data);

	for (i = 1; i <= 2 * n - 3; i++) {
		LONGUEUR[i - 1] = LONGUEUR[i];
	}
	for (i = 1; i <= 2 * (2 * n - 3); i++) {
		ARETE[i - 1] = ARETE[i];
	}

	(*kt) = 2 * n - 3 - (*kt);

	free(string);
	free(string1);
	free(string2);
	//free(string3);

	return n;
}

int readNextTree(int newickFormat, const char* infile, const char* outfile, int position) {

	FILE* in = fopen(infile, "r");
	FILE* out = fopen(outfile, "w+");
	int car, nbTree = 1, i, j;
	double val;
	int nbElt;
	char nom[30];

	if (newickFormat == 1) {
		//= rechercher la position
		while (nbTree < position) {
			car = fgetc(in);
			if (car == ';')
				nbTree++;
			if (car == EOF) {
				fclose(in);
				fclose(out);
				a = 0;
				return -1;
			}
		}
		//	printf("===>nbTree = %d",nbTree);

			//= rechercher la paranthese
		do {
			car = fgetc(in);
			if (car == EOF) break;
		} while (car != '(');

		//= copie de l'arbre
		if (car == '(') {
			do {
				if (car == -1) {
					fclose(in);
					fclose(out);
					a = 0;
					exit(-1);
				}
				if (((char)car != '\n') && ((char)car != '\r'))
					//fputc((char)car,out);
					fprintf(out, "%c", car);
				//			printf("%c",car);
				car = fgetc(in);
			} while (car != ';');
			fputc((char)car, out);
			fclose(in);
			fclose(out);
			//printf("\nApres fclose");
			a = 0;
			return 0;
		}
		else {
			fclose(in);
			fclose(out);
			a = 0;
			return -1;
		}
	}
	else {
		//= rechercher la position
		while (nbTree < position) {
			if (fscanf(in, "%d", &nbElt) == EOF) {
				fclose(in);
				fclose(out);
				a = 0;
				return -1;
			}
			else {
				//			printf("%d",nbElt);
				for (i = 0; i < nbElt; i++) {
					fscanf(in, "%s", nom);
					//printf("\n%s",nom);
					for (j = 0; j < nbElt; j++) {
						fscanf(in, "%lf", &val);
						//printf("%lf ",val);
					}
				}
			}
			nbTree++;
			//		printf("allo");
		}

		//= copie de l'arbre
		fscanf(in, "%d", &nbElt);
		fprintf(out, "%d\n", nbElt);
		for (i = 0; i < nbElt; i++) {
			fscanf(in, "%s", nom);
			fprintf(out, "%s\t", nom);
			for (j = 0; j < nbElt; j++) {
				fscanf(in, "%lf", &val);
				fprintf(out, "%lf ", val);
			}
		}
		fclose(in);
		fclose(out);
		a = 0;
	}
}

void TrierMatrices(double** DISS, char** NomsDISS, char** NomsADD, int n)
{
	int ligne, colonne, i, j;
	double** DISS_;
	char** NomsDISS_;
	char noms[50];
	int* table;
	int trouve;

	table = (int*)malloc((n + 1) * sizeof(int));
	DISS_ = (double**)malloc((n + 1) * sizeof(double*));
	NomsDISS_ = (char**)malloc((n + 1) * sizeof(char*));

	for (i = 0; i <= n; i++)
	{
		DISS_[i] = (double*)malloc((n + 1) * sizeof(double));
		NomsDISS_[i] = (char*)malloc((n + 1) * sizeof(15));
	}

	for (ligne = 1; ligne <= n; ligne++)
	{
		strcpy(noms, NomsADD[ligne - 1]);
		trouve = 0;
		for (colonne = 1; colonne <= n; colonne++)
		{
			if (strcmp(noms, NomsDISS[colonne - 1]) == 0)
			{
				trouve = 1;
				table[ligne] = colonne;
				strcpy(NomsDISS_[ligne - 1], noms);
			}
		}
		if (trouve == 0)
		{
			printf("\n%s %s", noms, "is not in the gene matrix.This program must stop");
			exit(5);
		}
	}

	for (i = 1; i <= n; i++)
		for (j = 1; j <= i; j++) DISS_[i][j] = DISS_[j][i] = DISS[table[i]][table[j]];

	for (i = 1; i <= n; i++)
	{
		strcpy(NomsDISS[i - 1], NomsDISS_[i - 1]);
		for (j = 1; j <= n; j++) DISS[i][j] = DISS_[i][j];
	}

	for (i = 0; i <= n; i++)
	{
		free(DISS_[i]);
		free(NomsDISS_[i]);
	}

	free(DISS_); free(NomsDISS_);
}

void CopieMatriceDouble(double** source, double** destination, int taille)
{
	int i, j;
	for (i = 0; i <= taille; i++)
		for (j = 0; j <= taille; j++)
			destination[i][j] = source[i][j];
}

void CopieArrayString(char** source, char** dest, int taille) {

	int i;

	for (i = 0; i < taille; i++)
		strcpy(dest[i], source[i]);
}

void readMatrix(char* tmp, double** DISS, char** NOMS, int* n) {

	double val;
	int nbElt, i, j;
	char nom[50];
	FILE* in = fopen(tmp, "r");

	fscanf(in, "%d", &nbElt);
	//printf("%d\n",nbElt);
	(*n) = nbElt;
	for (i = 1; i <= nbElt; i++) {
		fscanf(in, "%s", nom);
		strcpy(NOMS[i - 1], nom);
		for (j = 1; j <= nbElt; j++) {
			fscanf(in, "%lf", &val);
			DISS[i][j] = val;
		}
	}
	fclose(in);
}

void EcrireMatrice(FILE* outmat, double** DISS, char** NOMS, int n) {

	int i, j, k;
	fprintf(outmat, "%d", n);
	a++;
	if (a > SEUIL) {
		fclose(outmat);
		exit(-1);
	}
	for (i = 1; i <= n; i++) {
		fprintf(outmat, "\n%s", NOMS[i - 1]);
		a++;
		if (a > SEUIL) {
			fclose(outmat);
			exit(-1);
		}
		int mk = (10 - strlen(NOMS[i - 1]));
		if (mk > 0)
			for (k = 0; k <= mk; k++) {
				fprintf(outmat, " ");
				a++;
				if (a > SEUIL) {
					fclose(outmat);
					exit(-1);
				}
			}
		for (j = 1; j <= n; j++) {
			fprintf(outmat, "%lf ", DISS[i][j]);
			a++;
			if (a > SEUIL) {
				fclose(outmat);
				exit(-1);
			}
		}
	}
}


/* Main function */

//= format : rf.exe inputfile outputfile tmp matrixfile
int main(int argc, char** argv)
{
	a = 0;
	//= 1 : input.txt
	//= 2 : output.txt
	//= 3 : tmp.txt
	//= 4 : matrice.txt
	int i, j, k, ii, MatrixNumber, RF, ** B, ** BI, * PLACE1, * PLACE2, m, m1, n, n_, tailleMax = 0, finalTaille;
	double** DI, ** D, ** DISS, ** DISS_opt, ** DISSopt, ** DISScopy, ** DISS_;
	double* LONGUEURS, * LONGUEURS_;
	int long* ARETES, * ARETES_;
	char** NOMS, ** NOMS_, ** NOMSopt, ** NOMS_opt, ** NOMScopy;
	int newickFormat = 0, ** resultat;
	FILE* data;
	FILE* outmat;
	FILE* output;
	int position = 1;
	int bidon, kt1, kt2;

	if (argc != 5) {
		printf("Erreur, nombre d'arguments invalide, usage : \n");
		printf("./rf input.txt ouput.txt tmp.txt matrices.txt");
	}
	//if ((data=fopen(argv[1],"r"))==0) { printf("\nFile %s was not found ", argv[1]); exit(5); } 
	//printf("\n\n");

	outmat = fopen(argv[4], "w+");
	output = fopen(argv[2], "w+");


	fprintf(output, "* --------------------------------------------------*\n");
	fprintf(output, "* Computation of the Robinson and Foulds            *\n");
	fprintf(output, "* topological distance between two (or more) trees. *\n");
	fprintf(output, "* ------------------------------------------------- *\n");


	//====================================================
	//= lecture du nombres d'elements du premier arbre
	//====================================================
	if (fileInNewickFormat(argv[1]) == 1) {
		newickFormat = 1;
		n = nbSpeciesNewick(argv[1]);
	}
	else
		n = nbSpeciesPhylip(argv[1]);
	//printf("\nFormat fichier = %s, \nNombre d'elements=%d",(newickFormat==1)?"Newick":"Philip",n);

	//==========================================================================
	//= recherche de la taille du plus grand arbre pour l'allocation de mémoire
	//==========================================================================
	MatrixNumber = 0;
	i = 0;
	while (1) {
		if (readNextTree(newickFormat, argv[1], argv[3], i + 1) == -1)
			break;
		if (newickFormat == 1)
			n_ = nbSpeciesNewick(argv[3]);
		else
			n_ = nbSpeciesPhylip(argv[3]);
		if (tailleMax < n_)
			tailleMax = n_;
		MatrixNumber++;
		i++;
	}
	if (newickFormat == 0)
		MatrixNumber--;

	//printf("\nnombre de matrice = %d(%d)",MatrixNumber,tailleMax);	
	if (MatrixNumber < 2)
	{
		printf("You must analyse at least two additive distance matrices\n\n"); exit(5);
	}

	tailleMax++;

	//===========================================================
	//=allocation dynamique de la mémoire pour les autres arbres
	//===========================================================
	DISS_ = createDoubleMatrix(2 * tailleMax);
	LONGUEURS_ = (double*)malloc((2 * tailleMax) * sizeof(double));
	ARETES_ = (long int*)malloc(2 * (2 * tailleMax) * sizeof(long int));
	NOMS_ = createStringMatrix(tailleMax, 50);
	resultat = (int**)malloc((MatrixNumber + 1) * sizeof(int*));
	for (i = 0; i < MatrixNumber + 1; i++)
		resultat[i] = (int*)malloc((MatrixNumber + 1) * sizeof(int));

	//=========================================================
	//=allocation dynamique de la mémoire pour le premier arbre
	//=========================================================
	DISS = createDoubleMatrix(2 * tailleMax);
	DISS_opt = createDoubleMatrix(2 * tailleMax);
	DISSopt = createDoubleMatrix(2 * tailleMax);
	DISScopy = createDoubleMatrix(2 * tailleMax);
	LONGUEURS = (double*)malloc((2 * tailleMax) * sizeof(double));
	ARETES = (long int*)malloc(2 * (2 * tailleMax) * sizeof(long int));
	NOMS = createStringMatrix(tailleMax, 50);
	NOMScopy = createStringMatrix(tailleMax, 50);
	NOMS_opt = createStringMatrix(tailleMax, 50);
	NOMSopt = createStringMatrix(tailleMax, 50);

	//=========================================================
	//=allocation dynamique de la mémoire pour le calcul de RF
	//=========================================================
	B = (int**)malloc((2 * tailleMax - 2) * sizeof(int*));
	BI = (int**)malloc((2 * tailleMax - 2) * sizeof(int*));
	PLACE1 = (int*)malloc((2 * tailleMax - 2) * sizeof(int));
	PLACE2 = (int*)malloc((2 * tailleMax - 2) * sizeof(int));

	for (i = 0; i < 2 * tailleMax - 2; i++)
	{
		B[i] = (int*)malloc((2 * tailleMax - 2) * sizeof(int));
		BI[i] = (int*)malloc((2 * tailleMax - 2) * sizeof(int));
	}


	for (k = 1; k <= MatrixNumber; k++) {

		if (readNextTree(newickFormat, argv[1], argv[3], k) == -1) { printf("no more tree in file"); exit(5); }
		if (newickFormat == 1) {
			n = nbSpeciesNewick(argv[3]);
			lectureNewick(argv[3], ARETES, LONGUEURS, NOMS, &kt1);
			TreeMatrix(DISS, ARETES, LONGUEURS, n, kt1);
			EcrireMatrice(outmat, DISS, NOMS, n);
			fprintf(outmat, "\n\n");
		}
		else {
			n = nbSpeciesPhylip(argv[3]);
			readMatrix(argv[3], DISS, NOMS, &n);
		}

		//========================================================
		//= 1) Lecture des autres arbres en matrice de distance 
		//= 2) Filtrage des matrices par rapport à la premiere
		//= 3) Triage des matrices par rapport à la premiere
		//= 4) Calcul de la distance de Robinson and Foulds 
		//========================================================

		//position++;
		for (ii = k + 1; ii <= MatrixNumber; ii++) {
			//	printf("\nLecture arbre #%d :",ii+1);
			if (readNextTree(newickFormat, argv[1], argv[3], ii) == -1) {
				printf("probleme lecture arbre");
				exit(5);
			}
			//	printf("ok");	
			if (newickFormat == 1)
				n_ = nbSpeciesNewick(argv[3]);
			else
				n_ = nbSpeciesPhylip(argv[3]);


			if (newickFormat == 1) {
				lectureNewick(argv[3], ARETES_, LONGUEURS_, NOMS_, &kt2);
				TreeMatrix(DISS_, ARETES_, LONGUEURS_, n_, kt2);
				//fprintf(outmat,"\n\n");
				//EcrireMatrice(outmat,DISS_,NOMS_,n_);
			}
			else
				readMatrix(argv[3], DISS_, NOMS_, &n_);
			/*
			printf("\n\nMatrice lue");
				for(i=1;i<=n;i++){
					printf("\n%s\t",NOMS_[i-1]);
					for(j=1;j<=n;j++)
						printf("%lf ",DISS_[i][j]);
				}
			*/

			CopieMatriceDouble(DISS, DISScopy, n);
			CopieArrayString(NOMS, NOMScopy, n);
			/*
					printf("\n\nMatrice initiale apres copie");
					for(i=1;i<=n;i++){
						printf("\n%s\t",NOMScopy[i-1]);
						for(j=1;j<=n;j++)
							printf("%lf ",DISScopy[i][j]);
					}
			*/
			filtrerMatrice(DISScopy, DISS_, NOMScopy, NOMS_, n, n_);

			/*		printf("\n\nMatrice initiale apres filtre");
					for(i=1;i<=n;i++){
						printf("\n%s\t",NOMScopy[i-1]);
						for(j=1;j<=n;j++)
							printf("%lf ",DISScopy[i][j]);
					}
			*/
			finalTaille = ecrireMatrice2(DISScopy, DISSopt, n, NOMScopy, NOMSopt);
			/*
					printf("\n\nMatrice initiale finale (%d)",finalTaille);
					for(i=1;i<=finalTaille;i++){
						printf("\n%s\t",NOMSopt[i-1]);
						for(j=1;j<=finalTaille;j++)
							printf("%lf ",DISSopt[i][j]);
					}
			*/
			finalTaille = ecrireMatrice2(DISS_, DISS_opt, n_, NOMS_, NOMS_opt);
			/*		printf("\n\nMatrice lue apres filtre (%d)",finalTaille);
					for(i=1;i<=finalTaille;i++){
						printf("\n%s\t",NOMS_opt[i-1]);
						for(j=1;j<=finalTaille;j++)
							printf("%lf ",DISS_opt[i][j]);
					}
				*/
			TrierMatrices(DISS_opt, NOMS_opt, NOMSopt, finalTaille);
			/*
			printf("\n\nMatrice lu apres tri");
				for(i=1;i<=finalTaille;i++){
					printf("\n%s\t",NOMS_opt[i-1]);
					for(j=1;j<=finalTaille;j++)
						printf("%lf ",DISS_opt[i][j]);
				}
		*/
		//========================================
		//= calcul de la distance RF
		//========================================		
			for (i = 0; i <= 2 * tailleMax - 3; i++) {
				PLACE1[i] = 0;
				PLACE2[i] = 0;
				for (j = 0; j <= tailleMax; j++)
					BI[i][j] = 0;
				B[i][j] = 0;
			}

			m = Bipartition_Table(DISSopt, B, PLACE1, finalTaille);
			m1 = Bipartition_Table(DISS_opt, BI, PLACE2, finalTaille);
			if (m1 == 0)
				exit(1);
			RF = Table_Comparaison(B, BI, PLACE1, PLACE2, m, m1, finalTaille);

			fprintf(output, "\nRF Distance between Tree %d and Tree %d = %d", k, ii, RF);
			//printf("\nRF Distance between tree %d and tree %d = %d",k,ii,RF);
			resultat[k][ii] = resultat[ii][k] = RF;
			//position++;
		}
	}

	fprintf(output, "\n\n");
	a++;
	if (a > SEUIL) {
		fclose(output);
		exit(-1);
	}
	for (i = 1; i <= MatrixNumber; i++) {
		fprintf(output, "\n%d\t", i);
		a++;
		if (a > SEUIL) {
			fclose(output);
			exit(-1);
		}
		//	printf("\n%d\t",i);
		for (j = 1; j <= MatrixNumber; j++) {
			resultat[i][i] = 0;
			fprintf(output, "%d ", resultat[i][j]);
			a++;
			if (a > SEUIL) {
				fclose(output);
				exit(-1);
			}
			//	printf("%d ",resultat[i][j]);
		}
	}
	//printf("\n");
	fclose(outmat);
	fclose(output);
	a = 0;
	return 0;
}




/* Computing a circular order X of n objects of the distance matrix D starting from the objects i1 and i2 */

void odp(double** D, int* X, int* i1, int* j1, int n)
{
	double S1, S;
	int i, j, k, a, * Y1;

	Y1 = (int*)malloc((n + 1) * sizeof(int));

	for (i = 1; i <= n; i++)
		Y1[i] = 1;

	X[1] = *i1;
	X[n] = *j1;
	if (n == 2) return;
	Y1[*i1] = 0;
	Y1[*j1] = 0;
	for (i = 0; i <= n - 3; i++)
	{
		a = 2;
		S = 0;
		for (j = 1; j <= n; j++)
		{
			if (Y1[j] > 0)
			{
				S1 = D[X[n - i]][X[1]] - D[j][X[1]] + D[X[n - i]][j];
				if ((a == 2) || (S1 <= S))
				{
					S = S1;
					a = 1;
					X[n - i - 1] = j;
					k = j;
				}
			}
		}
		Y1[k] = 0;
	}
	free(Y1);
}


/* Computation of an ordered bipartition table B(2n-3,n) and its rank list PLACE(2n-3)
  for the tree associated with an additive distance matrix D */

int Bipartition_Table(double** D, int** B, int* PLACE, int n)
{

	int i, j, k, l, l1, * MaxCol, * X, EdgeNumberPath, m, uv, PlaceNumber, edge, * Path, M, F;
	double S, DIS, DIS1, * LengthPath;
	double EPS = 1.e-5;
	double EPS1 = 1.e-2;

	/* Memory allocation */

	MaxCol = (int*)malloc((2 * n - 2) * sizeof(int));
	X = (int*)malloc((n + 1) * sizeof(int));
	LengthPath = (double*)malloc((2 * n - 2) * sizeof(double));
	Path = (int*)malloc((2 * n - 2) * sizeof(int));

	/* Computation of a circular order X for D */

	i = 1; j = n; odp(D, X, &i, &j, n);

	/* Initialization */
	for (i = 1; i <= 2 * n - 3; i++)
	{
		MaxCol[i] = 0;
		PLACE[i] = 0;
		for (j = 1; j <= n; j++)
			B[i][j] = 0;
	}
	B[1][X[2]] = 1; MaxCol[1] = X[2]; Path[1] = 1; PlaceNumber = 1;
	PLACE[1] = 1; LengthPath[1] = D[X[1]][X[2]]; EdgeNumberPath = 1; m = 1;


	/* The main loop */
	for (k = 2; k <= n - 1; k++)
	{
		/* Point 2.1 of the algorithm (see the referenced article by Makarenkov and Leclerc) */

		DIS = (D[X[1]][X[k]] + D[X[k]][X[k + 1]] - D[X[1]][X[k + 1]]) / 2;
		DIS1 = (D[X[1]][X[k + 1]] + D[X[k]][X[k + 1]] - D[X[1]][X[k]]) / 2;

		if ((DIS <= -EPS1) || (DIS1 <= -EPS1)) {
			printf("\n This is not an additive distance \n");
			free(MaxCol); free(X); free(LengthPath); free(Path); return 0;
		}
		if (DIS <= EPS) DIS = 0.0; if (DIS1 <= EPS) DIS1 = 0.0;

		S = 0.0; i = EdgeNumberPath; if (LengthPath[i] == 0.0) i--;
		while (S <= DIS - EPS)
		{
			if (i == 0) { S = DIS; break; }  /* checking the limit */
			S = S + LengthPath[i];
			i--;
		}

		/* Point 2.2 of the algorithm */

		if (fabs(S - DIS) <= EPS)
		{
			M = m + 2; DIS = S;
			if (i == 0) F = 1;
			else if (i == EdgeNumberPath) F = 2;
			else { M--; F = 3; }
		}
		else { M = m + 2; F = 0; }


		if (M == m + 2)
		{
			if (F == 0) {
				uv = Path[i + 1]; EdgeNumberPath = i + 2; LengthPath[i + 1] = S - DIS; LengthPath[i + 2] = DIS1;
				Path[i + 1] = m + 2; Path[i + 2] = m + 1;
			}
			else if (F == 1) {
				uv = Path[1]; EdgeNumberPath = 2; LengthPath[1] = 0.0; LengthPath[2] = DIS1;
				Path[1] = m + 2; Path[2] = m + 1;
			}
			else if (F == 2) {
				uv = Path[EdgeNumberPath]; EdgeNumberPath = EdgeNumberPath + 1; LengthPath[EdgeNumberPath] = DIS1;
				Path[EdgeNumberPath - 1] = m + 2; Path[EdgeNumberPath] = m + 1;
			}

			for (j = 1; j <= n; j++)
				B[m + 2][j] = B[uv][j];
			MaxCol[m + 2] = MaxCol[uv];
		}

		else
		{
			EdgeNumberPath = i + 1; LengthPath[i + 1] = DIS1; Path[i + 1] = m + 1;
		}

		/* Point 2.3 of the algorithm */

		for (j = 1; j <= EdgeNumberPath; j++)
			B[Path[j]][X[k + 1]] = 1;

		/* Point 2.4 of the algorithm */

		for (j = 1; j <= EdgeNumberPath; j++)
			if (MaxCol[Path[j]] < X[k + 1]) MaxCol[Path[j]] = X[k + 1];

		/* Point 2.5 of the algorithm */

		for (j = PlaceNumber; j >= 1; j--)
			PLACE[j + 1] = PLACE[j];
		PLACE[1] = m + 1; PlaceNumber++;

		if (M == m + 2) {
			i = 2;
			while (PLACE[i] != uv)
				i++;
			for (j = PlaceNumber; j >= i + 1; j--)
				PLACE[j + 1] = PLACE[j];
			PLACE[i + 1] = m + 2; PlaceNumber++;
		}

		i = M - 1; edge = 2;
		do
		{
			if (PLACE[i] == Path[edge])
			{
				edge++; j = i + 1;
				while (X[k + 1] > MaxCol[PLACE[j]])
					j++;
				if (j > i + 1)
				{
					l1 = PLACE[i];
					for (l = i + 1; l <= j - 1; l++)
						PLACE[l - 1] = PLACE[l];
					PLACE[j - 1] = l1;

				}
			}
			i--;
		} while (i != 0);

		m = M;
	}
	/* memory liberation */

	free(MaxCol);
	free(X);
	free(LengthPath);
	free(Path);

	return m;

}

/* Integer function returning the value of the Robinson and Foulds topological distance
   between two trees given their ordered bipartition table B(2n-3,n) and B1(2n-3,n)
   and its rank lists PLACE(2n-3) and  PLACE1(2n-3) */


int Table_Comparaison(int** B, int** B1, int* PLACE, int* PLACE1, int m, int m1, int n)
{
	int RF = 0, i, p, p1;

	p = 1; p1 = 1;

	while ((p <= m) && (p1 <= m1))
	{
		i = n;
		while ((B[PLACE[p]][i] == B1[PLACE1[p1]][i]) && (i > 1))
			i--;
		if (i == 1) { RF = RF + 1; p++; p1++; }
		else if (B[PLACE[p]][i] > B1[PLACE1[p1]][i]) p1++;
		else p++;

	}
	RF = (m - RF) + (m1 - RF);

	return RF;
}










