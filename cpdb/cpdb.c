//
//  cpdb.c
//  cpdb
//
//  Created by Gokhan SELAMET on 28/09/2016.
//  Copyright © 2016 Gokhan SELAMET. All rights reserved.
//

#include <string.h>
#include "cpdb.h"
#include <stdlib.h>

#define INC_CAPACITY 15;

int extractStr(char *dest, const char *src, int begin, int end) ;
void getAtomType (const char *line, char *out);
void getResType (const char *line, char *out);
void getAtomElement (const char *line, char *out);
void getCoordinates (const char *line, v3 *out);
void getAtomId (const char *line, int *atomId);
void getResidueId (const char *line, int *residueId);
void getAlternativeLoc (const char *line, char *altLoc);
void getChainId (const char *line, char *chainId);
void getOccupancy (const char *line, float *occupancy);
void getTempFactor (const char *line, float *tempFactor);
void updateResiduePointers (chain *C);
void updateAtomPointers (residue *R);

pdb* initPDB() {
	return calloc(1, sizeof(pdb));
}

int parsePDB (char *pdbFilePath, pdb *P , char *options) {
    chain *C;
    char altLocFlag = 1; // Default: skip alternate location atoms on
    char hydrogenFlag = 1; // Default: skip hydrogene atoms on
	int errorcode = 0;
	C = calloc(1, sizeof (chain) );
    char line[128];
    FILE* pdbFile;
    char resType[5], atomType[5], atomElement[3], chainId, altLoc ;
    int resId, atomId, length;
    float tempFactor, occupancy;
	int resIdx = 0;
	int atomIdx = 0;
    v3 coor;
    residue *currentResidue = NULL;
    chain *currentChain = NULL;
    atom *currentAtom = NULL;
	
	// Setting Parsing Options
	int chindex = 0;
	char ch = options[chindex++];

	while (ch) {
		switch( ch ) {
			case 'h' :
				hydrogenFlag = 0;
				break;
			case 'l' :
				altLocFlag = 0;
				break;
			case 'a':
				hydrogenFlag = 0;
				altLocFlag = 0;
				break;
			default:
				fprintf(stderr, "invalid option: %c \n", ch );
				errorcode ++;
		}
		ch=options[chindex++];
	}

    pdbFile = fopen(pdbFilePath, "r");
    length = 0;
	
	if (pdbFile == NULL) {
		perror("pdb file can not read");
		exit (2);
	}

    while (fgets(line, sizeof(line), pdbFile)) {
        if (!strncmp(line, "ATOM  ", 6)) {
            getCoordinates(line, &coor);
            getAtomType(line, atomType);
            getResType(line, resType);
            getAtomElement(line, atomElement);
            getAtomId(line, &atomId);
            getResidueId(line, &resId);
            getChainId(line, &chainId);
            getAlternativeLoc(line, &altLoc);
            getTempFactor(line, &tempFactor);
            getOccupancy(line, &occupancy);
            if ( hydrogenFlag && !strncmp(atomElement, "H", 3) ) continue;
            if ( altLocFlag && !(altLoc == 'A' || altLoc == ' ' ) ) continue;
			
            // Chain Related
            if (currentChain == NULL || currentChain->id != chainId) {
                appendChaintoPdb(P, (chain){.id=chainId, .residues=NULL, .size=0,
											.__capacity=0});
                currentChain = &(P->chains[P->size - 1]);
                currentAtom = NULL;
                currentResidue = NULL;
            }
            
            // Residue Related
            if (currentResidue == NULL || currentResidue->id != resId) {
                appendResiduetoChain(currentChain, (residue)
									{.id=resId, .idx=resIdx++, .type={0},
									.atoms=NULL, .size=0, .__capacity=0,
									.next=NULL, .prev=currentResidue} );
                currentResidue = &(currentChain->residues[currentChain->size-1]);
                strncpy(currentResidue->type, resType, 5);
                if (currentResidue->prev) currentResidue->prev->next=currentResidue;
            }
            
            //Atom Related
            if (currentAtom == NULL || currentAtom->id != atomId) {
                appendAtomtoResidue(currentResidue, (atom)
										{.id=atomId, .idx=atomIdx++, .coor=coor,
										.type={0}, .element={0}, .tfactor=tempFactor,
										.occupancy=occupancy, .next=NULL, .prev=currentAtom,
										.res=currentResidue});
                currentAtom = &(currentResidue->atoms[currentResidue->size-1]);
                strncpy(currentAtom->type, atomType, 5);
                strncpy(currentAtom->element, atomElement, 3);
                if (currentAtom->prev) currentAtom->prev->next = currentAtom;
            }
            
        }
        else if (!strncmp(line, "ENDMDL", 6) || !strncmp(line,  "TER   ", 6))
            break;
    }
    fclose(pdbFile);
    return errorcode;
}

int writePDB (const char *filename, const pdb *P) {
	FILE *out;
	out = fopen(filename, "w");
	
	if (out == NULL) {
		perror("pdb file write error");
		exit(1);
	}
	writeFilePDB (out, P);
	fclose(out);
	return 0;
}

int writeFilePDB (FILE *F, const pdb *P) {
	chain *C = NULL;
	atom *A = NULL;
	int chainId;
	for (chainId = 0 ; chainId < P->size; chainId++) {
		C = &P->chains[chainId];
		A = &C->residues[0].atoms[0];
		while (A != NULL) {
			fprintf (F, "ATOM  %5d %-4s%c%-3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
					 A->id, A->type,' ', A->res->type, C->id, A->res->id,
					 A->coor.x, A->coor.y, A->coor.z, A->occupancy, A->tfactor, A->element);
			A = A->next;
		}
	}
	return 0;
}

int printPDB(pdb *P) {
	writeFilePDB(stdout, P);
	return 0;
}

void freePDB(pdb *P) {
	int i,j;
	chain *C;
	residue *R;
	for (i=0; i<P->size; i++) {
		C = & P->chains[i];
		for (j=0; j<C->size; j++) {
			R = & C->residues[j];
			free (R->atoms);
		}
		free(C->residues);
	}
	free(P->chains);
	free(P);
}

int extractStr(char *dest, const char *src, int begin, int end) {
    int length;
    while (src[begin] == ' ') begin++;
    while (src[end] == ' ') end--;
    length = end - begin + 1;
    strncpy(dest, src + begin, length);
    dest[length] = 0;
    return length;
}

void getCoordinates (const char *line, v3 *out) {
    char coorx[9] = {0};
    char coory[9] = {0};
    char coorz[9] = {0};
    strncpy (coorx, &line[30], 8);
    strncpy (coory, &line[38], 8);
    strncpy (coorz, &line[46], 8);
    out->x = atof (coorx);
    out->y = atof (coory);
    out->z = atof (coorz);
}

void getAtomType (const char *line, char *out) {
    extractStr(out, line, 12, 16);
}

void getResType (const char *line, char *out) {
    extractStr(out, line, 17, 19);
}

void getAtomElement (const char *line, char *out) {
    extractStr(out, line, 76, 78);
}

void getAtomId (const char *line, int *atomId) {
    char _temp[6]={0};
    extractStr(_temp, line, 6, 10);
    *atomId = atoi(_temp);
}

void getResidueId (const char *line, int *residueId) {
    char _temp[5]={0};
    extractStr(_temp, line, 22, 25);
    *residueId = atoi(_temp);
}

void getAlternativeLoc(const char *line, char *altLoc) {
    *altLoc = line[16];
}

void getChainId(const char *line, char *chainId) {
    *chainId = line[21];
}

void getOccupancy (const char *line, float *occupancy) {
    char _temp[7]={0};
    extractStr(_temp, line, 54, 59);
    *occupancy = atof(_temp);
}

void getTempFactor (const char *line, float *tempFactor){
    char _temp[7]={0};
    extractStr(_temp, line, 60, 65);
    *tempFactor = atof(_temp);
}

void appendChaintoPdb (pdb *P, chain newChain) {
	if ( P->__capacity == 0) {
		P->chains = calloc (1, sizeof(chain));
		P->__capacity++;
	}
	if (P->size == P->__capacity) {
        P->__capacity += INC_CAPACITY;
		chain *newChain;
		newChain = realloc (P->chains, P->__capacity * sizeof(chain));
		if (newChain == 0) {
			fprintf (stderr, "Allocation error \n");
			exit(1);
		}
		P->chains = newChain;
    }
    P->chains [P->size++] = newChain;
}

void appendResiduetoChain (chain *C, residue newResidue) {
	if ( C->__capacity == 0)
		C->residues = calloc (1, sizeof(residue));
	if ( C->size == C->__capacity ) {
        C->__capacity += INC_CAPACITY;
		residue *newresidues;
		newresidues = realloc(C->residues, C->__capacity * sizeof(residue));
		if (newresidues == 0) {
			fprintf (stderr, "Allocation error \n");
			exit(1);
		}
		C->residues = newresidues;
		memset(&C->residues[C->size], 0, (C->__capacity - C->size)*sizeof(residue));
		C->residues [C->size++] = newResidue;
        updateResiduePointers(C);
	} else {
    C->residues [C->size++] = newResidue;
	}
}

void appendAtomtoResidue (residue *R, atom newAtom) {
	if ( R->__capacity == 0)
		R->atoms = calloc (1, sizeof(atom));
	if ( R->size == R->__capacity) {
        R->__capacity += INC_CAPACITY;
		atom *newatoms;
        newatoms = realloc(R->atoms, R->__capacity * sizeof(atom));
		if (newatoms == 0) {
			fprintf (stderr, "Allocation error \n");
			exit(1);
		}
		R->atoms = newatoms;
		memset(&R->atoms[R->size], 0, (R->__capacity - R->size)*sizeof(atom));
        R->atoms[R->size++] = newAtom;
        updateAtomPointers(R);
    } else {
    R->atoms[R->size++] = newAtom;
    }
}

void updateResiduePointers (chain *C) {
	int i, j;
	residue *p = NULL;
	residue *n = NULL;
	atom *a = NULL;
	for (i=0; i<C->size; i++) {
		n = & C->residues[i];
		n->prev = p;
		if (p) p->next = n;
		p = n;
		for (j=0; j<n->size; j++) {
			a = & n->atoms[j];
			a->res = n;
		}
	}
	if (n) n->next = NULL;
    
}

void updateAtomPointers (residue *R) {
	int i;
	atom *p = NULL;
	atom *c = NULL;
	if (R->prev) p = & R->prev->atoms[R->prev->size-1];

	for (i=0; i<R->size; i++){
		c = & R->atoms[i];
		if (p) p->next = c;
		c->prev = p;
		p = c;
	}

	if (R->next) {
		c= & R->next->atoms[0];
		c->prev = p;
		p->next = c;
	}
	else  {
		if(c) c->next = NULL;
	}
}

// TODO
/*
void removeResidue();
void removeAtom();
void insertAtom();
void insertResidue();
 */
