/* Entete des utilitaires de PBM+ */
/* E.B.  12/98                   */

typedef unsigned char bit;
typedef double gray;

char pm_getc(FILE* file);
bit pm_getbit(FILE* file);
unsigned char pm_getrawbyte(FILE* file);
int pm_getint( FILE* file);

void pm_erreur(char *);
