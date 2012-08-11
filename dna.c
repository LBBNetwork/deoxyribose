/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

 
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

/* FFUUU WINDOWS */
#ifdef _WIN32
#ifdef ERROR
#undef ERROR
#endif
#endif

#define __DNA_VERSION__ "0.1pre"

/* debug definitions */
#define UNIMPLEMENTED    debug_printf(DBGFLTR_MISC, "WARNING:  %s at %s:%d is UNIMPLEMENTED!\n",__FUNCTION__,__FILE__,__LINE__);
#define FATAL(fmt, ...)  debug_printf(DBGFLTR_RELEASE, "(%s:%d) FATAL ERROR (Aborting): " fmt, __FILE__, __LINE__, ##__VA_ARGS__), exit(-1);
#define ERR(fmt, ...)    debug_printf(DBGFLTR_ERR, "(%s:%d) ERROR: " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define WARN(fmt, ...)   debug_printf(DBGFLTR_WARN, "(%s:%d) WARNING: " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define TRACE(fmt, ...)  debug_printf(DBGFLTR_TRACE, "(%s:%d) TRACE: " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define INFO(fmt, ...)   debug_printf(DBGFLTR_INFO, "(%s:%d) INFO: " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define DPRINT(fmt, ...) debug_printf(DBGFLTR_DPRINT, "(%s:%d) " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define STATUS(fmt, ...) debug_printf(DBGFLTR_RELEASE, fmt,  ##__VA_ARGS__)
#define DBGFLTR_MISC		0x8
#define DBGFLTR_TRACE		0x7
#define DBGFLTR_INFO		0x6
#define DBGFLTR_DPRINT		0x5
#define DBGFLTR_WARN		0x4
#define DBGFLTR_ERR			0x3
#define DBGFLTR_FATAL		0x2
#define DBGFLTR_RELEASE		0x1

typedef struct _RNA_CODON {
	char* codon;
	char* name;
} RNA_CODON, *PRNA_CODON;

RNA_CODON codon_table[64] = {
	{"UUU", "Phenylalaline"},
	{"UUC", "Phenylalaline"},
	{"UUA", "Leucine"},
	{"UUG", "Leucine"},
	{"UCU", "Serine"},
	{"UCC", "Serine"},
	{"UCA", "Serine"},
	{"UCG", "Serine"},
	{"UAU", "Tyrosine"},
	{"UAC", "Tyrosine"},
	{"UAA", "[Stop]"},
	{"UAG", "[Stop]"},
	{"UGU", "Cysteine"},
	{"UGC", "Cysteine"},
	{"UGA", "[Stop]"},
	{"UGG", "Tryptophan"},
	{"CUU", "Leucine"},
	{"CUC", "Leucine"},
	{"CUA", "Leucine"},
	{"CUG", "Leucine"},
	{"CCU", "Proline"},
	{"CCC", "Proline"},
	{"CCA", "Proline"},
	{"CCG", "Proline"},
	{"CAU", "Histidine"},
	{"CAC", "Histidine"},
	{"CAA", "Glutamine"},
	{"CAG", "Glutamine"},
	{"CGU", "Arginine"},
	{"CGC", "Arginine"},
	{"CGA", "Arginine"},
	{"CGG", "Arginine"},
	{"AUU", "Isoleucine"},
	{"AUC", "Isoleucine"},
	{"AUA", "Isoleucine"},
	{"AUG", "[Methionine (Start)]"},
	{"ACU", "Threonine"},
	{"ACC", "Threonine"},
	{"ACA", "Threonine"},
	{"ACG", "Threonine"},
	{"AAU", "Asparagine"},
	{"AAC", "Asparagine"},
	{"AAA", "Lysine"},
	{"AAG", "Lysine"},
	{"AGU", "Serine"},
	{"AGC", "Serine"},
	{"AGA", "Arginine"},
	{"AGG", "Arginine"},
	{"GUU", "Valine"},
	{"GUC", "Valine"},
	{"GUA", "Valine"},
	{"GUG", "Valine"},
	{"GCU", "Alanine"},
	{"GCC", "Alanine"},
	{"GCA", "Alanine"},
	{"GCG", "Alanine"},
	{"GAU", "Aspartate"},
	{"GAC", "Aspartate"},
	{"GAA", "Glutamate"},
	{"GAG", "Glutamate"},
	{"GGU", "Glycine"},
	{"GGC", "Glycine"},
	{"GGA", "Glycine"},
	{"GGG", "Glycine"},
};

/* global */
static bool generate_amino_acids = false;
static int dna_debug_level = DBGFLTR_MISC;

/*
 * \fn void debug_printf(int dbglevel, char *fmt, ...)
 * \brief Debug printf wrapper
 *
 * \param dbglevel Current debug message level
 * \param fmt Format string
 */
#define BUFSIZE		2048

void debug_printf(int dbglevel, char *fmt, ...) {
	char buf[BUFSIZE];		/* i know,  too large */
	va_list ap;

	va_start(ap, fmt);
	vsnprintf(buf, BUFSIZE, fmt, ap);
	va_end(ap);

	if(dbglevel <= dna_debug_level)
		printf("%s", buf);

	return;
}

#undef BUFSIZE

/*
 * \fn void mrna_to_amino_acid(char* mRNA)
 * \brief Transmute mRNA into chain of amino acids.
 *
 * \param mRNA DNA to be transformed.
 */
void mrna_to_amino_acid(char* mRNA) {
	int i, j, len = strlen(mRNA);
	char tmpbuf[4] = {0};
	
	STATUS("[*] Amino Acid Chain: ");
	
	for(i = 0; i <= (len / 3); i++) {
		memcpy(tmpbuf, mRNA, 3);
		for(j = 0; j <= 63; j++) {
			if(!strcmp(codon_table[j].codon, tmpbuf)) {
				printf("%s ", codon_table[j].name);
			}
		}
		mRNA += 3;
	}
	
	printf("\n");
}
 
bool mrna_shown_warn = false;
/*
 * \fn void dna_to_rna(char* mRNA)
 * \brief Transmute DNA into mRNA
 *
 * \param mRNA DNA to be transformed.
 */
void dna_to_rna(char* mRNA) {
	int i, len = strlen(mRNA);
	for(i = 0; i <= len; i++) {
		if((mRNA[i]) == 'A' || (mRNA[i]) == 'a') {
			mRNA[i] = 'U';
		} else if((mRNA[i]) == 'T' || (mRNA[i]) == 't') {
			mRNA[i] = 'A';
		} else if((mRNA[i]) == 'C' || (mRNA[i]) == 'c') {
			mRNA[i] = 'G';
		} else if((mRNA[i]) == 'G' || (mRNA[i]) == 'g') {
			mRNA[i] = 'C';
		} else if(mRNA[i] != '\0') {
			mRNA[i] = '?';
			if(mrna_shown_warn == false) {
				WARN("Unknown character in tRNA chain, disabling amino lookup.\n");
				mrna_shown_warn = true;
			}
			generate_amino_acids = false;
		}
	}
}

bool trna_shown_warn = false;
/*
 * \fn void rna_to_trna(char* tRNA)
 * \brief Transmute mRNA into tRNA
 *
 * \param tRNA RNA to be transformed.
 */
void rna_to_trna(char* tRNA) {
	int i, len = strlen(tRNA);
	for(i = 0; i <= len; i++) {
		if((tRNA[i]) == 'A' || (tRNA[i]) == 'a') {
			tRNA[i] = 'U';
		} else if((tRNA[i]) == 'U' || (tRNA[i]) == 'u') {
			tRNA[i] = 'A';
		} else if((tRNA[i]) == 'C' || (tRNA[i]) == 'c') {
			tRNA[i] = 'G';
		} else if((tRNA[i]) == 'G' || (tRNA[i]) == 'g') {
			tRNA[i] = 'C';
		} else if(tRNA[i] != '\0') {
			tRNA[i] = '?';
			if(trna_shown_warn == false) {
				WARN("Unknown character in tRNA chain, disabling amino lookup.\n");
				trna_shown_warn = true;
			}
			generate_amino_acids = false;
		}
	}
}

/*
 * \fn int main(int argc, char* argv[])
 * \brief Main function
 *
 * \param argc Argument count
 * \param argv Argument variables
 */
int main(int argc, char* argv[]) {
	char *mRNA = NULL, *tRNA = NULL;

	STATUS("Deoxyribose, version: " __DNA_VERSION__ "\n"
	       "Licensed under the GPLv3, enjoy!\n\n");
	
	/* sanity */
	if(argc != 2) {
		printf("usage: %s [DNA-Sequence]\n", argv[0]);
		exit(-1);
	}
	
	/* verify */
	if(strlen(argv[1]) % 3 == true) {
		WARN("Amino acid chain will not be generated.\n");
		generate_amino_acids = false;
	} else {
		DPRINT("Amino acid chain will be generated.\n");
		generate_amino_acids = true;
	}

	STATUS("[*] DNA Chain: %s\n", argv[1]);
	
	/* mRNA transform */
	mRNA = strdup(argv[1]);
	if(!mRNA) {
		FATAL("Memory for mRNA could not be allocated.\n");
	}
	
	dna_to_rna(mRNA);
	STATUS("[*] mRNA Chain: %s\n", mRNA);

	/* tRNA transform */
	tRNA = strdup(mRNA);
	if(!tRNA) {
		FATAL("Memory for tRNA could not be allocated.\n");
	}
	
	rna_to_trna(tRNA);
	STATUS("[*] tRNA Chain: %s\n", tRNA);
	
	/* amino acid lookup */
	if(generate_amino_acids == true)
		mrna_to_amino_acid(mRNA);
		
	return 0;
}