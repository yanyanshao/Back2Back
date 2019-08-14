/******************************************************************
    > File Name: Cutadapt.B2B.c
    >  Author: yys
    >  mail: shayy0919@163.com
    >  Created Time: 2019年03月07日 星期四 15时35分24秒
******************************************************************/

#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

#include "kseq.h"
#include "parse.h"

#define KEYLEN 128 // maxmum length for hash key
#define BARCODE 6
#define LOCNUM 8   // initial number for location index
#define READLEN 150
#define ISHEAD(_node) (!_node->next && !_node->num)

KSEQ_INIT(gzFile, gzread)

#define err_open(_fp, _fn, _mode) do {\
    _fp = fopen(_fn, _mode); \
    if (!_fp) { \
        fprintf(stderr, "\nerr: failed to open %s!\n", _fn); exit(-1); \
    } \
} while(0)

/*! @typedef loc_t
 @abstract structure for the location index
 @field frloc   forward or reverse index [0=>forward and 1=>reverse]
 @field sloc    int index of key in the index sequence
 */
typedef struct __loc_t {
    int frloc;
    int sloc;
}loc_t;

typedef struct __node_t node_t;
struct __node_t {
    int num;   
    char key[KEYLEN];
    loc_t *loc;
    node_t *next;
};

typedef struct __hash_t {
    uint32_t size;
    node_t *table;
} hash_t;

typedef struct __status {
    int find;
    node_t *node;
} status;

typedef struct __query_t {
    int qstart;
    int qend;
    int fr;
} query_t;

static int ProperSize(int size)
{
    int flag = 1;
    int tril = size;

    while(1) {
        tril++; flag = 1;
        for (int i=2; i<tril; i++){
            if (tril % i == 0) {
                flag = 0; break;
            }
        }
        if (flag)
            return tril;
    }
}

void PrintHash(hash_t *H)
{
    node_t *node;
    node_t *p;
    int num=0;

    for (int i=0; i<H->size; i++){
        node = (node_t *)malloc(sizeof(H->table[i]));
        node = &H->table[i];
        num = node->num;
        do{
            for (int j=0; j<num; j++)
                printf("%s\t%d\t%d\n",node->key,node->loc[j].frloc,node->loc[j].sloc);
        }while((node = node->next) != NULL);
    }
}

hash_t *InitHash(int size)
{
    hash_t *H;

    H = (hash_t *)malloc(sizeof(hash_t));
    if (!H)
        goto _memerror;
    
    H->size = ProperSize(size);
    H->table = (node_t *)calloc(H->size, sizeof(node_t));
    if (!H->table)
        goto _memerror;

    return H;
    
    _memerror:
        fprintf(stderr, "[Err::%s] \
                Failed to allocate memory\n", __func__); exit(-1);

}

int Hash(const char *key, int size)
{
    uint32_t val = 0;

    while (*key != '\0')
        val = (val << 5) -val + *key++;

    return val % size;
}

void Search(char *key, hash_t *H, status *s)
{
    node_t *p, *t;

    p = &H->table[Hash(key, H->size)];

    if(ISHEAD(p)) {
        s->find = 0; s->node = p;
        return ;
    }    
    do {
        t = p;
        if (!strcmp(key, p->key)) {
            s->find = 1; s->node = p;
            return ;
        } 
        p = p->next;
    } while(p);
    s->find = 0; s->node = t;
}

static node_t *NewNode(void)
{
    node_t *node;

    node = (node_t *)calloc(1, sizeof(node_t));
    if (!node)
        goto _memerror;
    node->loc = (loc_t *)calloc(1, sizeof(loc_t));
    if (!node->loc)
        goto _memerror;

    return node;

  _memerror:
      fprintf(stderr, "[Err::%s] \
          Failed to allocate memory\n", __func__); exit(-1);

}

void UpdateLoc(node_t *node, loc_t *loc)
{
    int locnum;

    locnum = node->num;
    if (locnum % LOCNUM == 0) {
        int newsize = locnum + LOCNUM;
        loc_t *tem = (loc_t *)realloc(node->loc, newsize*sizeof(loc_t));
        if (!tem)
            goto _memerror;
        node->loc = tem;
    }
    node->loc[locnum].frloc = loc->frloc;
    node->loc[locnum].sloc = loc->sloc;

    node->num++;

    return ;

  _memerror:
      fprintf(stderr, "[Err:%s] Failed to alloc memory\n", __func__);
}

void Insert(char *key, loc_t *loc, hash_t *H)
{
    status S;

    Search(key, H, &S);

    if (!S.find && ISHEAD(S.node)) {
        strcpy(S.node->key, key);
        UpdateLoc(S.node, loc);return ;
    }
    if (!S.find) {
        S.node->next = NewNode();
        strcpy(S.node->next->key, key);
        UpdateLoc(S.node->next, loc); 
        return ;
    }
    if (S.find) {
        UpdateLoc(S.node, loc); return ;
    }
}

void Index(char (*indexlist)[50], hash_t *T, int kmer)
{
    char slic[KEYLEN];
    loc_t loc;
    int ilen;

    for (int i=0; i<2; i++) {
        ilen = strlen(*(indexlist+i)) - kmer + 1;
        for (int k=0; k<ilen; k++) {
            strncpy(slic, *(indexlist+i) + k, kmer); 
            slic[kmer] = '\0';
            loc.frloc = i; loc.sloc = k;
            Insert(slic, &loc, T);  
        }
    }
    return ;
}

bool MisCheck(char *string1, char *string2, int mismatch)
{
    int len1 = strlen(string1);
    int len2 = strlen(string2);
    int minlen, misnum = 0;

    if (len1 < len2) {minlen = len1;} else {minlen = len2;}
    for (int i=0; i<minlen; i++){
        if (*string1++ != *string2++) {
            misnum++;

        if (misnum > mismatch) return false;
        }
    }
    return true;
}

static query_t SeqQuery(char *seq, hash_t *H, char (*indexlist)[50], int mis, int kmer)
{
    int seql;
    int num;
    int fr, sloc;
    char key[KEYLEN];
    status S;
    query_t Q = {0};

    seql = strlen(seq);
    for(int i=0; i<seql; i++){
        if (*(seq+i) == 'N') continue;
        
        strncpy(key,(seq+i),kmer); key[kmer] = '\0';
        Search(key, H, &S);

        if (S.find == 0) continue;

        num = S.node->num;
        for(int j=0; j<num; j++){
            fr = S.node->loc[j].frloc;
            sloc = S.node->loc[j].sloc;
            if ((i-sloc) < 0) break;
            if (MisCheck(seq+i-sloc,*(indexlist+fr),mis)) {
                Q.qstart = i-sloc;
                Q.qend = i-sloc + strlen(*(indexlist+fr));
                Q.fr = fr;
                return Q;  
            }

        }
    }
    return Q;
}

int CutBarcodeSeq(query_t *Q, kseq_t *seq, char (*indexlist)[50])
{
    char seqs1[READLEN] = "\0";
    char qual1[READLEN] = "\0";
    char barcode[BARCODE]= "\0";

    if (Q->fr == 0) {
        strncpy(barcode,(seq->seq.s)+Q->qend, BARCODE);
        //barcode[BARCODE] = '\0';
        strcpy(seqs1, (seq->seq.s)+Q->qend+BARCODE+strlen(*(indexlist+2))+strlen(*(indexlist+1)));
        strncat(seqs1,(seq->seq.s), Q->qstart);
        strcpy(seq->seq.s, seqs1);

        strcpy(qual1, (seq->qual.s)+Q->qend+BARCODE+strlen(*(indexlist+2))+strlen(*(indexlist+1)));
        strncat(qual1,(seq->qual.s), Q->qstart);
        strcpy(seq->qual.s, qual1);

        strcat(seq->name.s,"|");
        strcat(seq->name.s,barcode);
        return 0;
    }
    if (Q->fr == 1) {
        strncpy(seqs1, (seq->seq.s), Q->qstart);
        strcat(seqs1,(seq->seq.s)+Q->qend+strlen(*(indexlist+2))+strlen(*(indexlist+0)));
        strcpy(seq->seq.s, seqs1);

        strncpy(qual1, (seq->qual.s), Q->qstart);
        strcat(qual1,(seq->qual.s)+Q->qend+strlen(*(indexlist+2))+strlen(*(indexlist+0)));
        strcpy(seq->qual.s, qual1);
        return 1;
    }

}

void SeqWrite(kseq_t *seq1, kseq_t *seq2, FILE *fp1, FILE *fp2)
{
    fprintf(fp1, "%s %s\n%s\n+\n%s\n", \
            seq1->name.s, seq1->comment.s, seq1->seq.s, seq1->qual.s);
    fprintf(fp2, "%s %s\n%s\n+\n%s\n", \
            seq2->name.s, seq2->comment.s, seq2->seq.s, seq2->qual.s);
}

int main(int argc, char **argv)
{
    arg_t *args = ParseOpt(argc, argv);
    if (args->help) Usage();


    FILE *fqo1, *fqo2;
    char fname1[PATH_MAX];
    char fname2[PATH_MAX];
    char *bname1 = "_trimmed_R1.fq";
    char *bname2 = "_trimmed_R2.fq";
    sprintf(fname1, "%s%s", args->outdir, bname1);
    sprintf(fname2, "%s%s", args->outdir, bname2);
    err_open(fqo1, fname1, "w");
    err_open(fqo2, fname2, "w");
        
    char *f1;
    gzFile fp1, fp2;
    kseq_t *seq1, *seq2;
    char seqs1[READLEN] = "\0";
    char seqs2[READLEN] = "\0";
    char barcode[BARCODE] = "\0";
    int flag1, flag2;
    char FR[3][50] = {{0}};
    query_t Q1, Q2;

    hash_t *H;
    
    /*F : GCTCGGAGATGTGTATAAGAGACAGNNNNNNATTGGAGTCCT
     *    -------------------------      -----------
     *              FR[0]                    FR[3]
     *R : GGACTCCAATACACTCTATCGCTACACGACG
     *    -------------------------------
     *              FR[2]
    */
    f1 = strchr(args->adapt1, 'N');
    strncpy(FR[0], args->adapt1, (f1-args->adapt1)/sizeof(char));
    strcpy(FR[2], args->adapt1+strlen(FR[0])+6*sizeof(char));
    strcpy(FR[1], args->adapt2);
    
    H = InitHash(2 << 5);
    Index(FR, H, args->kmer);
    
    
    fp1 = gzopen(args->read1, "r");
    fp2 = gzopen(args->read2, "r");
    
    if (fp1 == NULL || fp2 == NULL){
        fprintf(stderr, "Fail to open file %s or %s!\n", __func__, args->read1, args->read2);
        return -1;
    }

    seq1 = kseq_init(fp1);
    seq2 = kseq_init(fp2);

    while((kseq_read(seq1) >=0) && (kseq_read(seq2) >= 0)){
        Q1 = SeqQuery(seq1->seq.s, H, FR, args->mismatch, args->kmer);
        Q2 = SeqQuery(seq2->seq.s, H, FR, args->mismatch, args->kmer);

        if (Q1.qend !=0 || Q2.qend !=0) {
            flag1 = CutBarcodeSeq(&Q1, seq1, FR);       
            flag2 = CutBarcodeSeq(&Q2, seq2, FR);

            if ((flag1 & flag2) == 0) {
                if(strcmp(seq1->name.s,seq2->name.s)>0)
                    strcpy(seq2->name.s, seq1->name.s);
                if(strcmp(seq1->name.s,seq2->name.s)<0)
                    strcpy(seq1->name.s, seq2->name.s);

                SeqWrite(seq1, seq2, fqo1, fqo2);
                continue;

                printf("name:%s\n",seq1->name.s);
                printf("comment:%s\n", seq1->comment.s);
                printf("sequence:%s\n", seq1->seq.s);
                printf("quality:%s\n", seq1->qual.s); 

                printf("name:%s\n",seq2->name.s);
                printf("comment:%s\n", seq2->comment.s);
                printf("sequence:%s\n", seq2->seq.s);
                printf("quality:%s\n", seq2->qual.s);
            }
            SeqWrite(seq1, seq2, fqo1, fqo2); continue;
        }
        SeqWrite(seq1, seq2, fqo1, fqo2); continue;
    }
    
    
    /*
    while((kseq_read(seq1) >= 0)){
        printf("name:%s\t%d\t%d\n",seq1->name.s,seq1->name.l,seq1->name.m);
        printf("comment:%s\t%d\t%d\n", seq1->comment.s,seq1->comment.l,seq1->comment.m);
        printf("sequence:%s\t%d\t%d\n", seq1->seq.s,seq1->seq.l,seq1->seq.m);
        printf("quality:%s\t%d\t%d\n", seq1->qual.s,seq1->qual.l,seq1->qual.m);
       

    }
    */
    
}
    
