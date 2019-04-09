/******************************************************************
    > File Name: SplitReads.c
    >  Author: yys
    >  mail: shayy0919@163.com
    >  Created Time: 2019年03月20日 星期三 10时44分42秒
******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <getopt.h>

#include "kseq.h"

#define READLEN 150
#define PATH_MAX 128

typedef struct __arg_t {
    int help;
    char read1[PATH_MAX];
    char read2[PATH_MAX];
    char process1[PATH_MAX];
    char outbase[PATH_MAX];
} arg_t;

void Usage(void)
{
    char *usage =
        "\nUsage: SplitReads [options]\n"
        "Version: 1.0\n"
        "\n"
        "Options:\n"
        "       -h|--help        print help infomation\n"
        "       -f|--read1       [required] raw read1 data\n"
        "       -r|--read2       [required] raw read2 data\n"
        "       -p|--process1    [required] processed read1 data include trimmed and untrimmed\n"
        "       -o|--outbase     [required] trimmed & untrimmed reads for fastq file basename [.fq|.gz]\n";

    fprintf(stderr, "%s", usage);
    exit(-1);
}

static const struct option long_options[] =
{
    { "help", no_argument, NULL, 'h' },
    { "read1", required_argument, NULL, 'f' },
    { "read2", required_argument, NULL, 'r' },
    { "process1", required_argument, NULL, 'p' },
    { "outbase", required_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

arg_t *ParseOpt( int argc, char **argv )
{
    int opt =0, opterr =0;
    arg_t *Arg;

    Arg = (arg_t *)calloc(1, sizeof(arg_t));
    while ( (opt = getopt_long(argc, argv, "f:r:p:o:h", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': Arg->help = 1; break;
            case 'f': strcpy(Arg->read1, optarg); break;
            case 'r': strcpy(Arg->read2, optarg); break;
            case 'p': strcpy(Arg->process1, optarg); break;
            case 'o': strcpy(Arg->outbase, optarg); break;
            case '?': fprintf(stderr, \
                            "[Err::%s::%d]  Option error occour!.\n", __func__, __LINE__);
                      Arg->help = 1;
        }
    }
    if (!Arg->read1[0] || !Arg->read2[0] || !Arg->process1[0] || !Arg->outbase[1]) {
        fprintf(stderr, \
            "[Err::%s::%d]  Please give the [requied] parmeters!\n", __func__, __LINE__);
        Arg->help = 1;
    }

    return Arg;
}

KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv)
{
    arg_t *args = ParseOpt(argc, argv);
    if (args->help) Usage();

    FILE *fp_trimmed_R1, *fp_trimmed_R2;
    FILE *fp_untrimmed_R1, *fp_untrimmed_R2;
    gzFile fp1_1, fp1_2, fp2_1;
    kseq_t *seq1_1, *seq1_2, *seq2_1;
    char trimmed_R1[PATH_MAX], trimmed_R2[PATH_MAX];
    char untrimmed_R1[PATH_MAX], untrimmed_R2[PATH_MAX];

    sprintf(trimmed_R1, "%s%s", args->outbase, "_trimmed_R1.fq");
    sprintf(trimmed_R2, "%s%s", args->outbase, "_trimmed_R2.fq");
    sprintf(untrimmed_R1, "%s%s", args->outbase, "_untrimmed_R1.fq");
    sprintf(untrimmed_R2, "%s%s", args->outbase, "_untrimmed_R2.fq");

    fp1_1 = gzopen(args->read1, "r");
    fp1_2 = gzopen(args->read2, "r");
    fp2_1 = gzopen(args->process1, "r");

    fp_trimmed_R1 = fopen(trimmed_R1, "w");
    fp_trimmed_R2 = fopen(trimmed_R2, "w");
    fp_untrimmed_R1 = fopen(untrimmed_R1, "w");
    fp_untrimmed_R2 = fopen(untrimmed_R2, "w");

    seq1_1 = kseq_init(fp1_1);
    seq1_2 = kseq_init(fp1_2);
    seq2_1 = kseq_init(fp2_1);

    while((kseq_read(seq1_1) >= 0 && kseq_read(seq1_2) >= 0) && kseq_read(seq2_1)) {
        if (strlen(seq1_1->seq.s) != strlen(seq2_1->seq.s)) {
            fprintf(fp_trimmed_R1, "@%s %s\n%s\n+\n%s\n", \
                    seq1_1->name.s, seq1_1->comment.s, seq1_1->seq.s, seq1_1->qual.s);
            fprintf(fp_trimmed_R2, "@%s %s\n%s\n+\n%s\n", \
                    seq1_2->name.s, seq1_2->comment.s, seq1_2->seq.s, seq1_2->qual.s);
        }
        else {
            fprintf(fp_untrimmed_R1, "@%s %s\n%s\n+\n%s\n", \
                    seq1_1->name.s, seq1_1->comment.s, seq1_1->seq.s, seq1_1->qual.s);
            fprintf(fp_untrimmed_R2, "@%s %s\n%s\n+\n%s\n", \
                    seq1_2->name.s, seq1_2->comment.s, seq1_2->seq.s, seq1_2->qual.s);
        }
    }

    return 0;
}
