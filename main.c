#include <stdio.h>
#include <stdlib.h>
#include <intrin.h>

#define MAX_SEQ_LENGTH 100

int max(int a, int b, int c) {
    if (a >= b && a >= c) return a;
    if (b >= a && b >= c) return b;
    return c;
}

void needlemanWunsch(char *seq1, char *seq2, int gapPenalty, int mismatchPenalty) {
    int lenSeq1 = strlen(seq1);
    int lenSeq2 = strlen(seq2);

    int scoreMatrix[MAX_SEQ_LENGTH + 1][MAX_SEQ_LENGTH + 1];

    // Initialize the score matrix
    for (int i = 0; i <= lenSeq1; i++) {
        for (int j = 0; j <= lenSeq2; j++) {
            if (i == 0) {
                scoreMatrix[i][j] = j * gapPenalty;
            } else if (j == 0) {
                scoreMatrix[i][j] = i * gapPenalty;
            } else {
                int match = scoreMatrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? 0 : mismatchPenalty);
                int delete = scoreMatrix[i - 1][j] + gapPenalty;
                int insert = scoreMatrix[i][j - 1] + gapPenalty;

                scoreMatrix[i][j] = max(match, delete, insert);
            }
        }
    }

    // Prints the alignment matrix
    for (int i = 0; i <= lenSeq1; i++) {
        for (int j = 0; j <= lenSeq2; j++) {
            printf("%d\t", scoreMatrix[i][j]);
        }
        printf("\n");
    }


    // Traceback to find the aligned sequences
    int i = lenSeq1;
    int j = lenSeq2;
    int alignedSeq1[MAX_SEQ_LENGTH * 2];  // Assuming the length won't exceed twice the original length
    int alignedSeq2[MAX_SEQ_LENGTH * 2];
    int alignedIndex = 0;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? 0 : mismatchPenalty)) {
            alignedSeq1[alignedIndex] = seq1[i - 1];
            alignedSeq2[alignedIndex] = seq2[j - 1];
            i--;
            j--;
        } else if (i > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gapPenalty) {
            alignedSeq1[alignedIndex] = seq1[i - 1];
            alignedSeq2[alignedIndex] = '-';
            i--;
        } else {
            alignedSeq1[alignedIndex] = '-';
            alignedSeq2[alignedIndex] = seq2[j - 1];
            j--;
        }
        alignedIndex++;
    }

    // Print the aligned sequences in reverse order
    printf("Aligned Sequence 1: ");
    for (int k = alignedIndex - 1; k >= 0; k--) {
        printf("%c", alignedSeq1[k]);
    }
    printf("\n");

    printf("Aligned Sequence 2: ");
    for (int k = alignedIndex - 1; k >= 0; k--) {
        printf("%c", alignedSeq2[k]);
    }
    printf("\n");
}

int main() {
    char seq1[] = "GATCACAG"; // First 8 characters of Human Genome Sequence
    char seq2[] = "TATGCAC"; //Generic 8 character Genomic Sequence

    int gapPenalty = -2;
    int mismatchPenalty = -1;

    needlemanWunsch(seq1, seq2, gapPenalty, mismatchPenalty);

    return 0;
}

