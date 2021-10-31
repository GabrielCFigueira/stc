#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


unsigned char *randomBitString(int n) {
    srand(clock());
    unsigned char *res = (unsigned char*) malloc(n * sizeof(unsigned char));
    for (int i = 0; i < n; i++)
        res[i] = rand() & 0x1;
    return res;
}

void printBitArray(unsigned char *array, int n) {
    for (int i = 0; i < n; i++)
        printf("%d", array[i]);
    printf("\n");
}

//--------------------------- STC parameters-----------------------------------
int h = 7;
int w = 4;
int H[] = {81, 95, 107, 121};
int H_hat[] = {81, 95, 107, 121};
int Ht[] = {15, 6, 4, 7, 13, 3, 15};


/*
int h = 7;
int w = 2; // 1 / alpha (alpha being the relative payload)
int H[] = {71, 109};
int H_hat[] = {71, 109};
int Ht[] = {3, 2, 3, 1, 0, 1, 3};
*/
/*
int h = 10;
int w = 2;
int H[] = {627, 997};
int H_hat[] = {627, 997};
int Ht[] = {3, 2, 1, 0, 2, 3, 3, 1, 1, 3};
*/

//----------------------------------------------------------------------------



int numBlocks = 200; // message size

int indx = 0;
int indm = 0;

int main() {
    int hpow = 1 << h;


    unsigned char *message = randomBitString(numBlocks);
//    for (int i = 0; i < numBlocks; i++)
//        message[i] = 2;
    unsigned char *cover = randomBitString(numBlocks * w);
    if (message == NULL || cover == NULL)
        printf("OOPS\n");

    float *wght = (float*) malloc(hpow * sizeof(float));
    wght[0] = 0;
    for (int i = 1; i < hpow; i++)
        wght[i] = INFINITY;

    unsigned char *path = malloc(hpow * w * numBlocks * sizeof(unsigned char));
    unsigned char *messagePath = malloc(hpow * numBlocks * sizeof(unsigned char));

    if (path == NULL || messagePath == NULL) {
        printf("ah oh\n");
        exit(1);
    }


    //Forward part of the Viterbi algorithm
    float w0, w1;
    float *newwght = (float*) malloc(hpow * sizeof(float));
    float *temp;

    clock_t start = clock();

    for (int i = 0; i < numBlocks; i++) {

        if (i >= numBlocks - (h - 1)) {
            for (int j = 0; j < w; j++)
                H[j] = H_hat[j] & ((1 << (numBlocks - i)) - 1);
            hpow = hpow >> 1;
        }

        for (int j = 0; j < w; j++) {
            for (int k = 0; k < hpow; k++) {
                w0 = wght[k] + cover[indx];
                w1 = wght[k ^ H[j]] + !cover[indx];

                path[indx * (1 << h) + k] = w1 < w0;
                newwght[k] = w1 < w0 ? w1 : w0;
            }
            indx++;
            temp = wght;
            wght = newwght;
            newwght = temp;

        }

        for (int j = 0; j < hpow >> 1; j++) {
            if(message[indm] == 2) {
                wght[j] = wght[j << 1] < wght[(j << 1) + 1] ? wght[j << 1] : wght[(j << 1) + 1];
                messagePath[indm * (1 << h) + j] = wght[j << 1] < wght[(j << 1) + 1];
            } else {
                wght[j] = wght[(j << 1) + message[indm]];
                messagePath[indm * (1 << h) + j] = message[indm];
            }
        }

        for (int j = hpow >> 1; j < hpow; j++)
            wght[j] = INFINITY;

        indm++;
    }

    //Backward part of the Viterbi algorithm

    float embeddingCost = wght[0];
    int state = 0;
    indx--;
    indm--;
    unsigned char stego[numBlocks * w];


    for (int i = numBlocks - 1; i >= 0; i--) {
        
        message[indm] = messagePath[indm * (1 << h) + state];
        state = (state << 1) + message[indm];
        indm--;
        
        for (int j = w - 1; j >= 0; j--) {
            stego[indx] = path[indx * (1 << h) + state];
            state = state ^ (stego[indx] ? H[j] : 0);
            indx--;
        }
        
        if (i >= numBlocks - (h - 1))
            for (int j = 0; j < w; j++)
                H[j] = H_hat[j] & ((1 << (numBlocks - i + 1)) - 1);

    }

    clock_t end = clock();

    printf("Message:\n");
    printBitArray(message, numBlocks);
    printf("\n");
    printf("Cover:\n");
    printBitArray(cover, numBlocks * w);
    printf("Steganogram:\n");
    printBitArray(stego, numBlocks * w);

    printf("Embedding cost: %d\n", (int) embeddingCost);


    printf("Embedding Time: %lf\n", ((double) end - start) / CLOCKS_PER_SEC);

    unsigned char *ext = (unsigned char*) malloc(numBlocks * sizeof(unsigned char));

    start = clock();

    int line = 0;
    for(int i = 0; i < h; i++)
        line += Ht[i] << (w * i);

    for(int i = 0; i < numBlocks; i++) {
        int mask = 1;
        int index = 0;
        int bit = 0;

        for(int j = w * (i + 1) - 1; j > w * (i + 1 - h) - 1; j--) {

            if (j < 0)
                break;

            bit ^= stego[j] & ((line & mask) >> index);

            index++;
            mask <<= 1;
        }

        ext[i] = bit;
    }

    end = clock();    

    printf("Extraction Time: %lf\n", ((double) end - start) / CLOCKS_PER_SEC);


    return 0;
}
