#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
unsigned long long SEED = 31122891;
unsigned long long RANV;
int RANI = 0;
double Ranq1();
void normal(double sig, double *n1, double *n2);
double CHK(double L1, double L2);
double table[8] = {0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, 0};

int main() {
    int i, j, k1, m; 
    int n, rc;              // n is column 9872 and rc is row 4936
    int dv1, dv2, dc1, dc2, dc3;              // dv: column have #1 and dc: row have #1
    int **L1 = NULL;                // check node connect 6 variable nodes
    int L1lenrow = rc/2;        // rc = 4936/2
    int L1lencolumn = dc1;     // dc = 12
    int **L2 = NULL; //y                // check node connect 6 variable nodes
    int L2lenrow = rc/2-617;        // rc = 4936/2-617=1851
    int L2lencolumn = dc2;     // dc = 5
    int **L3 = NULL; //y                // check node connect 6 variable nodes
    int L3lenrow = 617;        // rc = 617
    int L3lencolumn = dc3;     // dc = 6
    int **M1 = NULL;//y
    int M1lenrow = 5553;       // n = 9872/2
    int M1lencolumn = dv1;   // dv1 = 3
    int **M2 = NULL;//y
    int M2lenrow = 4319;       // n = 9872/2
    int M2lencolumn = dv2;   // dv2 = 6
    double **qij1 = NULL;    // from down to top
    int qij1row = dv1;        // qijrow = dv1 = 3
    int qij1column = 5553;      // qijcolumn = n = 7376/2
    double **qij2 = NULL;    // from down to top
    int qij2row = dv2;        // qijrow = dv2 = 6
    int qij2column = 4319;      // qijcolumn = n = 7376/2
    double **uij1 = NULL;    // from top to down
    int uij1row = rc/2;        // uijrow = rc = 3688/2
    int uij1column = dc1;     // uijcloumn = dc = 12
    double **uij2 = NULL;    // from top to down
    int uij2row = rc/2-617;        // uijrow = rc = 3688
    int uij2column = dc2;     // uijcloumn = dc = 5
    double **uij3 = NULL;    // from top to down
    int uij3row = 617;        // uijrow = rc = 3688
    int uij3column = dc3;     // uijcloumn = dc = 6

    int *codarray = NULL;          //codeword, Dynamic memory allocation, length = codarraylen
    int codarraylen = n;
    int *codearray = NULL;         //0->1; 1->-1,  Dynamic memory allocation, length = codearraylen
    int codearraylen = n;
    double *outp = NULL;           // codeword + noise, Dynamic memory allocation, length = outparray
    int outparray = n;
    int *output = NULL;            // result of interative algorithm decoding, Dynamic memory allocation, length = outputarray
    int outputarray = n;
    double *Lj = NULL;             // LLR
    int Ljlen = n;
    double tempqij1[11];
    double tempqij2[4]; 
    double tempqij3[5]; 
    double tempuij;
    double temp1uij1[3];
    double temp1uij2[6];
    double temp1qij; 
    int *comput = NULL;
    int computlen = n;      // computlen = 816
    int *comput1 = NULL;
    int comput1len = rc;    // comput1len = 408
    int *comput2 = NULL;
    int comput2len = rc;
    double *qj = NULL;
    int qjlen = n;
    int *checkbit = NULL;
    int checkbitlen = rc;
    int krc;
    int **G = NULL;
    int Glenrow = krc;
    int Glencolumn = n;
    int *u = NULL;
    int ulen = krc;

    int step;               // test 6 EbN0(dB) result
    int s = 0;              // receive 100 error block
    int num = 0;            // do compute block
    int error1;
    int error2;
    int totalerror1=0;
    int totalerror2=0;
    int restart = 0;
    double sigma;
    double ebn0;
    double ebn0s[6];
    double bers1[6];
    double bers2[6];
    //double berscompare[6];
    double x, y;
    int valL;
    int valL2;
    double app;
    double app1;
    int stp;
    for (i = 0; i < 6; i++) bers1[i] = 0;  
    for (i = 0; i < 6; i++) bers2[i] = 0;   
    ebn0s[0] = 0.2;
    ebn0s[1] = 0.6;
    ebn0s[2] = 0.8;
    ebn0s[3] = 1.0;
    ebn0s[4] = 1.2;
    // open file
    FILE *fpr;
    fpr=fopen("C3parchematrix.txt","r");
    fscanf(fpr,"%d",&n);
    fscanf(fpr,"%d",&rc);
    printf("column = %d\n", n);
    printf("row = %d\n", rc);
    fscanf(fpr,"%d",&dc1);
    fscanf(fpr,"%d",&dc2);
    fscanf(fpr,"%d",&dc3);
    fscanf(fpr,"%d",&dv1);
    fscanf(fpr,"%d",&dv2);
    printf("dc = %d\n", dc1); 
    printf("dc = %d\n", dc2);
    printf("dc = %d\n", dc3);
    printf("dv = %d\n", dv1);
    printf("dv = %d\n", dv2);
    L1lenrow = rc/2;
    L1lencolumn = dc1;
    L2lenrow = 1851;
    L2lencolumn = dc2;
    L3lenrow = 617;
    L3lencolumn = dc3;

    M1lenrow = 5553; 
    M1lencolumn = dv1;
    M2lenrow = 4319; 
    M2lencolumn = dv2;

    M1 = (int **)malloc(M1lenrow * sizeof(int *));
    for (i = 0; i < M1lenrow; i++) M1[i] = (int *)malloc(M1lencolumn * sizeof(int));
    M2 = (int **)malloc(M1lenrow * sizeof(int *));
    for (i = 0; i < M2lenrow; i++) M2[i] = (int *)malloc(M2lencolumn * sizeof(int));
    L1 = (int **)malloc(L1lenrow * sizeof(int *));
    for (i = 0; i < L1lenrow; i++) L1[i] = (int *)malloc(L1lencolumn * sizeof(int));
    L2 = (int **)malloc(L2lenrow * sizeof(int *));
    for (i = 0; i < L2lenrow; i++) L2[i] = (int *)malloc(L2lencolumn * sizeof(int));
    L3 = (int **)malloc(L3lenrow * sizeof(int *));
    for (i = 0; i < L3lenrow; i++) L3[i] = (int *)malloc(L3lencolumn * sizeof(int));

    for (j = 0; j < M1lenrow; j++) {
        for (i = 0; i < M1lencolumn; i++) {
            fscanf(fpr,"%d",&M1[j][i]);
        }
    }
    for (j = 0; j < M2lenrow; j++) {
        for (i = 0; i < M2lencolumn; i++) {
            fscanf(fpr,"%d",&M2[j][i]);
        }
    }
    int temp;
    for (j = 0; j < M1lenrow; j++) {
        for (i = 0; i < M1lencolumn; i++) {
            for (m = i; m < M1lencolumn; m++) {
                if (M1[j][m] < M1[j][i]) {
                    temp = M1[j][m];
                    M1[j][m] = M1[j][i];
                    M1[j][i] = temp;
                }
            }
        }
    }
    for (j = 0; j < M2lenrow; j++) {
        for (i = 0; i < M2lencolumn; i++) {
            for (m = i; m < M2lencolumn; m++) {
                if (M2[j][m] < M2[j][i]) {
                    temp = M2[j][m];
                    M2[j][m] = M2[j][i];
                    M2[j][i] = temp;
                }
            }
        }
    }
    for (i = 0; i < L1lenrow; i++) {
        for (j = 0; j < L1lencolumn; j++) {
            fscanf(fpr,"%d",&L1[i][j]);
        }
    }
    for (i = 0; i < L2lenrow; i++) {
        for (j = 0; j < L2lencolumn; j++) {
            fscanf(fpr,"%d",&L2[i][j]);
        }
    }
    for (i = 0; i < L3lenrow; i++) {
        for (j = 0; j < L3lencolumn; j++) {
            fscanf(fpr,"%d",&L3[i][j]);
        }
    }
    for (i = 0; i < L1lenrow; i++) {
        for (j = 0; j < L1lencolumn; j++) {
            for (m = j; m < L1lencolumn; m++) {
                if (L1[i][m] < L1[i][j]) {
                    temp = L1[i][m];
                    L1[i][m] =L1[i][j];
                    L1[i][j] = temp;
                }
            }
        }
    }
    for (i = 0; i < L2lenrow; i++) {
        for (j = 0; j < L2lencolumn; j++) {
            for (m = j; m < L2lencolumn; m++) {
                if (L2[i][m] < L2[i][j]) {
                    temp = L2[i][m];
                    L2[i][m] = L2[i][j];
                    L2[i][j] = temp;
                }
            }
        }
    }
    for (i = 0; i < L3lenrow; i++) {
        for (j = 0; j < L3lencolumn; j++) {
            for (m = j; m < L3lencolumn; m++) {
                if (L3[i][m] < L3[i][j]) {
                    temp = L3[i][m];
                    L3[i][m] = L3[i][j];
                    L3[i][j] = temp;
                }
            }
        }
    }
    fclose(fpr);
    FILE *fpr1;
    fpr1=fopen("generator1.txt","r");
    //fscanf(fpr1,"%d",&krc);
    //fscanf(fpr1,"%d",&n);
    n = 9872;
    krc = 4936;
    printf("column = %d\n", n);
    printf("row = %d\n", krc);
    Glenrow = krc;
    Glencolumn = n;
    G = (int **)malloc(Glenrow * sizeof(int *));
    for (i = 0; i < Glenrow; i++) G[i] = (int *)malloc(Glencolumn * sizeof(int));
    for (i = 0; i < Glenrow; i++) {
        for (j = 0; j < Glencolumn; j++) {
            fscanf(fpr1,"%d",&G[i][j]);
            //if (G[i][j] == 1) printf("%d ", j);
        }
        //printf("\n");
    }
    // for CODE part
    fclose(fpr1);
    codarraylen = n;
    codarray = (int *)malloc(codarraylen * sizeof(int));
    if( codarray == NULL ) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    codearraylen = n;
    codearray = (int *)malloc( codearraylen * sizeof(int) );
    if( codearray == NULL ) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    outparray = n;
    outp = (double *)malloc( outparray * sizeof(double) );
    if (outp == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    outputarray = n;
    output = (int *)malloc( outputarray * sizeof(int) );
    if (output == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    Ljlen = n;
    Lj = (double *)malloc(Ljlen * sizeof(double));
    if (Lj == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    qij1row = dv1;        // qijrow = dv = 3
    qij1column = 5553;      // qijcolumn = n = 816
    qij1 = (double **)malloc(qij1row * sizeof(double *));
    for (i = 0; i < qij1row; i++) qij1[i] = (double *)malloc(qij1column * sizeof(double));
    qij2row = dv2;        // qijrow = dv = 3
    qij2column = 4319;      // qijcolumn = n = 816
    qij2 = (double **)malloc(qij2row * sizeof(double *));
    for (i = 0; i < qij2row; i++) qij2[i] = (double *)malloc(qij2column * sizeof(double));
    
    uij1row = rc/2;        // uijrow = rc = 408
    uij1column = dc1;     // uijcolumn = dc = 6 
    uij1 = (double **)malloc(uij1row * sizeof(double *));
    for (i = 0; i < uij1row; i++) uij1[i] = (double *)malloc(uij1column * sizeof(double));
    uij2row = 1851;        // uijrow = rc = 408
    uij2column = dc2;     // uijcolumn = dc = 6 
    uij2 = (double **)malloc(uij2row * sizeof(double *));
    for (i = 0; i < uij2row; i++) uij2[i] = (double *)malloc(uij2column * sizeof(double));
    uij3row = 617;        // uijrow = rc = 408
    uij3column = dc3;     // uijcolumn = dc = 6 
    uij3 = (double **)malloc(uij3row * sizeof(double *));
    for (i = 0; i < uij3row; i++) uij3[i] = (double *)malloc(uij3column * sizeof(double));
    computlen = n;
    comput = (int *)malloc(computlen * sizeof(int));
    comput1len = rc;
    comput1 = (int *)malloc(comput1len * sizeof(int));
    comput2len = rc;
    comput2 = (int *)malloc(comput2len * sizeof(int));
    qjlen = n;
    qj = (double *)malloc(qjlen * sizeof(double));
    checkbitlen = rc;
    checkbit = (int *)malloc(checkbitlen * sizeof(int));
    ulen = krc;
    u = (int *)malloc(ulen * sizeof(int));
    for (step = 0; step < 1; step++) {
        s = 0;
        num = 0;
        totalerror1 = 0;
        totalerror2 = 0;
        while (s < 100) {
            for (i = 0; i < codarraylen; i++) {
                codarray[i] = 0;
            }
            if (num == 0) {
                u[0] = 1;
                u[1] = 0;
                u[2] = 0;
                u[3] = 0;
                u[4] = 0;
                u[5] = 0;
                for (i = 0; i < krc - 6; i++) u[i + 6] = (u[i + 1] + u[i]) % 2;
            }
            else {
                u[0] = (u[krc-5] + u[krc-6]) % 2;
                u[1] = (u[krc-4] + u[krc-5]) % 2;
                u[2] = (u[krc-3] + u[krc-4]) % 2;
                u[3] = (u[krc-2] + u[krc-3]) % 2;
                u[4] = (u[krc-1] + u[krc-2]) % 2;
                u[5] = (u[krc-0] + u[krc-1]) % 2;
                for (i = 0; i < krc - 6; i++) u[i + 6] = (u[i + 1] + u[i]) % 2;
            }
            num++;
            for (i = 0; i < ulen; i++) {
                if (u[i] == 1) {
                    for (j = 0; j < n; j++) codarray[j] = (codarray[j] + G[i][j]) % 2;
                }
            }
            free(u);
            for (i = 0; i < codearraylen; i++) {
                if (codarray[i] == 0) codearray[i] = 1;
                else codearray[i] = -1;
            }
            ebn0 = ebn0s[step];
            sigma = sqrt(1.0 / (pow(10, ebn0/10)));
            for(i = 0; i < rc; i++) {
                normal(sigma, &x, &y);
                outp[2 * i] = codearray[2 * i] + x;
                outp[2 * i + 1] = codearray[2 * i + 1] + y;
            }
            free(codearray);
            ebn0 = pow(10, ebn0/10);
            for(i = 0; i < Ljlen; i++) {
                Lj[i] =4 * 0.5 * ebn0 * outp[i];     //  0.5 * 1.2544 = Es/N0
            }
            free(outp);
            for (j = 0; j < qij1column; j++) {                               // initialization
                for (i = 0; i < qij1row; i++) {
                        qij1[i][j] = Lj[j];
                }  
            }
            for (j = 0; j < qij2column; j++) {                               // initialization
                for (i = 0; i < qij2row; i++) {
                        qij2[i][j] = Lj[5553 + j];
                }  
            }
            
            for (k1 = 0; k1 < 100 && restart != rc; k1++) {         // message passing, for predetermined threshold = 100
                restart = 0;  
                printf("yes");  
                for (i = 0; i < 11; i++) {                          // bottom-up
                    tempqij1[i] = 0.0;
                }
                printf("1 ");
                for (i = 0; i < 4; i++) {                          // bottom-up
                    tempqij2[i] = 0.0;
                }
                printf("1 ");
                for (i = 0; i < 5; i++) {                          // bottom-up
                    tempqij3[i] = 0.0;
                }
                for (i = 0; i < computlen; i++) {
                    comput[i] = 0;
                }
                printf("1 ");
                for (i = 0; i < rc; i++) {
                    if (i < (rc/2)) {
                        for (j = 0; j < 12; j++) {
                            for (m = 0; m < 11; m++) {
                                if (m < j) {
                                    valL = L1[i][m]-1;
                                    if (valL < 5553) tempqij1[m] = qij1[comput[valL]][valL];
                                    else tempqij1[m] = qij2[comput[valL]][valL - 5553];
                                } else if (m >= j) {
                                    valL = L1[i][m+1]-1;
                                    if (valL < 5553) tempqij1[m] = qij1[comput[valL]][valL];
                                    else tempqij1[m] = qij2[comput[valL]][valL - 5553];
                                }
                            }
                            tempuij = tempqij1[0];
                            for(m = 1; m < 11; m++) {
                                tempuij = CHK(tempuij, tempqij1[m]);
                            }
                            uij1[i][j] = tempuij;
                        }
                    } else if (i >=(rc/2) && i < 4319) {
                        for (j = 0; j < 5; j++) {
                            for (m = 0; m < 4; m++) {
                                if (m < j) {
                                    valL = L2[i-2468][m]-1;
                                    if (valL < 5553) tempqij2[m] = qij1[comput[valL]][valL];
                                    else tempqij2[m] = qij2[comput[valL]][valL-5553];
                                } 
                                else if (m >= j) {
                                    valL = L2[i-2468][m+1]-1;
                                    if (valL < 5553) tempqij2[m] = qij1[comput[valL]][valL];
                                    else tempqij2[m] = qij2[comput[valL]][valL-5553];
                                }
                            }
                            tempuij = tempqij2[0];
                            for(m = 1; m < 4; m++) {
                                tempuij = CHK(tempuij, tempqij2[m]);
                            }
                            uij2[i-2468][j] = tempuij;
                        }
                    } else {
                        for (j = 0; j < 6; j++) {
                            for (m = 0; m < 5; m++) {
                                if (m < j) {
                                    valL = L3[i-4319][m]-1;
                                    if (valL < 5553) tempqij3[m] = qij1[comput[valL]][valL];
                                    else tempqij3[m] = qij2[comput[valL]][valL-5553];
                                } 
                                else if (m >= j) {
                                    valL = L3[i-4319][m+1]-1;
                                    if (valL < 5553) tempqij3[m] = qij1[comput[valL]][valL];
                                    else tempqij3[m] = qij2[comput[valL]][valL-5553];
                                }
                            }
                            tempuij = tempqij3[0];
                            for(m = 1; m < 5; m++) {
                                tempuij = CHK(tempuij, tempqij3[m]);
                            }
                            uij3[i-4319][j] = tempuij;
                        }
                    }
                    if (i < (rc/2)) {
                        for (m = 0; m < 12; m++) {
                            comput[L1[i][m] - 1] += 1;
                        }
                    } else if (i >=(rc/2) && i < 4319) {
                        for (m = 0; m < 5; m++) {
                            comput[L2[i-2468][m] - 1] += 1;
                        }
                    } else {
                        for (m = 0; m < 6; m++) {
                            comput[L3[i-4319][m] - 1] += 1;
                        }
                    }
                }

                // top-down
                for(i = 0; i < 3; i++) {
                    temp1uij1[i] = 0.0;
                }
                for(i = 0; i < 6; i++) {
                    temp1uij2[i] = 0.0;
                }
                for (i = 0; i < comput1len; i++) comput1[i] = 0;
                for (j = 0; j < n; j++) {
                    if (j < 5553) {
                        for (i = 0; i < 3; i++) {
                            for (m = 0; m < 2; m++) {
                                if (m < i) { 
                                    valL = M1[j][m] - 1;
                                    if (valL < (rc/2)) temp1uij1[m] = uij1[valL][comput1[valL]];
                                    else if (valL >= (rc/2) && valL < 4319) temp1uij1[m] = uij2[valL-2468][comput1[valL]];
                                    else temp1uij1[m] = uij3[valL-4319][comput1[valL]];
                                }
                                else if (m >= i) {
                                    valL = M1[j][m + 1] - 1;
                                    if (valL < (rc/2)) temp1uij1[m] = uij1[valL][comput1[valL]];
                                    else if (valL >= (rc/2) && valL < 4319) temp1uij1[m] = uij2[valL-2468][comput1[valL]];
                                    else temp1uij1[m] = uij3[valL-4319][comput1[valL]];
                                }
                            }
                            temp1uij1[2] = Lj[j];
                            qij1[i][j] = temp1uij1[0] + temp1uij1[1] + temp1uij1[2];
                        }
                    } else {
                        for (i = 0; i < 6; i++) {
                            for (m = 0; m < 5; m++) {
                                if (m < i) { 
                                    valL = M2[j-5553][m] - 1;
                                    if (valL < (rc/2)) temp1uij2[m] = uij1[valL][comput1[valL]]; 
                                    else if (valL >= (rc/2) && valL < 4319) temp1uij2[m] = uij2[valL-2468][comput1[valL]];
                                    else temp1uij1[m] = uij3[valL-4319][comput1[valL]]; 
                                }
                                else if (m >= i) {
                                    valL = M2[j-5553][m + 1] - 1;
                                    if (valL < (rc/2)) temp1uij2[m] = uij1[valL][comput1[valL]]; 
                                    else if (valL >= (rc/2) && valL < 4319) temp1uij2[m] = uij2[valL-2468][comput1[valL]];
                                    else temp1uij1[m] = uij3[valL-4319][comput1[valL]];
                                }
                            }
                            temp1uij2[5] = Lj[j];
                            qij2[i][j-5553] = temp1uij2[0] + temp1uij2[1] + temp1uij2[2] +temp1uij2[3]+temp1uij2[4]+temp1uij2[5];
                        }
                    }
                    if (j < 5553) {
                        for (m = 0; m < 3; m++) {
                            comput1[M1[j][m] - 1] += 1;
                        }
                    } else {
                        for (m =0; m < 6; m++) {
                            comput1[M2[j-5553][m] - 1] += 1;
                        }
                    }
                }

                // decision
                for (i = 0; i < comput2len; i++) comput2[i] = 0;
                for (j = 0; j < n; j++) {
                    qj[j] = Lj[j];
                    if (j < 5553) {
                        for (i = 0; i < 3; i++) {
                            valL = M1[j][i] - 1;
                            if (valL < (rc/2)) qj[j] += uij1[valL][comput2[valL]];
                            else if (valL >= (rc/2) && valL < 4319) qj[j] += uij2[valL-2468][comput2[valL]];
                            else qj[j] += uij3[valL-4319][comput2[valL]];
                        }
                        if (qj[j] >= 0) output[j] = 0;
                        else if (qj[j] < 0) output[j] = 1;
                        for (i = 0; i < 3; i++) {
                            comput2[M1[j][i] - 1] += 1;
                        }
                    } else {
                        for (i = 0; i < 6; i++) {
                            valL = M2[j-5553][i] - 1;
                            if (valL < (rc/2)) qj[j] += uij1[valL][comput2[valL]];
                            else if (valL >= (rc/2) && valL < 4319) qj[j] += uij2[valL-2468][comput2[valL]];
                            else qj[j] += uij3[valL-4319][comput2[valL]];
                            //qj[j] += uij2[valL][comput2[valL]];
                        }
                        if (qj[j] >= 0) output[j] = 0;
                        else if (qj[j] < 0) output[j] = 1;
                        for (i = 0; i < 6; i++) {
                            comput2[M2[j-4936][i] - 1] += 1;
                        }
                    }
                }

                // to check Hx=0     
                for (i = 0; i < rc; i++) {
                    checkbit[i] = 0;
                    if (i < (rc/2)) {
                        for (j = 0; j < 12; j++) {
                            checkbit[i] += output[L1[i][j] - 1];
                        }
                        checkbit[i] = checkbit[i] % 2;
                    } else if (i >= (rc/2) && i < 4319) {
                        for (j = 0; j < 5; j++) {
                            checkbit[i] += output[L2[i-2468][j] - 1];
                        }
                        checkbit[i] = checkbit[i] % 2;
                    } else {
                        for (j = 0; j < 6; j++) {
                            checkbit[i] += output[L3[i-4319][j] - 1];
                        }
                        checkbit[i] = checkbit[i] % 2;
                    }
                }

                for (i = 0; i < rc; i++) {
                    if (checkbit[i] == 0) restart += 1; // restart = 408 is success
                }
                stp = 0;
                if (k1 == 99 && restart != rc) {
                    stp = 1;
                    s++;
                }
            }
            printf("yes\n");
            /*if(k == 100) */printf("s = %d; k[%d] = %d\n", s, num, k1);
            error1 = 0;
            error2 = 0;
            for(i = 0; i < 5553; i++) {
                if (output[i] != codarray[i]) {
                    error1 += 1;
                }
            }
            for(i = 5553; i < n; i++) {
                if (output[i] != codarray[i]) {
                    error2 += 1;
                }
            }
            if ((error1+error2) != 0 && stp == 0) s++;
            restart = 0;
            if(error1 != 0) printf("error1 = %d\n", error1);
            if(error2 != 0) printf("error2 = %d\n", error2);       
            totalerror1 += error1;
            totalerror2 += error2;
        }
        double ber1;
        double ber2;
        ber1 = (double)totalerror1 / (num * 5553);
        ber2 = (double)totalerror2 / (num * 4319);
        printf("totalerror1 = %d\n", totalerror1);
        printf("totalerror2 = %d\n", totalerror2);
        printf("BER1 = %g\n", ber1);
        printf("BER2 = %g\n", ber2);
        bers1[step] = ber1;
        bers2[step] = ber2;
    }
    free(codarray);
    free(codearray);
    free(outp);
    free(output);
    for (i = 0; i < L1lenrow; i++) free(L1[i]);
    free(L1);
    for (i = 0; i < M1lenrow; i++) free(M1[i]);
    free(M1);
    for (i = 0; i < L2lenrow; i++) free(L2[i]);
    free(L2);
    for (i = 0; i < L3lenrow; i++) free(L3[i]);
    free(L3);
    for (i = 0; i < M2lenrow; i++) free(M2[i]);
    free(M2);
    for (i = 0; i < qij1row; i++) free(qij1[i]);
    free(qij1);
    for (i = 0; i < uij1row; i++) free(uij1[i]);
    free(uij1);
    for (i = 0; i < qij2row; i++) free(qij2[i]);
    free(qij2);
    
    for (i = 0; i < uij2row; i++) free(uij2[i]);
    free(uij2);
    for (i = 0; i < uij3row; i++) free(uij3[i]);
    free(uij3);
    free(Lj);
    free(comput);
    free(comput1);
    free(comput2);
    for (i = 0; i < krc; i++) free(G[i]);
    free(G);
    free(u);


    return 0;
}

void normal(double sig, double *n1, double *n2)
{   
    double x1,x2;
    double s;
    //printf("sigma = %g\n", sig);
    do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);
    *n1 = sig * x1 * sqrt((-2.0 * log(s))/ s);
    *n2 = sig * x2 * sqrt((-2.0 * log(s))/ s);
    
}

double Ranq1() {
    if ( RANI == 0 ){
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;

    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

double CHK(double L1, double L2) 
{   
    double sgn1, sgn2, min;

    if(L1>0) sgn1 = 1;
    else if(L1 == 0) sgn1 = 0;
    else sgn1 = -1;
    if(L2>0) sgn2 = 1;
    else if(L2 == 0) sgn2 = 0;
    else sgn2 = -1;
    if(fabs(L1) >= fabs(L2)) min = fabs(L2);
    else min = fabs(L1);
    double temp1,temp2;
    double ope1,ope2;
    temp1 = L1+L2;
    temp2 = L1-L2;
    ope1 = fabs(temp1);
    ope2 = fabs(temp2);
    double lnup, lndown;
    if (ope1 >= 0 && ope1 < 0.196) lnup = table[0];
    else if (ope1 >= 0.196 && ope1 < 0.433) lnup = table[1];
    else if (ope1 >= 0.433 && ope1 < 0.71) lnup = table[2];
    else if (ope1 >= 0.71 && ope1 < 1.05) lnup = table[3];
    else if (ope1 >= 1.05 && ope1 < 1.508) lnup = table[4];
    else if (ope1 >= 1.508 && ope1 < 2.252) lnup = table[5];
    else if (ope1 >= 2.252 && ope1 < 4.5) lnup = table[6];
    else if (ope1 >= 4.5) lnup = table[7];

    if (ope2 >= 0 && ope2 < 0.196) lndown = table[0];
    else if (ope2 >= 0.196 && ope2 < 0.433) lndown = table[1];
    else if (ope2 >= 0.433 && ope2 < 0.71) lndown = table[2];
    else if (ope2 >= 0.71 && ope2 < 1.05) lndown = table[3];
    else if (ope2 >= 1.05 && ope2 < 1.508) lndown = table[4];
    else if (ope2 >= 1.508 && ope2 < 2.252) lndown = table[5];
    else if (ope2 >= 2.252 && ope2 < 4.5) lndown = table[6];
    else if (ope2 >= 4.5) lndown = table[7];

    return sgn1 * sgn2 * min +lnup-lndown /*answer*/;
}