#ifndef _NR_H_
    #define _NR_H_
    #ifndef _FCOMPLEX_DECLARE_T_
        typedef struct FCOMPLEX {float r,i;} fcomplex;
        #define _FCOMPLEX_DECLARE_T_
    #endif /* _FCOMPLEX_DECLARE_T_ */
    #ifndef _ARITHCODE_DECLARE_T_
        typedef struct {
            unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad; }arithcode;
        #define _ARITHCODE_DECLARE_T_
    #endif /* _ARITHCODE_DECLARE_T_ */
    #ifndef _HUFFCODE_DECLARE_T_
        typedef struct {
            unsigned long *icod,*ncod,*left,*right,nch,nodemax; } huffcode;
        #define _HUFFCODE_DECLARE_T_
    #endif /* _HUFFCODE_DECLARE_T_ */
    #include <stdio.h>
    #if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */
        void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
        float qgaus(float (*func)(float), float a, float b);

    #else /* ANSI */
/* traditional - K&R */
        void addint();
        void airy();
//Rest of traditional declarations are here on the diskette.
    #endif /* ANSI */
#endif /* _NR_H_ */