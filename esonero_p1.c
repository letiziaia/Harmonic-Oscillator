#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define A 30
#define B 20

/* Strutture per contenere i dati delle fasi */
struct sp_fasi{
       double x;
       double v;
       } xEv1, xEv2, init;

 
/* Prototipi delle funzioni */
double f2(double alpha, double x);            /*eq del moto*/
struct sp_fasi rungekutta(double alpha, double dt, struct sp_fasi old_xEv);  /*integra con rungekutta*/
struct sp_fasi verletauto(double alpha, double dt, struct sp_fasi old_xEv);  /*integra con verlet autosufficiente*/


int main(void){
    double alpha, x0, v0;
    double tempo, dt=0.01;
    int passi, i;
    FILE *fpt1, *fpt2, *fpv1, *fpv2, *cm;

	/* Inserimenti parametri da tastiera */
	
    printf("Studio di un moto unidimensionale\n\n");
    printf("inserire il tempo totale di integrazione: ");
    scanf("%lf", &tempo);
    printf("\n");

    printf("inserire il valore di x0 (opzione base -0.1): ");
    scanf("%lf", &x0);
    printf("\n");
    
    printf("inserire il valore di v0 (opzione base 0.0): ");
    scanf("%lf", &v0);
    printf("\n");

    printf("inserire il valore di alpha (opzione base 10.0): ");
    scanf("%lf", &alpha);
    printf("\n");

	/* Calcolo del numero di step */
    passi=(long int)(tempo/dt);
     
	/* Apertura dei file di risultati in scrittura */
    if((fpt1=fopen("t_rk.dat", "w"))==NULL){
        printf("Errore nell'apertura del file 't_rk.dat'\n");
        exit(EXIT_FAILURE);
    }
        
    if((fpt2=fopen("t_va.dat", "w"))==NULL){
        printf("Errore nell'apertura del file 't_va.dat'\n");
        exit(EXIT_FAILURE);
    }

    if((fpv1=fopen("v_rk.dat", "w"))==NULL){
        printf("Errore nell'apertura del file 'v_rk.dat'\n");
        exit(EXIT_FAILURE);
    }

    if((fpv2=fopen("v_va.dat", "w"))==NULL){
        printf("Errore nell'apertura del file 'v_va.dat'\n");
        exit(EXIT_FAILURE);
    }
	
	/* Settaggio valori iniziali */
    init.x=x0;
    init.v=v0;
    
    fprintf(fpt1, "%lf %lf\n", 0, init.x);
    fprintf(fpt2, "%lf %lf\n", 0, init.x);
    fprintf(fpv1, "%lf %lf\n", init.x, init.v);
	fprintf(fpv2, "%lf %lf\n", init.x, init.v);
	
    xEv1.x=init.x;
    xEv1.v=init.v;
    xEv2.x=init.x;
    xEv2.v=init.v;
	
	/* Ciclo di esecuzione degli algoritmi */
    for(i=0; i<passi; i++){
		/* Esecuzione e salvataggio risultati di Runge-Kutta */
        xEv1=rungekutta(alpha, dt, xEv1);
        fprintf(fpt1, "%lf %lf\n", i*dt, xEv1.x);
		fprintf(fpv1, "%lf %lf\n", xEv1.x, xEv1.v);
		/*Esecuzione e salvataggio risultati di Verlet */
        xEv2=verletauto(alpha, dt, xEv2);
        fprintf(fpt2, "%lf %lf\n", i*dt, xEv2.x);
		fprintf(fpv2, "%lf %lf\n", xEv2.x, xEv2.v);
	}

	/*Chiusura dei file */
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpv1);
    fclose(fpv2);

    printf("Eseguite le integrazioni. I risultati sono in salvati nei file .dat in formato compatibile con gnuplot.\n");    
    
	/* Creazione file script per GnuPlot */
     if((cm=fopen("comandignuplot.dat", "w"))==NULL){      
        printf("Errore nell'apertura del file 'comandignuplot.dat'\n");
        exit(EXIT_FAILURE);
        }
        
    fprintf(cm, "set terminal gif\n");
    fprintf(cm, "set output \"plot_trk.gif\"\n");
    fprintf(cm, "plot \"t_rk.dat\" with lines\n");
    fprintf(cm, "set output \"plot_tva.gif\"\n");
    fprintf(cm, "plot \"t_va.dat\" with lines\n");
    fprintf(cm, "set terminal gif\n");
    fprintf(cm, "set output \"plot_vrk.gif\"\n");
    fprintf(cm, "plot \"v_rk.dat\" with lines\n");
    fprintf(cm, "set output \"plot_vva.gif\"\n");
    fprintf(cm, "plot \"v_va.dat\" with lines\n");
    fclose(cm);

	/* Esecuzione di GnuPlot */
    system("gnuplot.exe comandignuplot.dat");           
        
    printf("\nFINE\n");
 
    return(0);
}


double f2(double alpha, double x){
       return(-alpha*x+A*x*x-B*x*x*x);
}

struct sp_fasi rungekutta(double alpha, double dt, struct sp_fasi old_xEv){
       struct sp_fasi newrk_xEv, star_xEv;
       double c;
       
       star_xEv.x=(old_xEv.v)*dt;
       c=old_xEv.x+0.5*(star_xEv.x);
       star_xEv.v=f2(alpha,c)*dt; /* ?? */
       
        
       newrk_xEv.x=old_xEv.x+(old_xEv.v + 0.5*star_xEv.v)*dt;
       newrk_xEv.v=old_xEv.v+f2(alpha,c)*dt;
           
       return(newrk_xEv);
       }
       
       
struct sp_fasi verletauto(double alpha, double dt, struct sp_fasi old_xEv){
       struct sp_fasi newv_xEv;
       double c, c1;
       
       c=old_xEv.x;
       newv_xEv.x=old_xEv.x+old_xEv.v*dt+0.5*(f2(alpha,c))*dt*dt;
	   
       c1=newv_xEv.x;
       newv_xEv.v=old_xEv.v+dt*0.5*((f2(alpha,c))+(f2(alpha,c1)));
       
       return(newv_xEv);
       }
