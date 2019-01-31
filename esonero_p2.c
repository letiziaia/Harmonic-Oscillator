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
    double alpha= 10.0, x0= -0.1, v0= 0.0;
    double tempo, dt=0.01, dx=0.001;
    int passi, i, j=0, solved=0;
    FILE *fpt1, *fpt2, *fpv1, *fpv2, *cm;
	double x0_lim_rk = 0, x0_lim_va = 0;
	double x_lim = 0.5;

	/* Inserimenti parametri da tastiera */
    printf("Studio di un moto unidimensionale\n\n");
    printf("inserire il tempo totale di integrazione: ");
    scanf("%lf", &tempo);
    printf("\n");

    printf("inserire la x limite (opzione base calcolata: 0.5) : ");
    scanf("%lf", &x_lim);
    printf("\n");

	/* Calcolo del numero di step */
    passi=(long int)(tempo/dt);
	
	/* Settaggio valori iniziali */
	init.x=x0;
	init.v=v0;
	
	/* Decremento x0 e controllo finchè non supera il massimo locale */
	while ((solved<2) && (j<1000)) {
		/* Settaggio valori iniziali */
		j++;
		init.x=x0-j*dx;
		init.v=v0;
		printf("x0: %lf\n", init.x);
    
		xEv1.x=init.x;
		xEv1.v=init.v;
		xEv2.x=init.x;
		xEv2.v=init.v;
	
		/* Ciclo di esecuzione degli algoritmi */
		for(i=0; i<passi; i++){
			/* Esecuzione Runge-Kutta */
			xEv1=rungekutta(alpha, dt, xEv1);
			if (xEv1.x>x_lim) {
				solved++;
				x0_lim_rk = init.x;
			}
			/* Esecuzione Verlet */
			xEv2=verletauto(alpha, dt, xEv2);
			if (xEv2.x>x_lim) {
				solved++;
				x0_lim_va = init.x;
			}
		}
	}


    printf("Eseguite le integrazioni.\n");
	
	printf("Valore limite calcolato con Runge-Kutta: %lf.\n", x0_lim_rk);
	printf("Valore limite calcolato con Verlet     : %lf.\n", x0_lim_rk);
            
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
       star_xEv.v=f2(alpha,c)*dt; 
       
        
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
