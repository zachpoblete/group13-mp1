#include <stdio.h>

// void can be replaced if there's a return value.
// Add parameters to the functions as needed:
void UserInput();
void SystemOfODEs();
void Differentation();
void Integration();
void RootFinding();

int main() {
    printf("Hello groupmates!");
    return 0;
}

void UserInput() {
    
        double i, choice;
        double h, L, P, E, I;
        double z1, z2, z3, z4;
        double f1, f2, f3, f4;
        double g1, g2, g3, g4;
        double h1, h2, h3 ,h4;
        double i1, i2, i3, i4;
    
        printf("Enter P: ");
        scanf("%lf", &P);
        printf("Enter L: ");
        scanf("%lf", &L);
        printf("Enter E: ");
        scanf("%lf", &E);
        printf("Enter I: ");
        scanf("%lf", &I);
        printf("Enter h: ");
        scanf("%lf", &h);
    
        do{
            system("cls");
            printf("Choose a method:\n1. Euler's \n2. 4th Order Runge-Kutta \nEnter your choice: ");
            scanf("%lf", &choice);
        } while(choice != 1 && choice !=2);
    
        z1 = 0, z2 = 0, z3 = (P/1000.0)*L, z4 = -(P/1000.0);
    
        if(choice == 1){
          for (i=0; i <= (L-h); i+=h){
              z1 = z1 + h*z2;
              z2 = z2 + h*z3;
              z3 = z3 + h*z4;
            }
    
          if (i == L){
                z4 = z4 + h*(P/(E*(I/pow(10,6))));
                }
        }
    
        else if(choice == 2){
            for (i; i <= (L-h); i+=h){
             f1 = h*z2;
             g1 = h*z3;
             h1 = h*z4;
    
             f2 = h*(z2+0.5*f1);
             g2 = h*(z3+0.5*g1);
             h2 = h*(z4+0.5*h1);
    
             f3 = h*(z2+0.5*f2);
             g3 = h*(z3+0.5*g2);
             h3 = h*(z4+0.5*h2);
    
             f4 = h*(z2+f3);
             g4 = h*(z3+g3);
             h4 = h*(z4+h3);
    
             z1 = z1 + (1.0/6.0)*(f1+2*f2+2*f3+f4);
             z2 = z2 + (1.0/6.0)*(g1+2*g2+2*g3+g4);
             z3 = z3 + (1.0/6.0)*(h1+2*h2+2*h3+h4);
    
        }
        if (i == L){
                f1 = h*z2;
                g1 = h*z3;
                h1 = h*z4;
                i1 = h*(P/(E*(I/pow(10,6))));
    
                f2 = h*(z2+0.5*f1);
                g2 = h*(z3+0.5*g1);
                h2 = h*(z4+0.5*h1);
                i2 = h*((P+0.5*i1)/(E*(I/pow(10,6))));
    
                f3 = h*(z2+0.5*f2);
                g3 = h*(z3+0.5*g2);
                h3 = h*(z4+0.5*h2);
                i3 = h*((P+0.5*i2)/(E*(I/pow(10,6))));
    
                f4 = h*(z2+f3);
                g4 = h*(z3+g3);
                h4 = h*(z4+h3);
                i4 = h*((P+i3)/(E*(I/pow(10,6))));
    
                z1 = z1 + (1.0/6.0)*(f1+2*f2+2*f3+f4);
                z2 = z2 + (1.0/6.0)*(g1+2*g2+2*g3+g4);
                z3 = z3 + (1.0/6.0)*(h1+2*h2+2*h3+h4);
                z4 = z4 + (1.0/6.0)*(i1+2*i2+2*i3+i4);
                }
    
        }
    
    
        printf("\ny = %lf", z1/(E*(I/pow(10,6))));
        printf("\ntheta = %lf", z2*(M_PI/180));
        printf("\nM = %lf", z3);
        printf("\nV = %lf", z4*1000);
    }


// Assignee:
void SystemOfODEs() {

}

// Assignee:
void Differentation() {

}

// Assignee:
void Integration() {

}

// I'm putting RootFinding() in its own file, RootFinding.c
// Testing RootFinding.c and this current file separately will lessen the chance of bugs.
// Run test_RootFinding.c to test RootFinding.c

// Assignee: Raoul
#include "RootFinding.c"

// Before we submit, we'll replace `#include "RootFinding.c"` with the actual contents
// of RootFinding.c because the Machine Problem needs to be submitted as one file.
