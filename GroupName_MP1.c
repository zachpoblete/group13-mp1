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
 
    }


// Assignee:
void SystemOfODEs() {
    double i, choice;
    double h, L, P, E, I, E_I;
    double z1, z2, z3, z4;
    double f1, f2, f3, f4;
    double g1, g2, g3, g4;
    double h1, h2, h3 ,h4;
    double i1, i2, i3, i4;

    FILE *fp = fopen("rk4.csv", "w");
    if (!fp) {
        perror("Error opening file");
        return;
    }

    fprintf(fp, "x (m), Deflection y (m), Slope θ (rad), Moment M (Nm), Shear V (N)\n");

    printf("Input the height of the column (in m): ");
    scanf("%lf", &P);
    printf("Input the concentrated load P (in N): ");
    scanf("%lf", &L);
    printf("Input the modulus of elasticity (in GPa): ");
    scanf("%lf", &E);
    printf("Input the moment of inertia (in mm4): ");
    scanf("%lf", &I);
    printf("Input the step size: ");
    scanf("%lf", &h);

    do{
        system("cls");
        printf("Choose a method:\n1. Euler's \n2. 4th Order Runge-Kutta \nEnter your choice: ");
        scanf("%lf", &choice);
    } while(choice != 1 && choice !=2);

    z1 = 0, z2 = 0, z3 = (P/1000.0)*L, z4 = -(P/1000.0), E_I = E*I*1e-6;

    if(choice == 1){
      for (i=0; i <= L; i+=h){
          z1 = z1 + h*z2;
          z2 = z2 + h*z3/E_I;
          z3 = z3 + h*z4;
        }
    }

    else if(choice == 2){
        for (i=0; i <= L; i+=h){
         fprintf(fp, "%.3e, %.3e, %.3e, %.3e, %.3e\n", i, z1, z2, z3, z4);

         f1 = h*z2;
         g1 = h*z3/(E_I);
         h1 = h*z4;

         f2 = h*(z2+0.5*g1);
         g2 = h*(z3+0.5*h1)/E_I;
         h2 = h*(z4);

         f3 = h*(z2+0.5*g2);
         g3 = h*(z3+0.5*h2)/E_I;
         h3 = h*(z4);

         f4 = h*(z2+g3);
         g4 = h*(z3+h3)/E_I;
         h4 = h*(z4);

         z1 = z1 + (f1+2*f2+2*f3+f4)/6;
         z2 = z2 + (g1+2*g2+2*g3+g4)/6;
         z3 = z3 + (h1+2*h2+2*h3+h4)/6;

    }

    }


    printf("\ny = %lf", z1);
    printf("\ntheta = %lf", z2);
    printf("\nM = %lf", z3+h);
    printf("\nV = %lf", z4*1000);
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
