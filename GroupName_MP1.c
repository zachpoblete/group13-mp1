#include <stdio.h>
#include <math.h>

// void can be replaced if there's a return value.
// Add parameters to the functions as needed:
void UserInput();
void SystemOfODEs();
void Differentation();
void Integration();
void RootFinding();

int main() {
    UserInput();
}

void UserInput() {
    // ! TO DO: Input validation
    double L, P, E, I;
    printf("Input the column height L (in m): ");
    scanf("%lf", &L);
    printf("Input the concentrated load P (in N): ");
    scanf("%lf", &P);
    printf("Input the modulus of elasticity E (in GPa): ");
    scanf("%lf", &E);
    printf("Input the moment of inertia I (in mm⁴): ");
    scanf("%lf", &I);
    SystemOfODEs(L, P, E, I);
}

void SystemOfODEs(double L, double P, double E, double I) {
    double x, choice;
    double h, EI;
    double y, theta, M, V;
    double f1, f2, f3, f4;
    double g1, g2, g3, g4;
    double h1, h2, h3 ,h4;
    double i1, i2, i3, i4;

    printf("Input the step size h: ");
    scanf("%lf", &h);

    do{
        printf("\nChoose a method:\n1. Euler's \n2. 4th Order Runge-Kutta \n\nEnter your choice: ");
        scanf("%lf", &choice);
    } while(choice != 1 && choice !=2);

    FILE *fPtr = fopen("data.csv", "w");
    fprintf(fPtr, "n,x (m),y (m),θ (rad),M (Nm),V (N)\n");

    x = 0, y = 0, theta = 0, M = P*L, V = -P, EI = E*I*1e-3;

    int N = round((double) L/h);
    int n = 0;
    if (choice == 1){
        for (n = 0; n < N; n++){
            x = n * h;
            fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);

            y = y + h*theta;
            theta = theta + h*M/EI;
            M = M + h*V;
        }
        fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);
    }
    else if (choice == 2){
        for (n = 0; n < N; n++){
            x = n * h;
            fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);

            f1 = h*theta;
            g1 = h*M/EI;
            h1 = h*V;

            f2 = h*(theta+0.5*g1);
            g2 = h*(M+0.5*h1)/EI;
            h2 = h*V;

            f3 = h*(theta+0.5*g2);
            g3 = h*(M+0.5*h2)/EI;
            h3 = h*V;

            f4 = h*(theta+g3);
            g4 = h*(M+h3)/EI;
            h4 = h*V;

            y = y + (f1+2*f2+2*f3+f4)/6;
            theta = theta + (g1+2*g2+2*g3+g4)/6;
            M = M + (h1+2*h2+2*h3+h4)/6;
        }

        fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);
    }

    printf("\ny = %lf", y);
    printf("\ntheta = %lf", theta);
    printf("\nM = %lf", M);
    printf("\nV = %lf", V);

    fclose(fPtr);
    printf("\nThe results are successfully saved to data.csv!\n");
}

void Differentation() {

}

void Integration() {

}

// I'm putting RootFinding() in its own file, RootFinding.c
// Testing RootFinding.c and this current file separately will lessen the chance of bugs.
// Run test_RootFinding.c to test RootFinding.c

#include "RootFinding.c"

// Before we submit, we'll replace `#include "RootFinding.c"` with the actual contents
// of RootFinding.c because the Machine Problem needs to be submitted as one file.
