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
    double L, P, E, I;
    do {
        printf("Input the column height L (in m): ");
        scanf("%lf", &L);
        printf("Input the concentrated load P (in N): ");
        scanf("%lf", &P);
        printf("Input the modulus of elasticity E (in GPa): ");
        scanf("%lf", &E);
        printf("Input the moment of inertia I (in mm^4): ");
        scanf("%lf", &I);

        if (L > 0 && P > 0 && E > 0 && I > 0) {
            break;
        }
        printf("\nAll values must be positive. Enter again.\n");
    } while(1);

    SystemOfODEs(L, P, E, I);
}

void SystemOfODEs(double L, double P, double E, double I) {
    double h, response;

    do {
        printf("Input the step size h: ");
        scanf("%lf", &h);

        if (h > 0) {
            break;
        }
        printf("\nAll values must be positive. Enter again.\n");
    } while(1);

    printf("\nChoose a method:\n");
    printf("(1) Euler's\n");
    printf("(2) 4th Order Runge-Kutta\n");
    do {
        printf("Respond 1 or 2: ");
        scanf("%lf", &response);

        if (response == 1 || response == 2) {
            break;
        }
        printf("\nWrong Input. Enter again.\n");
    } while(1);

    FILE *fPtr = fopen("data.csv", "w");
    fprintf(fPtr, "n,x (m),y (m),θ (rad),M (Nm),V (N)\n");

    double x = 0;
    double y = 0;
    double theta = 0;
    double M = P*L;
    double V = -P;
    double EI = E*I*1e-3;

    int N = round((double) L/h);
    int n = 0;
    if (response == 1){
        for (n = 0; n < N; n++){
            x = n * h;
            fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);

            y = y + h*theta;
            theta = theta + h*M/EI;
            M = M + h*V;
        }
        fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);
    }
    else if (response == 2){
        for (n = 0; n < N; n++){
            x = n * h;
            fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);

            double f1 = h * theta;
            double g1 = h * M/EI;
            double h1 = h * V;

            double f2 = h * (theta + 0.5*g1);
            double g2 = h * (M + 0.5*h1)/EI;
            double h2 = h * V;

            double f3 = h * (theta + 0.5*g2);
            double g3 = h * (M + 0.5*h2)/EI;
            double h3 = h * V;

            double f4 = h * (theta + g3);
            double g4 = h * (M + h3)/EI;
            double h4 = h * V;

            y = y + (f1 + 2*f2 + 2*f3 + f4)/6;
            theta = theta + (g1 + 2*g2 + 2*g3 + g4)/6;
            M = M + (h1 + 2*h2 + 2*h3 + h4)/6;
        }
        fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);
    }
    fclose(fPtr);

    printf("\ny = %lf", y);
    printf("\ntheta = %lf", theta);
    printf("\nM = %lf", M);
    printf("\nV = %lf", V);
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
