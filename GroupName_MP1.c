#include <stdio.h>
#include <math.h>

// void can be replaced if there's a return value.
// Add parameters to the functions as needed:
void UserInput();
void SystemOfODEs(double L, double P, double E, double I);
void Differentation(double L, double P, double E, double I);
void Integration(double L, double P, double E, double I);
void RootFinding();

double Delta(double L, double P, double E, double I);

int main() {
    // Uncomment for quick testing (don't have to manually input values):
    // double L = 5;
    // double P = 1000;
    // double E = 200;
    // double I = 50000000;
    // Integration(L, P, E, I);

    UserInput();
    return 0;
}

void UserInput() {
    printf("USER INPUT:\n");
    double L, P, E, I;
    do {
        printf("Enter the column height, L (in m): ");
        scanf("%lf", &L);
        printf("Enter the concentrated load, P (in N): ");
        scanf("%lf", &P);
        printf("Enter the modulus of elasticity, E (in GPa): ");
        scanf("%lf", &E);
        printf("Enter the moment of inertia, I (in mm^4): ");
        scanf("%lf", &I);

        if (L > 0 && P > 0 && E > 0 && I > 0) {
            break;
        }
        printf("\nAll values must be positive. Enter again.\n");
    } while(1);

    SystemOfODEs(L, P, E, I);
    Integration(L, P, E, I);
}

void SystemOfODEs(double L, double P, double E, double I) {
    printf("\nSYSTEM OF ODES:\n");
    printf("Choose:\n");
    printf("(1) Euler Method\n");
    printf("(2) RK4\n\n");

    int response = 0;
    do {
        printf("Respond 1 or 2: ");
        scanf("%d", &response);

        if (response == 1 || response == 2) {
            break;
        }
        printf("\nWrong Input. Enter again.\n");
    } while(1);

    double h = 0;
    do {
        printf("Enter the step size, h: ");
        scanf("%lf", &h);

        if (h > 0) {
            break;
        }
        printf("\nh must be positive. Enter again.\n");
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
    if (response == 1) {  // Euler Method
        for (n = 0; n < N; n++) {
            x = n * h;
            fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);

            y = y + h*theta;
            theta = theta + h*M/EI;
            M = M + h*V;
        }
        fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);
    } else {  // RK4
        for (n = 0; n < N; n++) {
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
    printf("\nOutput written to data.csv.\n");
}

void Differentation(double L, double P, double E, double I) {

}

void Integration(double L, double P, double E, double I) {
    printf("\nNUMERICAL INTEGRATION:\n");
    int response = 0;
    do {
        printf("Choose:\n");
        printf("(1)  Trapezoidal Rule\n");
        printf("(2)  Simpson's Rule\n\n");
        printf("Respond 1 or 2: ");
        scanf("%d", &response);

        if (response == 1 || response == 2) {
            break;
        }
        printf("\nWrong Input. Enter again.\n");
    } while(1);

    int n = 0;
    do {
        printf("Enter the number of intervals, n: ");
        scanf("%d", &n);

        if (n > 0) {
            break;
        }
        printf("\nn must be positive. Enter again.\n");
    } while (1);

    double maxP = P;
    double h = maxP/n;
    double area = 0;

    if (response == 1) {
        double endNodesSum = Delta(L, maxP, E, I) + Delta(L, 0, E, I);
        double innerNodesSum = 0;
        for (int k = 1; k < n; k++) {
            P = h*k;
            innerNodesSum += Delta(L, P, E, I);
        }

        double endNodesTerm = h/2 * endNodesSum;
        double innerNodesTerm = h * innerNodesSum;
        area = endNodesTerm + innerNodesTerm;
    } else {
        double endNodesSum = Delta(L, maxP, E, I) + Delta(L, 0, E, I);
        double evenNodesSum = 0;
        double oddNodesSum = 0;
        for (int k = 1; k < n; k++) {
            P = h*k;
            if (k % 2 == 0) {
                evenNodesSum += Delta(L, P, E, I);
            } else {
                oddNodesSum += Delta(L, P, E, I);
            }
        }

        double endNodesFactor = (double) h / 3;
        double evenNodesFactor = 2 * (double) h / 3;
        double oddNodesFactor = 4 * (double) h / 3;

        double endNodesTerm = endNodesFactor * endNodesSum;
        double evenNodesTerm = evenNodesFactor * evenNodesSum;
        double oddNodesTerm = oddNodesFactor * oddNodesSum;

        area = endNodesTerm + evenNodesTerm + oddNodesTerm;
    }

    double exactArea = 0.5 * maxP * Delta(L, maxP, E, I);
    printf("\nCalculated Area: %lf\n", area);
    printf("Exact Area: %lf", exactArea);
}

// I'm putting RootFinding() in its own file, RootFinding.c
// Testing RootFinding.c and this current file separately will lessen the chance of bugs.
// Run test_RootFinding.c to test RootFinding.c

#include "RootFinding.c"

// Before we submit, we'll replace `#include "RootFinding.c"` with the actual contents
// of RootFinding.c because the Machine Problem needs to be submitted as one file.

double Delta(double L, double P, double E, double I) {
    double EI = E*I*1e-3;
    return P*L*L*L / (3*EI);
}
