#include <stdio.h>
#include <math.h>

#define tolerance 1e-10

void UserInput();
void SystemOfODEs(double L, double P, double E, double I);
void Differentation();

void Integration(double L, double P, double E, double I);
double Delta(double L, double P, double E, double I);

void RootFinding(double L, double P, double E, double I);
double y(double x, double P, double L, double d, double E, double I);
double y_prime(double x, double P, double L, double d, double E, double I);
double y_double_prime(double x, double P, double L, double d, double E, double I);
double newton_raphson(double P, double L, double d, double E, double I, double x0);
double regula_falsi(double P, double L, double d, double E, double I, double a, double b);

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

//= ============================================================================
//= User Input
//= ============================================================================

void UserInput() {
    printf("USER INPUT:\n");
    double L, P, E_GPa, I_mm4;
    do {
        printf("Enter the column height, L (in m): ");
        scanf("%lf", &L);
        printf("Enter the concentrated load, P (in N): ");
        scanf("%lf", &P);
        printf("Enter the modulus of elasticity, E (in GPa): ");
        scanf("%lf", &E_GPa);
        printf("Enter the moment of inertia, I (in mm^4): ");
        scanf("%lf", &I_mm4);

        if (L > 0 && P > 0 && E_GPa > 0 && I_mm4 > 0) {
            break;
        }
        printf("\nAll values must be positive. Enter again.\n");
    } while(1);

    double E = E_GPa * 1e9;
    double I = I_mm4 * 1e-12;

    SystemOfODEs(L, P, E, I);
    Integration(L, P, E, I);
    RootFinding(L, P, E, I);
}

//= ============================================================================
//= System of ODEs
//= ============================================================================

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
        printf("Enter the step size, h (m): ");
        scanf("%lf", &h);

        if (0 < h && h < L) {
            break;
        }
        printf("\nh must be between 0 and L. Enter again.\n");
    } while(1);

    FILE *fPtr = fopen("data.csv", "w");
    fprintf(fPtr, "n,x (m),y (m),θ (rad),M (Nm),V (N)\n");

    double x = 0;
    double y = 0;
    double theta = 0;
    double M = P*L;
    double V = -P;

    int N = round((double) L/h);
    int n = 0;
    if (response == 1) {
        // Euler Method:
        for (n = 0; n < N; n++) {
            x = n * h;
            fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);

            y = y + h*theta;
            theta = theta + h*M/(E*I);
            M = M + h*V;
        }
        fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);
    } else {
        // RK4:
        for (n = 0; n < N; n++) {
            x = n * h;
            fprintf(fPtr, "%d,%.4e,%.4e,%.4e,%.4e,%.4e\n", n, x, y, theta, M, V);

            double f1 = h * theta;
            double g1 = h * M/(E*I);
            double h1 = h * V;

            double f2 = h * (theta + 0.5*g1);
            double g2 = h * (M + 0.5*h1)/(E*I);
            double h2 = h * V;

            double f3 = h * (theta + 0.5*g2);
            double g3 = h * (M + 0.5*h2)/(E*I);
            double h3 = h * V;

            double f4 = h * (theta + g3);
            double g4 = h * (M + h3)/(E*I);
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

//= ============================================================================
//= Differentiation
//= ============================================================================

void Differentation() {

}

//= ============================================================================
//= Integration
//= ============================================================================

void Integration(double L, double P, double E, double I) {
    printf("\nNUMERICAL INTEGRATION:\n");
    printf("Choose:\n");
    printf("(1)  Composite Trapezoidal Rule\n");
    printf("(2)  Composite Simpson's Rule\n\n");

    int response = 0;
    do {
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

        if (n > 0 && n % 2 == 0) {
            break;
        }
        printf("\nn must be positive and even. Enter again.\n");
    } while (1);

    double maxP = P;
    double h = maxP/n;
    double area = 0;

    if (response == 1) {
        // Composite Trapezoidal Rule:

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
        // Composite Simpson's Rule:

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
    printf("Exact Area: %lf\n", exactArea);
}

//== ===========================================================================
//== Integration Helper Functions
//== ===========================================================================

double Delta(double L, double P, double E, double I) {
    return P*L*L*L / (3*E*I);
}

//= ============================================================================
//= Root Finding
//= ============================================================================

void RootFinding(double L, double P, double E, double I) {
    printf("\nROOT FINDING:\n");

    double d = 0;
    do {
        printf("Enter the pin-support distance, d (m): ");
        scanf("%lf", &d);

        if (0 < d && d < L) {
            break;
        }
        printf("\nd must be between 0 and L. Enter again.\n");
    } while(1);

    int response = 0;
    printf("Choose:\n");
    printf("(1)  Newton-Raphson\n");
    printf("(2)  Regula Falsi\n\n");
    do {
        printf("Respond 1 or 2: ");
        scanf("%d", &response);

        if (response == 1 || response == 2) {
            break;
        }
        printf("\nWrong Input. Enter again.\n");
    } while(1);

    double delta_at_free_end = fabs((P * (L - d)) / (E * I) * (-pow(L, 2) / 4.0 + pow(L, 3) / (4.0 * d) - (pow(L - d, 2) * (3 * L - d) / (12.0 * d))));
    // printf("\nMaximum deflection at the free end (in meters): %.6f\n", delta_free_end);
    double delta_max_within_d = 0;

    if (response == 1) {
        // printf("\nRunning Newton-Raphson for three initial guesses:\n");

        double delta = 0;
        double delta_new = 0;

        double guesses[3] = {0, d/2, d};
        for (int i = 0; i < 3; i++) {
            double x = newton_raphson(P, L, d, E, I, guesses[i]);
            delta_new = fabs(y(x, P, L, d, E, I));
            if (delta_new > delta) {
                delta_max_within_d = delta_new;
            }
        }
        // printf("Deflection at x = %.6f is %.6f m\n", x, delta);

        // x = newton_raphson(P, L, d, E, I, d / 2);
        // printf("Maximum deflection at x = %.6f: %.6f m\n", x, y(x, P, L, d, E, I));

        // x = newton_raphson(P, L, d, E, I, d);
        // printf("Maximum deflection at x = %.6f: %.6f m\n", x, y(x, P, L, d, E, I));
    } else {
        double a = d / 100.0;
        double b = d;

        double x = regula_falsi(P, L, d, E, I, a, b);
        delta_max_within_d = fabs(y(x, P, L, d, E, I));

        // if (x != -1) {
        //     printf("\nLocation of maximum deflection within [0, d] (in meters): %.6f\n", x);
        //     printf("Maximum deflection: %.6f\n", y(x, P, L, d, E, I));
        // } else {
        //     printf("\nCalculation failed.\n");
        // }
    }

    printf("\nMax deflection within d: %.6lf m\n", delta_max_within_d);
    printf("Max deflection at free end: %.6lf m\n", delta_at_free_end);

    if (delta_max_within_d > delta_at_free_end) {
        printf("Max deflection occurs within d\n");
    } else {
        printf("Max deflection occurs at free end\n");
    }
}

//== ===========================================================================
//== Root Finding Helper Functions
//== ===========================================================================

double y(double x, double P, double L, double d, double E, double I) {
    return (P * (L - d)) / (4 * E * I) * (-pow(x, 2) + (1.0 / d) * pow(x, 3));
}

double y_prime(double x, double P, double L, double d, double E, double I) {
    return (P * (L - d)) / (4 * E * I) * (-2 * x + (3.0 / d) * pow(x, 2));
}

double y_double_prime(double x, double P, double L, double d, double E, double I) {
    return (P * (L - d)) / (4 * E * I) * (-2 + (6.0 / d) * x);
}

double newton_raphson(double P, double L, double d, double E, double I, double x0) {
    double x = x0;
    double fx, fprime_x;
    int iter = 0;
    double x_old = 0;

    // printf("\nNewton-Raphson Method with initial guess x0 = %.6f:\n", x0);
    // printf("Iteration\t x\t\t y'(x)\t\t y''(x)\n");
    // printf("----------------------------------------------------------------\n");

    do {
        fx = y_prime(x, P, L, d, E, I);
        fprime_x = y_double_prime(x, P, L, d, E, I);

        // printf("%d\t\t %.6f\t %.6e\t %.6e\n", iter, x, fx, fprime_x);

        if (fabs(fprime_x) < tolerance) {
            // printf("Derivative near zero, cannot proceed with the current method. Please try again!\n");
            return -1;
        }
        x_old = x;
        x -= fx / fprime_x;
        iter++;
    } while (fabs(x - x_old) > tolerance && iter < 100);

    // printf("----------------------------------------------------------------\n");
    // printf("Converged to x = %.6f after %d iterations.\n", x, iter);
    return x;
}

double regula_falsi(double P, double L, double d, double E, double I, double a, double b) {
    double fa = y_prime(a, P, L, d, E, I);
    double fb = y_prime(b, P, L, d, E, I);

    if (fa * fb > 0) {
        // printf("Function does not change sign in the interval [%.6f, %.6f]. Please try again!\n", a, b);
        // printf("y'(a) = %.6e, y'(b) = %.6e\n", fa, fb);
        return -1;
    }

    double c;
    int iter = 0;

    // printf("\nRegula Falsi Method in interval [%.6f, %.6f]:\n", a, b);
    // printf("\nIteration\t a\t\t b\t\t c\t\t y'(c)\n");
    // printf("-------------------------------------------------------------------------------\n");

    do {
        c = (a * fb - b * fa)/ (fb - fa);
        double fc = y_prime(c, P, L, d, E, I);

        // printf("%d\t\t %.6f\t %.6f\t %.6f\t %.6e\n", iter, a, b, c, fc);

        if (fabs(fc) < tolerance) break;
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
        iter++;
    } while (fabs(b - a) > tolerance && iter < 100);

    // printf("-------------------------------------------------------------------------------\n");
    // printf("Converged to x = %.6f after %d iterations.\n", c, iter);
    return c;
}
