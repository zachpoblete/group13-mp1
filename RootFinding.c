#include <math.h>
#include <stdio.h>
#define tolerance 1e-10

void calculate_deflection() {
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

        printf("\nNewton-Raphson Method with initial guess x0 = %.6f:\n", x0);
        printf("Iteration\t x\t\t y'(x)\t\t y''(x)\n");
        printf("--------------------------------------------------------------\n");

        do {
            fx = y_prime(x, P, L, d, E, I);
            fprime_x = y_double_prime(x, P, L, d, E, I);

            printf("%d\t\t %.6f\t %.6e\t %.6e\n", iter, x, fx, fprime_x);

            if (fabs(fprime_x) < tolerance) {
                printf("Derivative near zero, cannot proceed.\n");
                return -1;
            }
            x -= fx / fprime_x;
            iter++;
        } while (fabs(fx) > tolerance && iter < 100);

        printf("--------------------------------------------------------------\n");
        printf("Converged to x = %.6f after %d iterations.\n", x, iter);
        return x;
    }

    double regula_falsi(double P, double L, double d, double E, double I, double a, double b) {
        double fa = y_prime(a, P, L, d, E, I);
        double fb = y_prime(b, P, L, d, E, I);

        if (fa * fb > 0) {
            printf("Function does not change sign in the interval [%.6f, %.6f].\n", a, b);
            printf("y'(a) = %.6e, y'(b) = %.6e\n", fa, fb);
            return -1;
        }

        double c;
        int iter = 0;

        printf("\nRegula Falsi Method in interval [%.6f, %.6f]:\n", a, b);
        printf("\nIteration\t a\t\t b\t\t c\t\t y'(c)\n");
        printf("--------------------------------------------------------------\n");

        do {
            c = a - (fa * (b - a)) / (fb - fa);
            double fc = y_prime(c, P, L, d, E, I);

            printf("%d\t\t %.6f\t %.6f\t %.6f\t %.6e\n", iter, a, b, c, fc);

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

        printf("--------------------------------------------------------------\n");
        printf("Converged to x = %.6f after %d iterations.\n", c, iter);
        return c;
    }

    double P, L, d, E_GPa, I_mm4;

    // Input parameters
    printf("Enter the load P (in Newtons): ");
    scanf("%lf", &P);
    printf("Enter the length L (in meters): ");
    scanf("%lf", &L);
    printf("Enter the distance d (in meters): ");
    scanf("%lf", &d);
    printf("Enter the modulus of elasticity E (in GPa): ");
    scanf("%lf", &E_GPa);
    printf("Enter the moment of inertia I (in mm^4): ");
    scanf("%lf", &I_mm4);

    double E = E_GPa * 1e9;
    double I = I_mm4 * 1e-12;

    double delta_free_end = (P * (L - d)) / (E * I) * (-pow(L, 2) / 4.0 + pow(L, 3) / (4.0 * d) - (pow(L - d, 2) * (3 * L - d) / (12.0 * d)));
    printf("\nMaximum deflection at the free end (in meters): %.6f\n", delta_free_end);

    int choice;
    printf("\nChoose a method:\n1. Newton-Raphson (three guesses)\n2. Regula Falsi\nEnter your choice: ");
    scanf("%d", &choice);

    double x_max;
    if (choice == 1) {
        printf("\nRunning Newton-Raphson for three initial guesses:\n");

        x_max = newton_raphson(P, L, d, E, I, 0.0);
        printf("Maximum deflection at x = %.6f: %.6f m   (did not converge)\n", x_max, y(x_max, P, L, d, E, I));

        x_max = newton_raphson(P, L, d, E, I, d / 2);
        printf("Maximum deflection at x = %.6f: %.6f m\n", x_max, y(x_max, P, L, d, E, I));

        x_max = newton_raphson(P, L, d, E, I, d);
        printf("Maximum deflection at x = %.6f: %.6f m\n", x_max, y(x_max, P, L, d, E, I));
    } else if (choice == 2) {
        double a = d / 100.0;
        double b = d;

        x_max = regula_falsi(P, L, d, E, I, a, b);

        if (x_max == x_max) {
            printf("\nLocation of maximum deflection within [0, d] (in meters): %.6f\n", x_max);
            printf("Maximum deflection: %.6f\n", y(x_max, P, L, d, E, I));
        } else {
            printf("\nCalculation failed.\n");
        }
    } else {
        printf("Invalid choice.\n");
    }
}

int main() {
    calculate_deflection();
    return 0;
}
