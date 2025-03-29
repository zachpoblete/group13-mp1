#include <math.h>
#include <stdio.h>

void RootFinding() {

#define TOLERANCE 1e-6
#define MAX_ITER 100

  double P, L, E, I, d;

  double deflection(double x) {
    return (P * (L - d) / (4 * E * I)) * (-x * x + (1 / d) * x * x * x);
  }

  double deflection_derivative(double x) {
    return (P * (L - d) / (4 * E * I)) * (-2 * x + (3 / d) * x * x);
  }

  double newton_raphson(double x0) {
    double x = x0;
    for (int i = 0; i < MAX_ITER; i++) {
      double fx = deflection_derivative(x);
      double fpx = (deflection_derivative(x + TOLERANCE) - fx) / TOLERANCE;
      double x_new = x - fx / fpx;
      if (fabs(x_new - x) < TOLERANCE) {
        return x_new;
      }
      x = x_new;
    }
    return x;
  }

  double regula_falsi(double a, double b) {
    double fa = deflection_derivative(a);
    double fb = deflection_derivative(b);
    for (int i = 0; i < MAX_ITER; i++) {
      double c = (a * fb - b * fa) / (fb - fa);
      double fc = deflection_derivative(c);
      if (fabs(fc) < TOLERANCE) {
        return c;
      }
      if (fc * fa < 0) {
        b = c;
        fb = fc;
      } else {
        a = c;
        fa = fc;
      }
    }
    return (a + b) / 2;
  }

  int main() {
    printf("Enter the distance d (0 < d < L): ");
    scanf("%lf", &d);

    printf("Enter the values for P, L, E, I: ");
    scanf("%lf %lf %lf %lf", &P, &L, &E, &I);

    int method;
    printf("Choose method (1 for Newton-Raphson, 2 for Regula Falsi): ");
    scanf("%d", &method);

    double x_max;
    if (method == 1) {
      x_max = newton_raphson(d / 2);
    } else {
      x_max = regula_falsi(0, d);
    }

    double y_max = deflection(x_max);
    double delta_free_end =
        (P * (L - d) / (E * I)) * (-L * L / 4 + L * L * L / (4 * d) -
                                   (L - d) * (L - d) * (3 * L - d) / (12 * d));

    printf("Location of maximum deflection within x = 0 to x = d: %lf\n",
           x_max);
    printf("Maximum deflection within x = 0 to x = d: %lf\n", y_max);
    printf("Maximum deflection at the free end: %lf\n", delta_free_end);
  }
}
