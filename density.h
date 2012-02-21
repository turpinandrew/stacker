
#define DENSITY_DATA_LEN 35

#define FEQUAL(_x,_y) (fabs((_x) - (_y)) < 0.0001)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double interp(double x1, double y1, double x2, double y2, double x);
double get_density_axis(double *xs, double *ys, double x);
double find_density(double x, double y);
