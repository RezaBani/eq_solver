#include <iostream>
#include <cmath>

using namespace std;

// for fixed point method input is g(x) that is, x = g(x) so we didn't use it here
int baniv_sign(double x); // for baniv method (return -1 for x = 0)
int sign(double x);
double equation(double x);
double derivative(double x);
double newton_raphson(double (*eq)(double), double (*de)(double), double initial_val, double rhs = 0.0, double max_err = 0.001, size_t max_trial = 10000);
double secant_method(double (*eq)(double), double initial_val1 , double initial_val2, double rhs = 0.0, double max_err = 0.001, size_t max_trial = 10000);
double baniv_method(double (*eq)(double), bool forward, double initial_val, double rhs = 0.0, double step = 0.001, size_t max_trial = 10000);
double bisection_method(double (*eq)(double), double initial_val1 , double initial_val2, double rhs = 0.0, double max_err = 0.001, size_t max_trial = 10000);
double false_position_method(double (*eq)(double), double initial_val1 , double initial_val2, double rhs = 0.0, double max_err = 0.001, size_t max_trial = 10000);
double ridders_method(double (*eq)(double), double initial_val1 , double initial_val2, double rhs = 0.0, double max_err = 0.001, size_t max_trial = 10000);
double fixed_point_method(double (*eq)(double), double initial_val, double rhs = 0.0, double max_err = 0.001, size_t max_trial = 10000);


int main() {
    const double x0 = 1.01;
    const double x1 = 20.0;
    const double y = 26.37;
    double result = 0.0;

    cout << "Newton-Raphson Method (X0 = " << x0 << "):" << endl;
    result = newton_raphson(equation, derivative, x0, y);
    cout << "X = " << result << endl << endl;

    cout << "Newton-Raphson Method (X1 = " << x1 << "):" << endl;
    result = newton_raphson(equation, derivative, x1, y);
    cout << "X = " << result << endl << endl;

    cout << "Secant Method:" << endl;
    result = secant_method(equation, x0, x1, y);
    cout << "X = " << result << endl << endl;

    cout << "Baniv Method Forward:" << endl;
    result = baniv_method(equation, true, x0, y);
    cout << "X = " << result << endl << endl;

    cout << "Baniv Method Backward:" << endl;
    result = baniv_method(equation, false, x1, y);
    cout << "X = " << result << endl << endl;

    cout << "Bisection Method:" << endl;
    result = bisection_method(equation, x0, x1, y);
    cout << "X = " << result << endl << endl;

    cout << "False Position Method:" << endl;
    result = false_position_method(equation, x0, x1, y);
    cout << "X = " << result << endl << endl;

    cout << "Ridders' Method:" << endl;
    result = ridders_method(equation, x0, x1, y);
    cout << "X = " << result << endl << endl;

    return 0;
}

int sign(double x){
    if (x < 0){
        return -1;
    } else {
        return 1;
    }
}

int baniv_sign(double x){
    if (x > 0){
        return 1;
    } else {
        return -1;
    }
}

double equation(double x){
    const double gamma = 1.4;
    const double pi = 3.141592653589793;
    const double a = sqrt((gamma + 1)/(gamma - 1));
    const double b = (gamma - 1)/(gamma + 1);
    return a*atan(sqrt(b*(x*x-1)))*180/pi-atan(sqrt(x*x-1))*180/pi;
}

double derivative(double x){
    const double gamma = 1.4;
    const double pi = 3.141592653589793;
    const double a = sqrt((gamma + 1)/(gamma - 1));
    const double b = (gamma - 1)/(gamma + 1);
    return 180*a*x*b/(pi*(1+b*(x*x-1))*(sqrt(b*(x*x-1)))) + 180/(pi*x*sqrt(x*x-1));
}

double newton_raphson(double (*eq)(double), double (*de)(double), double initial_val, double rhs, double max_err, size_t max_trial){
    size_t trial = 0;
    double x = initial_val;
    double y = 0;
    double dy = 0;
    do {
        y = eq(x) - rhs;
        dy = de(x);
        x = x - y/dy;
        trial++;
    } while ((abs(y) > max_err) && (trial < max_trial));
    cout << "Trials: " << trial << endl;
    cout << "Error: " << abs(y) << endl;
    return x;
}

double secant_method(double (*eq)(double), double initial_val1, double initial_val2, double rhs, double max_err, size_t max_trial){
    size_t trial = 0;
    double x0 = initial_val1;
    double x1 = initial_val2;
    double x2 = 0;
    double y0 = 0;
    double y1 = 0;
    do {
        y0 = eq(x0)-rhs;
        y1 = eq(x1)-rhs;
        x2 = x1 - y1*((x1-x0)/(y1-y0));
        x0 = x1;
        x1 = x2;
        trial++;
    } while ((abs(y1) > max_err) && (trial < max_trial));
    cout << "Trials: " << trial << endl;
    cout << "Error: " << abs(y1) << endl;
    return x1;
}

double baniv_method(double (*eq)(double), bool forward, double initial_val, double rhs, double step, size_t max_trial){
    size_t trial = 0;
    double x = initial_val;
    double y = 0;
    const int sgn = baniv_sign(forward);
    x -= sgn*step;
    do {
        x += sgn*step;
        y = eq(x) - rhs;
        trial++;
    } while ((sgn*y < 0) && (trial < max_trial));
    cout << "Trials: " << trial << endl;
    cout << "Error: " << abs(y) << endl;
    return x;
}

double bisection_method(double (*eq)(double), double initial_val1 , double initial_val2, double rhs, double max_err, size_t max_trial){
    size_t trial = 0;
    double x0 = initial_val1;
    double x1 = initial_val2;
    double xm = 0;
    double y0 = 0;
    double y1 = 0;
    double ym = 0;
    do {
        xm = (x0 + x1) / 2;
        y0 = eq(x0) - rhs;
        y1 = eq(x1) - rhs;
        ym = eq(xm) - rhs;
        if (signbit(ym) == signbit(y0)){
            x0 = xm;
        } else {
            x1 = xm;
        }
        trial++;
    } while ((abs(ym) > max_err) && (trial < max_trial));
    cout << "Trials: " << trial << endl;
    cout << "Error: " << abs(ym) << endl;
    return xm;
}

double false_position_method(double (*eq)(double), double initial_val1 , double initial_val2, double rhs, double max_err, size_t max_trial){
    size_t trial = 0;
    double x0 = initial_val1;
    double x1 = initial_val2;
    double y0 = 0;
    double y1 = 0;
    double xm = 0;
    double ym = 0;
    do {
        y0 = eq(x0) - rhs;
        y1 = eq(x1) - rhs;
        xm = (x0 * y1 - x1 * y0) / (y1 - y0);
        ym = eq(xm) - rhs;
        if (signbit(ym) == signbit(y0)){
            x0 = xm;
        } else {
            x1 = xm;
        }
        trial++;
    } while ((abs(ym) > max_err) && (trial < max_trial));
    cout << "Trials: " << trial << endl;
    cout << "Error: " << abs(ym) << endl;
    return xm;
}

double ridders_method(double (*eq)(double), double initial_val1 , double initial_val2, double rhs, double max_err, size_t max_trial){
    size_t trial = 0;
    double x0 = initial_val1;
    double x2 = initial_val2;
    double y0 = 0;
    double y2 = 0;
    double x1 = 0;
    double y1 = 0;
    double x3 = 0;
    double y3 = 0;
    do {
        y0 = eq(x0) - rhs;
        y2 = eq(x2) - rhs;
        x1 = (x0 + x2) / 2;
        y1 = eq(x1) - rhs;
        x3 = x1 + (x1 - x0)*(sign(y0)*y1)/(sqrt(y1*y1 - y0*y2));
        y3 = eq(x3) - rhs;
        if (signbit(y1) != signbit(y3)){
            x0 = x1;
            x2 = x3;
        } else {
            if (signbit(y0) != signbit(y3)){
                x2 = x3;
            } else {
                x0 = x3;
            }
        }
        trial++;
    } while ((abs(y3) > max_err) && (trial < max_trial));
    cout << "Trials: " << trial << endl;
    cout << "Error: " << abs(y3) << endl;
    return x3;
}

double fixed_point_method(double (*eq)(double), double initial_val, double rhs, double max_err, size_t max_trial){
    size_t trial = 0;
    double x = initial_val;
    double y = 0;
    do {
        y = eq(x) - rhs;
        x = y;
        trial++;
    } while ((abs(y - x) > max_err) && (trial < max_trial));
    cout << "Trials: " << trial << endl;
    cout << "Error: " << abs(y - x) << endl;
    return x;
}