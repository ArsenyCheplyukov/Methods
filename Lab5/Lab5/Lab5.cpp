#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

double function(const double& x) {
    return (1 + x + pow(x, 2))/pow((pow(x, 3) - 1), 0.5);
}

double function2(const double& x, const double& y) {
    return 1 / pow((x + y), 2);
}

double find_definite_integral_with_trapezoid_method(const double& first_span, const double& second_span, const double& accuracy) {
    int n = 2;
    double step1 = (second_span - first_span) / n;
    double Integral = function(first_span) + function(second_span);
    double Integral2;
    double current_x = first_span;
    for (int i = 0; i < n - 1; ++i) {
        current_x += step1;
        Integral += 2 * function(current_x);
    }
    Integral *= step1 / 2;
    do {
        n *= 2;
        step1 = (second_span - first_span) / n;
        Integral2 = Integral;
        Integral = function(first_span) + function(second_span);
        current_x = first_span;
        for (int i = 0; i < n - 1; ++i) {
            current_x += step1;
            Integral += 2 * function(current_x);
        }
        Integral *= step1 / 2;
    } while (abs(Integral - Integral2) > 3 * accuracy);
    double error;
    error = abs((Integral2 - Integral) / (pow(0.5, 2) - 1));
    std::cout << Integral << " (error " << error << ")\n";
    return Integral;
}

double find_definite_integral_with_Simpsons_method(const double& first_span, const double& second_span, const double& accuracy) {
    int n = 4;
    double step1 = (second_span - first_span) / n;
    double Integral = function(first_span) + function(second_span);
    double Integral2;
    double current_x = first_span;
    for (int i = 0; i < n - 1; i += 2) {
        current_x += step1;
        Integral += 4 * function(current_x);
        current_x += step1;
        Integral += 2 * function(current_x);
    }
    Integral *= step1 / 3;
    do {
        n *= 2;
        step1 = (second_span - first_span) / n;
        Integral2 = Integral;
        Integral = function(first_span) + function(second_span);
        current_x = first_span;
        for (int i = 0; i < n - 1; i += 2) {
            current_x += step1;
            Integral += 4 * function(current_x);
            current_x += step1;
            Integral += 2 * function(current_x);
        }
        Integral *= step1 / 3;
    } while (abs(Integral - Integral2) > 15 * accuracy);

    double error = abs((Integral2 - Integral) / (pow(0.5, 4) - 1));

    std::cout << Integral << " (error " << error << ")\n";
    return Integral;
}

double find_definite_integral_with_cube_Simpsons_method(const std::vector<double>& first_span, const std::vector<double>& second_span, const double& accuracy) {
    int n = 4, m = n;
    double step1 = (first_span.at(1) - first_span.at(0)) / (2 * n), step2 = (second_span.at(1) - second_span.at(0)) / (2 * m);
    double integral = 0, integral2;

    std::vector<double> x = {first_span.at(0)};
    std::vector<double> y = {second_span.at(0)};

    for (int i = 1; i <= 2 * (n+2); ++i) {
        x.push_back(first_span.at(0) + i * step1);
        y.push_back(second_span.at(0) + i * step2);
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            integral += function2(x.at(2 * i), y.at(2 * j)) + function2(x.at(2 * i + 2), y.at(2 * j)) + function2(x.at(2 * i + 2), y.at(2 * j + 2)) + function2(x.at(2 * i), y.at(2 * j + 2)) +
                4 * (function2(x.at(2 * i + 1), y.at(2 * j)) + function2(x.at(2 * i + 2), y.at(2 * j + 1)) + function2(x.at(2 * i + 1), y.at(2 * j + 2)) + function2(x.at(2 * i), y.at(2 * j + 1))) +
                16 * function2(x.at(2 * i + 1), y.at(2 * j + 1));
        }
    }
    integral *= step1 * step2 / 9;

    do {
        n *= 2;
        step1 = (first_span.at(1) - first_span.at(0)) / (2 * n);
        step2 = (second_span.at(1) - second_span.at(0)) / (2 * m);
        integral2 = integral;
        integral = 0;
        std::vector<double> x = { first_span.at(0) };
        std::vector<double> y = { second_span.at(0) };

        for (int i = 1; i <= 2 * (n+2); ++i) {
            x.push_back(first_span.at(0) + i * step1);
            y.push_back(second_span.at(0) + i * step2);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                integral += function2(x.at(2 * i), y.at(2 * j)) + function2(x.at(2 * i + 2), y.at(2 * j)) + function2(x.at(2 * i + 2), y.at(2 * j + 2)) + function2(x.at(2 * i), y.at(2 * j + 2)) +
                    4 * (function2(x.at(2 * i + 1), y.at(2 * j)) + function2(x.at(2 * i + 2), y.at(2 * j + 1)) + function2(x.at(2 * i + 1), y.at(2 * j + 2)) + function2(x.at(2 * i), y.at(2 * j + 1))) +
                    16 * function2(x.at(2 * i + 1), y.at(2 * j + 1));
            }
        }
        integral *= step1 * step2 / 9;
    } while (abs(integral - integral2) > 15 * accuracy);

    std::cout << integral;
    return integral;
}

int main()
{
    double first_span = 1.0;
    double second_span = 2.631;
    std::vector<double> first_span_1 = {3.0, 4.0};
    std::vector<double> second_span_1 = {1.0, 2.0};
    double epsilon_1 = pow(10, -4);
    double epsilon_2 = pow(10, -5);
    std::cout << " epsilon 10^(-4):\n";
    std::cout << "Trapezoid method: " << std::endl;
    find_definite_integral_with_trapezoid_method(first_span + epsilon_1, second_span, epsilon_1);
    std::cout << "Simpsons method: "  << std::endl;
    find_definite_integral_with_Simpsons_method(first_span + epsilon_1, second_span, epsilon_1);
    std::cout << "Simpsons cube method: " << std::endl;
    find_definite_integral_with_cube_Simpsons_method(first_span_1, second_span_1, epsilon_1);
    std::cout << "\n epsilon 10^(-5):\n";
    std::cout << "Trapezoid method: " << std::endl;
    find_definite_integral_with_trapezoid_method(first_span + epsilon_2, second_span, epsilon_2);
    std::cout << "Simpsons method: " << std::endl;
    find_definite_integral_with_Simpsons_method(first_span + epsilon_2, second_span, epsilon_2);
    std::cout << "Simpsons cube method: " << std::endl;
    find_definite_integral_with_cube_Simpsons_method(first_span_1, second_span_1, epsilon_2);
    std::cout << std::endl;
    system("pause");
    return 0;
}
