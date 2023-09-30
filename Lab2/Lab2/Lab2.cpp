#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>

const size_t outwidth = 20;
const size_t precision = 10;
const double round_error = 0.000000001;

template<typename T>
std::vector<std::vector<T>> initialization_matrix(size_t& n = 0) {
    if (n) {
        std::cout << "Write amount of elements in square matrix\n";
        std::cin >> n;
    }
    std::vector<std::vector<T>> matrix;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cout << i + 1 << " " << j + 1 << std::endl;
            std::cin >> matrix.at(i).at(j);
        }
    }
    return matrix;
}

template<typename T>
std::vector<T> initialize_vector(size_t& n = 0) {
    if (n) {
        std::cout << "Write amount of elements in vector\n";
        std::cin >> n;
    }
    std::vector<T> const_vector;
    for (size_t i = 0; i < n; ++i) {
        std::cout << i + 1 << std::endl;
        std::cin >> const_vector.at(i);
    }
}

template<typename T>
void print_matrix(const std::vector<std::vector<T>>& matrix) {
    const size_t rows_num = matrix.size();
    const size_t column_num = matrix.at(0).size();
    for (size_t i = 0; i < rows_num; ++i) {
        for (size_t j = 0; j < column_num; ++j) {
            std::cout << std::setw(outwidth) << matrix.at(i).at(j);
        }
        std::cout << std::endl;
    }
}

template<typename T>
void print_vector(const std::vector<T>& const_vector) {
    size_t length = const_vector.size();
    for (size_t i = 0; i < length; ++i) {
        std::cout << std::setw(outwidth) << const_vector.at(i);
    }
    std::cout << std::endl;
}

template<typename T>
std::vector<std::vector<T>> fill_rectangle_matrix(const std::vector<std::vector<T>>& matrix, const std::vector<T>& const_vector) {
    std::vector<std::vector<T>> solvation_matrix;
    const size_t rows_num = matrix.size();
    // filling matrix with square matrix and vector as last column
    for (size_t i = 0; i < rows_num; ++i) {
        std::vector<T> row;
        for (size_t j = 0; j <= rows_num; ++j) {
            if (j < rows_num) {
                row.push_back(matrix.at(i).at(j));
            }
            else {
                row.push_back(const_vector.at(i));
            }
        }
        solvation_matrix.push_back(row);
    }
    return solvation_matrix;
}

template<typename T>
void swap_rows_for_not_null_diagonal(std::vector<std::vector<T>>& matrix) {
    size_t matrix_size = matrix.size();
    for (size_t i = 0; i < matrix_size; ++i) {
        if (abs(matrix.at(i).at(i)) <= round_error) {
            if (i > 0) {
                std::swap(matrix.at(i), matrix.at(i - 1));
            }
            else {
                std::swap(matrix.at(i), matrix.at(i + 1));
            }
        }
    }
}

template<typename T>
std::vector<T> solve_system_with_gauss_method(std::vector<std::vector<T>>& solvation, const int& epsilon = 0) {
    swap_rows_for_not_null_diagonal(solvation);
    const size_t heigth = solvation.size();
    const size_t width = solvation.at(0).size();
    // forward pass in gauss method
    for (size_t i = 0; i < heigth; ++i) {
        // get diagonal element
        T aii = solvation.at(i).at(i);
        // divide row by it's first element (not zero in this algorithm)
        for (size_t j = i; j < width; ++j) {
            solvation.at(i).at(j) /= aii;
        }
        // perform step subtract one equation to next
        for (size_t j = i + 1; j < heigth; ++j) {
            T aji = solvation.at(j).at(i);
            for (size_t k = i + 1; k < width; ++k) {
                solvation.at(j).at(k) -= solvation.at(i).at(k) * aji;
            }
            // made zerolike triangle bottom in matrix
            solvation.at(j).at(i) = 0;
        }
    }
    // backward pass in gauss method
    for (int i = width - 2; i > 0; --i) {
        // perform step subtract one equation to next
        for (int j = i - 1; j >= 0; --j) {
            T aji = solvation.at(j).at(i);
            for (int k = width - 1; k >= i; --k) {
                solvation.at(j).at(k) -= solvation.at(i).at(k) * aji;
            }
        }
    }
    std::vector<T> answer;
    if (epsilon) {
        for (size_t i = 0; i < heigth; ++i) {
            answer.push_back((std::round((solvation.at(i).at(width - 1)) * pow(10, epsilon))) / pow(10, epsilon));
        }
    }
    else {
        for (size_t i = 0; i < heigth; ++i) {
            answer.push_back(solvation.at(i).at(width - 1));
        }
    }
    return answer;
}

// second task

double first_equasion(const std::vector<double>& x) {
    return 1.5 * pow(x.at(0), 3) - pow(x.at(1), 2) - 1;
}

double second_equasion(const std::vector<double>& x) {
    return x.at(0) * pow(x.at(1), 3) - x.at(1) - 4;
}

double diff_first_x1(const std::vector<double>& x) {
    return 4.5 * pow(x.at(0), 2);
}

double diff_first_x2(const std::vector<double>& x) {
    return -2 * x.at(1);
}

double diff_second_x1(const std::vector<double>& x) {
    return pow(x.at(1), 3);
}

double diff_second_x2(const std::vector<double>& x) {
    return 3 * x.at(0) * pow(x.at(1), 2) - 1;
}

std::vector<std::vector<std::function<double(std::vector<double>)>>> jacobian = {
    {diff_first_x1, diff_first_x2},
    {diff_second_x1, diff_second_x2}
};

std::vector<std::vector<double>> find_jacobian_1(const std::vector<double>& x) 
{
    std::vector<std::vector<double>> jacobian_1(x.size(), std::vector<double>(x.size(), 0));
    jacobian_1.at(0).at(0) = diff_first_x1(x);
    jacobian_1.at(0).at(1) = diff_first_x2(x);
    jacobian_1.at(1).at(0) = diff_second_x1(x);
    jacobian_1.at(1).at(1) = diff_second_x2(x);
    return jacobian_1;
}

std::vector<std::vector<double>> find_jacobian_2(const std::vector<double>& x, const double& M) 
{
    std::vector<std::vector<double>> jacobian_1(x.size(), std::vector<double>(x.size(), 0));

    std::vector<double> x_1 = { x.at(0) * (1 + M), x.at(1) };
    std::vector<double> x_2 = { x.at(0), x.at(1) * (1 + M) };

    jacobian_1.at(0).at(0) = (first_equasion(x_1) - first_equasion(x)) / (M*x.at(0));
    jacobian_1.at(0).at(1) = (first_equasion(x_2) - first_equasion(x)) / (M*x.at(1));
    jacobian_1.at(1).at(0) = (second_equasion(x_1) - second_equasion(x)) / (M*x.at(0));
    jacobian_1.at(1).at(1) = (second_equasion(x_2) - second_equasion(x)) / (M*x.at(1));
    return jacobian_1;
}

// dx1
double find_d1(const std::vector<double>& x)
{
    double x1 = 1.5 * pow(x.at(0), 3) - pow(x.at(1), 2) - 1;
    double x2 = x.at(0) * pow(x.at(1), 3) - x.at(1) - 4;
    return (abs(x1) > abs(x2)) ? abs(x1) : abs(x2);
}

// dx2
double find_d2(const std::vector<double>& x, const std::vector<double>& x1)
{
    int n = x.size();
    double max = 0;
    for (size_t i = 0; i < n; ++i) {
        if (abs(x1.at(i)) < 1) {
            if (abs(x1.at(i) - x.at(i)) > max) {
                max = abs(x1.at(i) - x.at(i));
            }
        }
        else {
            if (abs((x1.at(i) - x.at(i)) / x1.at(i)) > max) {
                max = abs((x1.at(i) - x.at(i)) / x1.at(i));
            }
        }
    }
    return max;
}

std::vector<double> solve_system_with_Newtons_method_two_equations(const std::vector<double>& epsilons, const std::vector<double>& start_x, const double& M=0.01, const size_t& number_of_itterations=100)
{
    size_t n = start_x.size();
    std::vector<double> x = start_x;
    size_t k = 1;
    double d1 = 1, d2 = 1;
    std::vector<double> F(2, 0);
    std::vector<std::vector<double>> J;
    std::vector<double> deltax;
    std::vector<double> x1(2, 0);
    while ((d1 > epsilons.at(0) || d2 > epsilons.at(1)) && k < number_of_itterations) {
        //find_residual_vector
        F.at(0) = -first_equasion(x);
        F.at(1) = -second_equasion(x);
        //J = find_jacobian_1(x);
        J = find_jacobian_2(x, M);
        std::vector<std::vector<double>> rectangle = fill_rectangle_matrix(J, F);
        deltax = solve_system_with_gauss_method(rectangle);
        for (size_t i = 0; i < n; ++i) {
            x1.at(i) = x.at(i) + deltax.at(i);
        }
        d1 = find_d1(x);
        d2 = find_d2(x, x1);
        std::cout << std::setw(outwidth) << std::setprecision(precision) << k << "   " << d1 << "   " << d2 << std::endl;
        k++;
        x = x1;
        if (k >= number_of_itterations) {
            std::cout << "ERROR: IER = 2" << std::endl;
        }
    }
    return x;
}

int main() 
{
    std::vector<double> epsilons(2, 0.000000001);
    std::vector<double> x(2, 1);
    print_vector(solve_system_with_Newtons_method_two_equations(epsilons, x, 0.05));
    system("pause");
    return 0;
}