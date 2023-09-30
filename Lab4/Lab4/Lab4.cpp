#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const size_t outwidth = 20;
const double round_error = pow(10, -6);

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


// N - number of experiments, m - power of approximating polynomial
template<typename T>
std::vector<T> minimal_square_method(const std::vector<T>& x, std::vector<T>& y, const size_t& M) {
    size_t m = M + 1;
    size_t N = x.size();
    std::vector<T> POWERX;
    for (size_t i = 1; i < 2 * m; ++i) {
        T staff = 0;
        for (size_t j = 0; j < N; ++j) {
            staff += pow(x.at(j), i);
        }
        POWERX.push_back(staff);
    }
    print_vector(POWERX);
    std::vector<std::vector<T>> SUMX(m, std::vector<T>(m, 0));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < m; ++j) {
            if (!i && !j) {
                SUMX.at(0).at(0) = N;
            }
            else {
                SUMX.at(i).at(j) = POWERX.at(i + j - 1);
            }
        }
    }
    print_matrix(SUMX);
    std::vector<T> right_part;
    for (size_t i = 0; i < m; ++i) {
        T staff = 0;
        for (size_t j = 0; j < N; ++j) {
            staff += y.at(j) * pow(x.at(j), i);
        }
        right_part.push_back(staff);
    }
    print_vector(right_part);
    std::vector<std::vector<T>> rectangle = fill_rectangle_matrix(SUMX, right_part);
    std::vector<T> answers = solve_system_with_gauss_method(rectangle);
    T remaining_variance = 0;
    for (size_t i = 0; i < N; ++i) {
        T staff = y.at(i);
        for (size_t j = 0; j < m; ++j) {
            staff -= answers.at(j)*pow(x.at(i), j);
        }
        remaining_variance += pow(staff, 2);
    }
    std::cout << sqrt(remaining_variance / (N - m - 1)) << std::endl;
    return answers;
}

int main()
{
    std::vector<double> x = {0.164, 0.328, 0.656, 0.984, 1.312, 1.640};
    // linearization
    size_t size = x.size();
    for (size_t i = 0; i < size; ++i) {
        x.at(i) = 1.0/x.at(i);
    }
    std::vector<double> y = {0.448, 0.432, 0.421, 0.417, 0.414, 0.412};
    print_vector(minimal_square_method(x, y, 2));
    system("pause");
    return 0;
}