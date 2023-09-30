#include<iostream>
#include<vector>
#include<iomanip>
#include<cmath>
#include<algorithm>
#include<functional>
#include<string>

const size_t outwidth = 25;
const double round_error = pow(10, -9);

// first task

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
                std::swap(matrix.at(i), matrix.at(i-1));
            }
            else {
                std::swap(matrix.at(i), matrix.at(i+1));
            }
        }
    }
}

template <typename T>
std::vector<size_t> find_maximum_index(const std::vector<std::vector<T>>& matrix) {
    std::vector<size_t> answer(2, 0);
    int rows = matrix.size();
    int columns = matrix.at(0).size();
    T max_element = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            if (matrix.at(i).at(j) > max_element) {
                max_element = matrix.at(i).at(j);
                answer.at(0) = i;
                answer.at(1) = j;
            }
        }
    }
    return answer;
}

template<typename T>
void swap_columns(std::vector<std::vector<T>>& matrix, const size_t& column1, const size_t& column2) {
    size_t size = matrix.size();
    for (size_t i = 0; i < size; ++i) {
        std::swap(matrix.at(i).at(column1), matrix.at(i).at(column2));
    }
}

template<typename T>
std::string made_max_element_first(std::vector<std::vector<T>>& matrix, std::vector<T>& bias) {
    std::vector<size_t> index_of_max = find_maximum_index(matrix);
    size_t max_index = matrix.size() - 1;
    std::swap(matrix.at(0), matrix.at(index_of_max.at(0)));
    std::swap(bias.at(0), bias.at(index_of_max.at(0)));
    swap_columns(matrix, 0, index_of_max.at(1));
    std::string answer = std::to_string(index_of_max.at(1));
    for (size_t i = 0; i <= max_index; ++i) {
        if (i != index_of_max.at(1)) {
            answer += " ";
            answer += std::to_string(i);
        }
    }
    return answer;
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
        for (size_t j = i+1; j < heigth; ++j) {
            T aji = solvation.at(j).at(i);
            for (size_t k = i+1; k < width; ++k) {
                solvation.at(j).at(k) -= solvation.at(i).at(k) * aji;
            }
            // made zerolike triangle bottom in matrix
            solvation.at(j).at(i) = 0;
        }
    }
    // backward pass in gauss method
    for (int i = width - 2; i > 0; --i) {
        // perform step subtract one equation to next
        for (int j = i-1; j >= 0; --j) {
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

template<typename T>
std::vector<std::vector<T>> get_cholesky_matrix(const std::vector<std::vector<T>>& matrix) {
    size_t matrix_size = matrix.size();
    // fill two-dimensional array with zeroes
    std::vector<std::vector<T>> answer(matrix_size, std::vector<T>(matrix_size, 0));
    // Decomposing a matrix into Lower Triangular
    for (size_t i = 0; i < matrix_size; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            T sum = 0;
            // summation for diagonals 
            if (j == i) {
                for (size_t k = 0; k < j; ++k) {
                    sum += pow(answer.at(j).at(k), 2);
                }
                answer.at(j).at(j) = sqrt(matrix.at(j).at(j) - sum);
            }
            else {
                // Evaluating L(i, j) using L(j, j)
                for (size_t k = 0; k < j; ++k) {
                    sum += (answer.at(i).at(k) * answer.at(j).at(k));
                }
                answer.at(i).at(j) = (matrix.at(i).at(j) - sum) / answer.at(j).at(j);
            }
        }
    }
    return answer;
}

// can get iagonal elements only in square matrix
template<typename T>
std::vector<T> get_diagonal_elements(const std::vector<std::vector<T>>& matrix) {
    size_t matrix_size = matrix.size();
    std::vector<T> answer;
    for (size_t i = 0; i < matrix_size; ++i) {
        answer.push_back(matrix.at(i).at(i));
    }
    return answer;
}

template <typename T>
std::vector <std::vector<T>> multiply(const std::vector <std::vector<T>>& a, const std::vector <std::vector<T>>& b)
{
    const size_t n = a.size();     // a rows
    const size_t m = a.at(0).size();  // a cols
    const size_t p = b.at(0).size();  // b cols

    std::vector <std::vector<T>> c(n, std::vector<T>(p, 0));
    for (size_t j = 0; j < p; ++j) {
        for (size_t k = 0; k < m; ++k) {
            for (size_t i = 0; i < n; ++i) {
                c.at(i).at(j) += a.at(i).at(k) * b.at(k).at(j);
            }
        }
    }
    return c;
}

template<typename T>
std::vector<T> find_residual_vector(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
    std::vector<std::vector<T>> rectangle = fill_rectangle_matrix(A, b);
    const size_t height = rectangle.size();
    std::vector<T> x_ = solve_system_with_gauss_method(rectangle);

    std::vector <std::vector<T>> x(height, std::vector<T>(1, 0));
    for (size_t i = 0; i < height; ++i) {
        x.at(i).at(0) = x_.at(i);
    }

    std::vector<std::vector<T>> A_mul_x = multiply(A, x);
    print_matrix(A_mul_x);

    std::vector<T> answer(height, 0);
    for (size_t i = 0; i < height; ++i) {
        answer.at(i) = (A_mul_x.at(i).at(0) - b.at(i));
    }
    return answer;
}

template<typename T>
T find_norm(const std::vector<T>& answer) {
    size_t answer_size = answer.size();
    T max = abs(answer.at(0));
    for (size_t i = 1; i < answer_size; ++i) {
        if (max < abs(answer.at(i))) {
            max = abs(answer.at(i));
        }
    }
    return max;
}

template<typename T>
double find_relative_error(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
    std::vector<std::vector<T>> rectangle = fill_rectangle_matrix(A, b);
    std::vector<T> x = solve_system_with_gauss_method(rectangle);
    size_t height = A.size();
    // made vector look like tensor
    std::vector<std::vector<T>> x_(height, std::vector<T>(1, 0));
    for (size_t i = 0; i < height; ++i) {
        x_.at(i).at(0) = x.at(i);
    }
    // multiply this vector to matrix A
    std::vector<std::vector<T>> A_mul_x_ = multiply(A, x_);
    // change result to normal vector
    std::vector<T> A_mul_x(height);
    for (size_t i = 0; i < height; ++i) {
        A_mul_x.at(i) = A_mul_x_.at(i).at(0);
    }
    // find answer in A*__X = A*_X
    std::vector<std::vector<T>> rectangle_1 = fill_rectangle_matrix(A, A_mul_x);
    std::vector<T> final_x = solve_system_with_gauss_method(rectangle_1);
    // find max (0<=i<=n) |__Xi - _Xi|
    T numerator = 0;
    for (size_t i = 0; i < height; ++i) {
        if (abs(final_x.at(i) - x.at(i)) > numerator) {
            numerator = abs(final_x.at(i) - x.at(i));
        }
    }
    // find maximum of first x
    T descriminator = *max_element(x.begin(), x.end());
    return numerator / descriminator;
}

/*template<typename T>
std::vector<T> solve_system_with_LDLT_method(const std::vector<std::vector<T>> matrix) {

}*/

int main() {
    std::vector<std::vector<double>> matrix = { {2.50, -3.00, 4.60},
                                                {-3.50, 2.60, 1.50},
                                                {-6.50, -3.50, 7.30} };
    std::vector<double> bias = { -1.05, -14.46, -17.73 };
    std::cout << "The new indexes are : \n" << made_max_element_first(matrix, bias) << std::endl;
    print_matrix(matrix);
    std::vector<std::vector<double>> rectangle = fill_rectangle_matrix(matrix, bias);
    print_matrix(rectangle);
    std::cout << std::endl << std::endl;
    print_vector(solve_system_with_gauss_method(rectangle, 2));
    std::vector<double> residual_vector = find_residual_vector(matrix, bias);
    print_vector(residual_vector);
    std::cout << find_norm(residual_vector) << std::endl;
    std::cout << find_relative_error(matrix, bias) << std::endl;
    system("pause");
    return 0;
}
