#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

typedef double precision_type;
const size_t x_steps_one_segment = 15;
const size_t x_nodes_number_all = x_steps_one_segment * 3 + 1;
const size_t y_steps_one_segment = 12;
const size_t y_nodes_number_all = y_steps_one_segment * 3 + 1;
const size_t max_repetitions = 3000;
const precision_type epsilon = 1.e-5;
const precision_type itteration_parameter = 1.9f;
const precision_type temperature1 = 5.0f;
const precision_type temperature2 = 15.0f;
const precision_type grid_step_x = 0.3f;
const precision_type grid_step_y = 0.2f;

void maxpvr(precision_type& t1, precision_type& delta, precision_type& max_delta)
{
    precision_type d = abs(delta) / abs(t1);
    if (d > max_delta) {
        max_delta = d;
    }
}

int main()
{
    std::ofstream out_dT("dT.txt", std::ios_base::out | std::ios_base::trunc);
    int x1, x2, x3, y1, y2, y3;
    int k = 0;
    precision_type h = grid_step_x, r = grid_step_y;
    std::vector<std::vector<precision_type>> T(x_nodes_number_all, std::vector<precision_type>(y_nodes_number_all, 0));
    precision_type lambda = itteration_parameter;
    int prz = 1;
    int nT = 0;
    precision_type alpha1 = -h / r;
    precision_type alpha2 = -r / h;
    precision_type alpha3 = alpha2 * 0.5f;
    precision_type alpha4 = alpha1 * 0.5f;
    precision_type gamma1 = -2.f * (alpha1 + alpha2);
    precision_type gamma2 = -1.5f * (alpha1 + alpha2);
    precision_type gamma3 = -(alpha1 + alpha2);
    precision_type gamma4 = -(alpha3 + alpha4);
    x1 = x_steps_one_segment;
    x2 = x1 * 2;
    x3 = x1 * 3;
    y1 = y_steps_one_segment;
    y2 = y1 * 2;
    y3 = y1 * 3;
    precision_type t0, t1, delta, max_delta = 0.0f;
    while (k < max_repetitions && prz == 1) {
        k++;

        // № 9
        //
        // went into all segments and skip that zones that we are not facing
        //  y, j
        //  -------K1
        //  _______________________
        // |       |       |       |
        // |  +++  |  +++  |  +++  |
        // |_______|_______|_______|
        // |       |       |       |
        // |  ---  |  +++  |  +++  |
        // |_______|_______|_______|
        // |       |       |       |
        // |  ---  |  ---  |  +++  |
        // |_______|_______|_______|   x, i
        //                  -------K2

        // Fill start temperature in first and last pins
        for (int i = 0; i <= x1; ++i) {
            T.at(i).at(y3) = temperature1;
        }
        for (int i = x2; i <= x3; ++i) {
            T.at(i).at(0) = temperature2;
        }
        for (int x = 0; x <= x3; ++x) {
            for (int y = 0; y < y3; ++y) {
                t0 = T.at(x).at(y);
                precision_type tx = 0;
                bool state = false;
                if (x == 0 && y > y2 && y < y3) {
                    tx = -(alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1)) / gamma3;
                    state = true;
                }
                else if (x == 0 && y == y2) {
                    tx = -(alpha3 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1)) / gamma4;
                    state = true;
                }
                else if (x <= x1 && y == y2) {
                    tx = -(alpha3 * T.at(x - 1).at(y) + alpha3 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1)) / gamma3;
                    state = true;
                }
                else if (x == x1 && y == y2) {
                    tx = -(alpha4 * T.at(x).at(y - 1) + alpha3 * T.at(x - 1).at(y) + alpha2 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1)) / gamma2;
                    state = true;
                }
                else if (x == x1 && y > y1 && y < y2) {
                    tx = -(alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1)) / gamma3;
                    state = true;
                }
                else if (x == x1 && y == y1) {
                    tx = -(alpha3 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1)) / gamma4;
                    state = true;
                }
                else if (y == y1 && x > x1 && x < x2) {
                    tx = -(alpha3 * T.at(x - 1).at(y) + alpha3 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1)) / gamma3;
                    state = true;
                }
                else if (x == x2 && y == y1) {
                    tx = -(alpha4 * T.at(x).at(y - 1) + alpha3 * T.at(x - 1).at(y) + alpha2 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1)) / gamma2;
                    state = true;
                }
                else if (x == x2 && y > y1) {
                    tx = -(alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1)) / gamma3;
                    state = true;
                }
                else if (x == x3 && y > 0 && y < y3) {
                    tx = -(alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x - 1).at(y) + alpha4 * T.at(x).at(y + 1)) / gamma3;
                    state = true;
                }
                else if (x == x3 && y == y3) {
                    tx = -(alpha4 * T.at(x).at(y - 1) + alpha3 * T.at(x - 1).at(y)) / gamma4;
                    state = true;
                }
                else if (y == y3 && x > x1 && x != x3) {
                    tx = -(alpha3 * T.at(x - 1).at(y) + alpha3 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1)) / gamma3;
                    state = true;
                }
                else if ((x > 0 && x < x3 && y > y2 && y < y3) || (x > x1 && x < x3 && y > y1 && y < y2) || (x > x2 && x < x3 && y > 0 && y < y1)) {
                    tx = -(alpha1 * T.at(x).at(y - 1) + alpha2 * T.at(x - 1).at(y) + alpha2 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1)) / gamma1;
                    state = true;
                }
                if (state) {
                    delta = lambda * (tx - t0);
                    t1 = t0 + delta;
                    T.at(x).at(y) = t1;
                    maxpvr(t1, delta, max_delta);
                }
            }
        }
        nT++;
        precision_type w = max_delta;
        out_dT << w << " ";
        if (max_delta < epsilon) {
            prz = 0;
            max_delta = 0.0f;
        }
    }
    out_dT.close();
    std::ofstream out_nT("nT.txt", std::ios_base::out | std::ios_base::trunc);
    out_nT << nT;
    out_nT.close();
    std::ofstream out_field("field.txt", std::ios_base::out | std::ios_base::trunc);
    for (int i = 0; i < x_nodes_number_all; ++i) {
        for (int j = 0; j < y_nodes_number_all; ++j) {
            out_field << T.at(i).at(j) << " ";
        }
        out_field << "\n";
    }
    out_field.close();
    std::ofstream out_parameters("parameters.txt", std::ios_base::out | std::ios_base::trunc);
    out_parameters << x_nodes_number_all << " " << y_nodes_number_all;
    out_parameters.close();
    system("pause");
    return 0;
}