#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

typedef double precision_type;
const size_t x_steps_one_segment = 15;
const size_t x_nodes_number_all = x_steps_one_segment * 3 + 1;
const size_t y_steps_one_segment = 12;
const size_t y_nodes_number_all = y_steps_one_segment * 3 + 1;
const size_t max_repetitions = 30000;
const size_t time_steps_reload = 1000;
const precision_type termal_diffusity = 1.9f;
const precision_type temperature1 = 5.0f;
const precision_type temperature2 = 15.0f;
const precision_type grid_step_x = 0.3f;
const precision_type grid_step_y = 0.2f;

int main()
{
    int x1, x2, x3, y1, y2, y3;
    precision_type h = grid_step_x, r = grid_step_y;
    std::vector<std::vector<precision_type>> T(x_nodes_number_all, std::vector<precision_type>(y_nodes_number_all, 0));
    double tau = 0.25f * pow(std::min(h, r), 2) / termal_diffusity;
    double alpha1 = -h / r;
    double alpha2 = -r / h;
    double alpha3 = 0.5f * alpha2;
    double alpha4 = 0.5f * alpha1;
    double beta1 = termal_diffusity * tau / (h * r);
    double beta2 = 2.0f * beta1;
    double beta3 = 4.0f * beta1;
    double beta4 = 4.0f * beta1 / 3;
    double gamma1 = -2.f * (alpha1 + alpha2);
    double gamma2 = -1.5 * (alpha1 + alpha2);
    double gamma3 = -(alpha1 + alpha2);
    double gamma4 = -(alpha3 + alpha4);
    x1 = x_steps_one_segment;
    x2 = x1 * 2;
    x3 = x1 * 3;
    y1 = y_steps_one_segment;
    y2 = y1 * 2;
    y3 = y1 * 3;
    size_t current_graph = 0;
    precision_type t0, t1, delta, max_delta = 0.0f;
    for (size_t k = 0; k < max_repetitions; ++k) {
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
        for (int i = x2 - 1; i <= x3; ++i) {
            T.at(i).at(0) = temperature2;
        }
        for (int x = 0; x <= x3; ++x) {
            for (int y = 0; y < y3; ++y) {
                // +
                if (x == 0 && y > y2 && y < y3) {
                    T.at(x).at(y) -= beta2 * (alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1) + gamma3 * T.at(x).at(y));
                }
                else if (x == 0 && y == y2) {
                    T.at(x).at(y) -=  beta1 * (alpha3 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1) + gamma4 * T.at(x).at(y));
                }
                // +
                else if (x <= x1 && y == y2) {
                    T.at(x).at(y) -= beta2 * (alpha3 * T.at(x - 1).at(y) + alpha3 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1) + gamma3 * T.at(x).at(y));
                }
                else if (x == x1 && y == y2) {
                    T.at(x).at(y) -= beta3 * (alpha4 * T.at(x).at(y - 1) + alpha3 * T.at(x - 1).at(y) + alpha2 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1) + gamma2 * T.at(x).at(y));
                }
                // +
                else if (y == y3 && x > x1 && x != x3) {
                    T.at(x).at(y) -= beta2 * (alpha3 * T.at(x - 1).at(y) + alpha3 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1) + gamma3 * T.at(x).at(y));
                }
                // +
                else if (x == x3 && y > 0 && y < y3) {
                    T.at(x).at(y) -= beta2 * (alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x - 1).at(y) + alpha4 * T.at(x).at(y + 1) + gamma3 * T.at(x).at(y));
                }
                // +
                else if (x == x2 && y > y1) {
                    T.at(x).at(y) -= beta2 * (alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1) + gamma3 * T.at(x).at(y));
                }
                // +
                else if (x == x1 && y > y1 && y < y2) {
                    T.at(x).at(y) -= beta2 * (alpha4 * T.at(x).at(y - 1) + alpha2 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1) + gamma3 * T.at(x).at(y));
                }
                else if (x == x1 && y == y1) {
                    T.at(x).at(y) -= beta1 * (alpha3 * T.at(x + 1).at(y) + alpha4 * T.at(x).at(y + 1) + gamma4 * T.at(x).at(y));
                }
                else if (y == y1 && x > x1 && x < x2) {
                    T.at(x).at(y) -= beta2 * (alpha3 * T.at(x - 1).at(y) + alpha3 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1) + gamma3 * T.at(x).at(y));
                }
                else if (x == x2 && y == y1) {
                    T.at(x).at(y) -= beta3 * (alpha4 * T.at(x).at(y - 1) + alpha3 * T.at(x - 1).at(y) + alpha2 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1) + gamma2 * T.at(x).at(y));
                }
                else if (x == x3 && y == y3) {
                    T.at(x).at(y) -= beta1 * (alpha4 * T.at(x).at(y - 1) + alpha3 * T.at(x - 1).at(y) + gamma4 * T.at(x).at(y));
                }
                else if ((x > 0 && x < x3 && y > y2 && y < y3) || (x > x1 && x < x3 && y > y1 && y < y2) || (x > x2 && x < x3 && y > 0 && y < y1)) {
                    T.at(x).at(y) -= beta4 * (alpha1 * T.at(x).at(y - 1) + alpha2 * T.at(x - 1).at(y) + alpha2 * T.at(x + 1).at(y) + alpha1 * T.at(x).at(y + 1) + gamma1 * T.at(x).at(y));
                }
            }
        }
        if (k % time_steps_reload == 0 || k == 0) {
            std::ofstream out_dT(std::string("./graphs/T") + std::to_string(current_graph) + std::string(".txt"), std::ios_base::out | std::ios_base::trunc);
            for (int i = 0; i < x_nodes_number_all; ++i) {
                for (int j = 0; j < y_nodes_number_all; ++j) {
                    out_dT << T.at(i).at(j) << " ";
                }
                out_dT << "\n";
            }
            out_dT.close();
            current_graph++;
        }
    }
    system("pause");
    return 0;
}