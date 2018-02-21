#include <matrix.hpp>

using namespace chislaki;


int main() {
    std::vector<int> arr = { 1, 2, 3, 4, 5, 6 };
    matrix<int> m(2, 3, arr);
    std::cout << m;
    std::cout << "max = " << m.col_max_value(0) << std::endl;
}
