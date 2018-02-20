#include <matrix.hpp>

int main() {
    int arr[] = { 1, 2, 3, 4, 5 , 6};

    auto m1 = chislaki::matrix<int>(4, 5);
    auto m2 = chislaki::matrix<int>(2, 3, arr);
}
