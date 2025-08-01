#include <iostream>

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

void print(Kokkos::View<int**, Kokkos::SharedSpace> &a)
{
    std::cout << "Contents: ";
    for (size_t i = 0; i < a.extent(0); i++)
    {
        for (size_t j = 0; j < a.extent(1); j++)
        {
            std::cout << a(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print(Kokkos::View<int*, Kokkos::SharedSpace> &a)
{
    std::cout << "Contents: ";
    for (size_t i = 0; i < a.size(); i++)
    {
        std::cout << a(i) << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);
    Kokkos::View<int**, Kokkos::SharedSpace> nums("nums",10,10);
    for (size_t i = 0; i < 100; i++)
    {
        nums((i/10), (i%10)) = i;
    }
    print(nums);
    Kokkos::View<int*, Kokkos::SharedSpace> num(nums, 4, Kokkos::ALL);
    print(num);
    num(6) = 200;
    print(num);
    print(nums);
    std::cout << nums.rank() << std::endl;
    std::cout << num.rank() << std::endl;
    int* it = &num(0);
    it += 4;
    std::cout << *it << std::endl;
    *it = 553;
    std::cout << *it << std::endl;
    print(nums);
    print(num);
    return 0;
}
