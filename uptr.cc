#include <iostream>
#include <memory>
#include <cstring> 
int main() 
{
    const int size = 10; 
    std::unique_ptr<int[]> fact(new int[size]);
 
    for (int i = 0; i < size; ++i) {
        fact[i] = (i == 0) ? 1 : i * fact[i-1];
    }
 
    for (int i = 0; i < size; ++i) {
        std::cout << i << ": " << fact[i] << '\n';
    }

    double * x = new double[10];
    for(size_t i =0;i<10;++i) x[i]=double(i);

    std::unique_ptr<double[]> ux(x);

    for(size_t i =0;i<10;++i) std::cout << ux[i] << '\n';
    
    //will cause a segfault
    //delete [] x;
    //will also segfault
    //ux.release();
    //for(size_t i =0;i<10;++i) std::cout << ux[i] << '\n';
    //for(size_t i =0;i<10;++i) std::cout << x[i] << '\n';

    char cx[] = "0.12345";

    double huh = std::stod(std::string(cx));//reinterpret_cast<double>(&cx[0]);

    std::cout << huh << '\n';
}
