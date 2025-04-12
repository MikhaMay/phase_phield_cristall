#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <sstream>
#include <cmath>

constexpr double pi = 3.14159265358979323846;

class MyArray{
    int size_ = 0;
    double * data = nullptr;
    double l = 0.0;
    double r = 0.0;
public:
    MyArray(const MyArray &) = delete;
    const MyArray& operator= (const MyArray& ) = delete;

    MyArray() {}

    MyArray(int N, double left, double right) : size_(N), l(left), r(right) {
        data = new double[N];
    }

    void swap (MyArray& arr){  
        std::swap(size_, arr.size_);
        std::swap(data, arr.data);
    }

    double& operator[] (int i) {
        if (i == -1) return data[size_-1]; //   Периодические ГУ
        else if (i == size_) return data[0]; 
        // if (i == -1) return l; //  ГУ Дирихле и Неймана
        // else if (i == size_) return r;
        return data[i];
    }

    ~MyArray() {
        delete[] data;
    }

    void show(const std::string& path) {
        // Открываем файл для записи в бинарном режиме
        std::ofstream outfile(path, std::ios::binary);

        // Записываем размер массива
        outfile.write(reinterpret_cast<const char *>(&size_), sizeof(size_));

        // Записываем данные массива
        outfile.write(reinterpret_cast<const char *>(data), size_ * sizeof(double));
        
        // Закрываем файл
        outfile.close();
    }

    void value_assert(double val) {
        for (int i = 0; i < size_; ++i)
            if (std::abs(data[i]) > val || data[i] != data[i])
                throw std::runtime_error("value assert");
    }

    void set_random_condition(double a, double b) {
        std::random_device rd;  // Источник случайных чисел
        std::mt19937 gen(rd()); // Генератор случайных чисел (Mersenne Twister)
        std::uniform_real_distribution<double> dis(a, b); // Равномерное распределение от 0 до 1
        for (size_t i = 0; i < size_; ++i) { // Заполнение массива случайными числами
            data[i] = dis(gen);
        }
    }

    void set_linear_condition() {
        for (int i = 0; i < size_; ++i)
            data[i] = -1.0 + 2.0 * i / (size_-1);   // линейные НУ
    }

    void set_step_condition() {
        for (int i = 0; i < size_/2; ++i)
            data[i] = -1.0;                     // ну в виде ступеньки
        for (int i = size_/2; i < size_; ++i)
            data[i] = 1.0;
    }

    void set_double_condition(const double& avg) {
        std::ifstream file("bin_data/for_initial.bin", std::ios::binary);
        int tmp_N;
        file.read(reinterpret_cast<char*>(&tmp_N), sizeof(tmp_N));
        double * new_data = new double[tmp_N];
        file.read(reinterpret_cast<char*>(new_data), tmp_N*sizeof(double));
        for (int i = 0; i < tmp_N; ++i)
            data[i] = avg + new_data[i];
        for (int i = tmp_N; i < size_; ++i)
            data[i] = avg;
        delete[] new_data;
    }

    void set_sinusoidal_condition(double, double, double);

    void set_book_condition();

    void do_laplacian(MyArray & arr);
};



//--------------------Настройка параметров--------------------

constexpr double L = 80;
constexpr int N = 200;
constexpr double h = L / N;
constexpr double dt = 1e-4;
constexpr int T_steps = 1'000'000;
constexpr double eps = 0.4;
constexpr int moment_in_frame = 10'000;
constexpr double total_time = dt*T_steps; 

//------------------------------------------------------------
void write_params(const char * path){
    std::ofstream os(path);
    os << T_steps<<std::endl;
    os << moment_in_frame<<std::endl;
    os << eps<<std::endl;
}

void MyArray::set_sinusoidal_condition(double q, double amplitude, double phase = 0.0) {
    for (int i = 0; i < size_; ++i){
        data[i] = amplitude*sin(phase);
        phase += q*h;
    }
    if (data[0] != amplitude*sin(phase))
        std::cout << "Conditions is not periodical";
}

void MyArray::set_book_condition() {
    double phi_0 = -0.275;
    double A = 0.5;
    for (int i = 0; i <= N/2; ++i)
        data[i] = phi_0 + 2*A/L * i*h - A/2;
    for (int i = N/2+1; i < N; ++i)
        data[i] = phi_0 - 2*A/L * i*h + 3*A/2;
}

void MyArray::do_laplacian(MyArray & arr){
    for (int i = 0; i < size_; ++i)
        data[i] = (arr[i-1] - 2*arr[i] + arr[i+1]) / (h * h);
}

int main() {
    write_params("bin_data/params.txt");
    MyArray phi(N, 0.0, 0.0); // первый параметр - размер массива - второй и третий параметры - значения на границах
    MyArray phi_next(N, 0.0, 0.0);
    MyArray hi(N, 0.0, 0.0);
    MyArray tmp(N, 0.0, 0.0);
    MyArray energies((T_steps/moment_in_frame) + 1, 0.0, 0.0); // запись энергии в разные моменты времени

    // phi.set_sinusoidal_condition(1.2, 0.4, pi/2);
    // phi.set_book_condition();
    //phi.set_random_condition(-0.2, 0.2);
    phi.set_double_condition(0.45);

    for (int step = 0; step <= T_steps; ++step){

        //======= рассчет следующего слоя =======

        for (int i = 0; i < N; ++i)
            hi[i] = phi[i] + (phi[i-1] - 2 * phi[i] + phi[i+1])/(h*h);

        for (int i = 0; i < N; ++i)
            tmp[i] = (-eps*phi[i] + hi[i] + (hi[i-1] - 2*hi[i] + hi[i+1])/(h*h) + phi[i]*phi[i]*phi[i]);

        for (int i = 0; i < N; ++i)
            phi_next[i] = phi[i] + dt*(tmp[i-1] - 2*tmp[i] + tmp[i+1])/(h*h);

        //=======================================

        if (step % moment_in_frame == 0) {

            phi.value_assert(10.0);

            //======= рассчет энергии слоя =======

            double energy = 0.0;
            for (int i = 0; i < N; ++i){
                double t = tmp[i]*phi[i]/2 - phi[i]*phi[i]*phi[i]*phi[i]/4;
                energy += t * h;
            }
            energies[step/moment_in_frame] = energy;

            //====================================

            phi.show("bin_data/data_" + std::to_string(step) + ".bin");

            std::cout << 100 * step / T_steps << '%' << std::endl;
        }

        phi.swap(phi_next);
    }

    energies.show(std::string("bin_data/energies.bin"));

    return 0;
}
