#ifndef __CALIBRATION_H__
#define __CALIBRATION_H__

#include <cctype>
#include <complex>
#include <exception>
#include <astroio.hpp>


struct header_t {
    char intro[8];
    uint32_t fileType;
    uint32_t structureType;
    uint32_t intervalCount, antennaCount, channelCount, polarizationCount;
    double startTime, endTime;
};


template
<typename T>
class JonesMatrix {

    public:

    std::complex<T> XX;
    std::complex<T> XY;
    std::complex<T> YX;
    std::complex<T> YY;

    inline friend std::istream& operator>>(std::istream& is, JonesMatrix& m){
        is.read(reinterpret_cast<char*>(&m.XX), sizeof(std::complex<T>));
        is.read(reinterpret_cast<char*>(&m.XY), sizeof(std::complex<T>));
        is.read(reinterpret_cast<char*>(&m.YX), sizeof(std::complex<T>));
        is.read(reinterpret_cast<char*>(&m.YY), sizeof(std::complex<T>));
        return is;
    }


    template <typename T1, typename T2>
    static JonesMatrix<T1> from_array(T2 *data){
        JonesMatrix<T1> res;
        res.XX = {static_cast<T1>(data[0]),  static_cast<T1>(data[1])};
        res.XY = {static_cast<T1>(data[2]),  static_cast<T1>(data[3])};
        res.YX = {static_cast<T1>(data[4]),  static_cast<T1>(data[5])};
        res.YY = {static_cast<T1>(data[6]),  static_cast<T1>(data[7])};
        return res;
    }

    bool operator==(const JonesMatrix& other){
        return XX == other.XX && XY == other.XY && YX == other.YX && YY == other.YY;
    }

    bool operator!=(const JonesMatrix& other){
        return !(*this == other);
    }


    JonesMatrix operator*(const JonesMatrix& other) const {
        JonesMatrix res;
        res.XX = XX * other.XX + XY * other.YX;
        res.XY = XX * other.XY + XY * other.YY;
        res.YX = YX * other.XX + YY * other.YX;
        res.YY = YX * other.XY + YY * other.YY;
        return res;
    }



    JonesMatrix operator-(const JonesMatrix& other) const {
        JonesMatrix res;
        res.XX = XX - other.XX;
        res.XY = XY - other.XY;
        res.YX = YX - other.YX;
        res.YY = YY - other.YY;
        return res;
    }

    double max_abs() const {
        double max {std::abs(XX)};
        if(std::abs(XY) > max) max = std::abs(XY);
        if(std::abs(YX) > max) max = std::abs(XY);
        if(std::abs(YY) > max) max = std::abs(XY);
        return max;
    }

    JonesMatrix conjtrans() const {
        JonesMatrix res;
        res.XX = std::conj(XX);
        res.XY = std::conj(YX);
        res.YY = std::conj(YY);
        res.YX = std::conj(XY);
        return res;
    }

    friend std::ostream& operator<<(std::ostream& of, const JonesMatrix& m){
        std::cout << "[" << m.XX << ", " << m.XY  << ", " << m.YX << ", " << m.YY << "]" << std::endl;
        return of;
    }
};


struct solution_t {
    header_t header;
    JonesMatrix<double> *data;
};



void read_solution(std::string filename, solution_t& sol);



template <typename T>
inline JonesMatrix<T> calibrate_visibility(const JonesMatrix<T>& vis, const JonesMatrix<T>& solA, const JonesMatrix<T>& solB){
    JonesMatrix<T> tmp = vis * solB.conjtrans();
    return solA * tmp;
}



void apply_solution(Visibilities &vis, solution_t& sol, unsigned int coarse_channel_index);

#endif
