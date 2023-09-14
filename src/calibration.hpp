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

template <typename T>
inline JonesMatrix<T> calibrate_visibility(const JonesMatrix<T>& vis, const JonesMatrix<T>& solA, const JonesMatrix<T>& solB){
    JonesMatrix<T> tmp = vis * solB.conjtrans();
    return solA * tmp;
}


void read_solution(std::string filename, solution_t& sol);

template <typename T>
void apply_solution(Visibilities<T> &vis, solution_t& sol, unsigned int coarse_channel_index){
    if(sol.header.channelCount % 24 != 0) throw std::exception();
    if(sol.header.antennaCount != vis.obsInfo.nAntennas) throw std::exception();
    
    const unsigned int n_solfreq_per_band = sol.header.channelCount / 24u;
    const unsigned int solutions_offset = coarse_channel_index * n_solfreq_per_band;
    
    int channelRatio;
    if(n_solfreq_per_band > vis.nFrequencies)
        channelRatio = n_solfreq_per_band / vis.nFrequencies;
    else
        channelRatio =  vis.nFrequencies / n_solfreq_per_band;

    const unsigned int n_baselines {(vis.obsInfo.nAntennas + 1) * (vis.obsInfo.nAntennas / 2)};
    const size_t matrix_size {n_baselines * vis.obsInfo.nPolarizations * vis.obsInfo.nPolarizations * 2};
    const size_t n_interval_values {matrix_size * vis.nFrequencies};

    #pragma omp parallel for schedule(static) if(vis.integration_intervals() <= 10)
    for(unsigned int nInterval = 0; nInterval < vis.integration_intervals(); nInterval++){
       #pragma omp parallel for schedule(static) if(vis.integration_intervals() <= 10)
        for(unsigned int baseline = 0; baseline < n_baselines; baseline++){
            for(unsigned int ch {0}; ch < vis.nFrequencies; ch++){
                unsigned int solChannel;
                if(n_solfreq_per_band > vis.nFrequencies)
                    solChannel = (ch + solutions_offset) * channelRatio;
                else
                    solChannel = (ch + solutions_offset) / channelRatio;
            
                unsigned int a1 {static_cast<unsigned int>(-0.5 + std::sqrt(0.25 + 2*baseline))};
                unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
                JonesMatrix<double> &solA = reinterpret_cast<JonesMatrix<double>&>(sol.data[a1 * sol.header.channelCount + solChannel]);
                JonesMatrix<double> &solB = reinterpret_cast<JonesMatrix<double>&>(sol.data[a2 * sol.header.channelCount + solChannel]);
                T *data {reinterpret_cast<T*>(vis.data) + nInterval * n_interval_values + matrix_size * ch + 8 * baseline};
                
                JonesMatrix<double> visData {JonesMatrix<double>::from_array<double, T>(data)};
                JonesMatrix<double> tmp = calibrate_visibility(visData, solA, solB);

                data[0] = static_cast<T>(tmp.XX.real());
                data[1] = static_cast<T>(tmp.XX.imag());
                data[2] = static_cast<T>(tmp.XY.real());
                data[3] = static_cast<T>(tmp.XY.imag());
                data[4] = static_cast<T>(tmp.YX.real());
                data[5] = static_cast<T>(tmp.YX.imag());
                data[6] = static_cast<T>(tmp.YY.real());
                data[7] = static_cast<T>(tmp.YY.imag());
            }
        }
    }
}


#endif
