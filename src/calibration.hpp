#ifndef __CALIBRATION_H__
#define __CALIBRATION_H__

//#include <cctype>
#include <jones_matrix.hpp>
#include <memory_buffer.hpp>
#include <astroio.hpp>



class CalibrationSolutions : public MemoryBuffer<JonesMatrix<double>> {
    public:
    class Header {
        public:
        char intro[8];
        uint32_t file_type;
        uint32_t structure_type;
        uint32_t interval_count, antenna_count, channel_count, polarization_count;
        double start_time, end_time;
    };

    CalibrationSolutions::Header header;

    CalibrationSolutions(MemoryBuffer<JonesMatrix<double>>&& data, CalibrationSolutions::Header& header) : MemoryBuffer {data}, header {header} {};
    static CalibrationSolutions from_file(const std::string& filename);
};


template <typename T>
inline JonesMatrix<T> calibrate_visibility(const JonesMatrix<T>& vis, const JonesMatrix<T>& solA, const JonesMatrix<T>& solB){
    JonesMatrix<T> tmp = vis * solB.conjtrans();
    return solA * tmp;
}

void apply_solutions(Visibilities &vis, CalibrationSolutions& sol, unsigned int coarse_channel_index);



void apply_solutions_cpu(Visibilities &vis, const CalibrationSolutions& sol, unsigned int coarse_channel_index);

#endif
