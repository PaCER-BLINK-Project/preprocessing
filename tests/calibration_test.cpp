#include <exception>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <cstdio>

#include "../src/calibration.hpp"
#include "common.hpp"


/**
 * @brief Compute the conjugate cross multiply (ccm) between complex numbers `a` and `b`
 * and add the result to the `res` accumulator.
 * 
 * @tparam T1 
 * @tparam T2 
 * @param a 
 * @param b 
 * @param res 
 */
template <typename T1, typename T2>
void ccm(const std::complex<T1>& a, const std::complex<T1>& b, std::complex<T2>& res){
    res.real(res.real() + static_cast<T2>(a.real()) * b.real() + static_cast<T2>(a.imag()) * b.imag());
    res.imag(res.imag() + static_cast<T2>(a.imag()) * b.real() - static_cast<T2>(a.real()) * b.imag());
}


// Stores the path to the directory containing the test data.
std::string dataRootDir;


template <typename T>
void read_fits_file(std::string filename, T*& buffer, size_t& length){
    FITS my_fits {FITS::from_file(filename)};
    auto hdu = my_fits[0];
    length = hdu.get_xdim() * hdu.get_ydim();
    T *tmp = new T[length];
    memcpy(tmp, hdu.get_image_data(), length * sizeof(T));
    buffer = tmp;
}

/*
vis_from_files

Reads visibility matrix from a set of files as produced by Marcin's code to "dump a CASA measurement set.".
*/
template <typename T>
Visibilities<T> vis_from_files(const std::string basepath, const std::string obs_id){
    T *buffer (nullptr), *data {nullptr};
    size_t length {0};
    const int nAntennas {128}, nPolarisations {4};
    const int nBaselines {(nAntennas / 2) * (nAntennas + 1)};
    for(int cmplx {0}; cmplx < 2; cmplx++){
        for(int pol {0}; pol < nPolarisations; pol++){
            std::stringstream filename;
            filename << basepath << "/" << obs_id << "_vis_" << (cmplx ? "imag" : "real") << "_channel000_time000000_pol" << pol << ".fits";
            read_fits_file(filename.str(), buffer, length);
            if(!data){
                data = new T[nBaselines * nPolarisations * 2];
            }
            // Data is stored in upper triangular form, but we want it in lower triangular form
            // as an array of Jones matrices indexed by baseline.
            for(int a1 {0}; a1 < nAntennas; a1++){
                for(int a2 {a1}; a2 < nAntennas; a2++){
                    int baseline {(a2 * (a2 + 1)) / 2 + a1}; 
                    data[baseline * nPolarisations * 2 + pol * 2 + cmplx] = buffer[a1 * nAntennas + a2];
                }
            }
            delete[] buffer;
            buffer = nullptr;
        }
    }
    // Each Jones matrix has to be conj-transposed because we went from upper to lower triangular
    // form. 
    for(int c {0}; c < nBaselines; c++){
        JonesMatrix<T>* m = reinterpret_cast<JonesMatrix<T>*>(&data[c * 8]);
        *m = m->conjtrans();        
    }
    return Visibilities<T>{reinterpret_cast<std::complex<T>*>(data), VCS_OBSERVATION_INFO, VCS_OBSERVATION_INFO.nTimesteps, VCS_OBSERVATION_INFO.nFrequencies};
}



void test_read_solution(){
    // as read by aocal.py
    double rawTestData[] {-0.7241842844376969, 1.1922007333545488, 0.0, 0.0, 0.0, 0.0, -0.14706941378717175, 0.9448835886397396};
    JonesMatrix<double>& jTestData {reinterpret_cast<JonesMatrix<double>&>(rawTestData)};
    
    solution_t sol;
    std::string solpath {dataRootDir + "/mwa/1103645160/CASA/calibrated/1103645160.bin"};
    read_solution(solpath, sol);
    if(sol.data[0] != jTestData) throw TestFailed("test_read_solution failed.");
    std::cout << "'test_read_solution' passed." << std::endl;
}


/*
Tests the calibration of a single Jones matrix.
*/
void test_apply_one_solution(){
    // antenna 0
    double _solA[] {-0.72418428, 1.19220073, 0., 0, 0., 0., -0.14706941, 0.94488359};
    JonesMatrix<double>& solA {reinterpret_cast<JonesMatrix<double>&>(_solA)};
    // antenna 1
    double _solB[] {0.57376861, 1.29451097, 0., 0., 0., 0., 0.24960554, 0.88270966};
    JonesMatrix<double>& solB {reinterpret_cast<JonesMatrix<double>&>(_solB)};

    double _uncalibrated[] {-10.99070549, -62.03197861, 90.57805634, -176.05194092,
                            -236.68218994, -140.27902222, 549.13391113, -132.24552917};
    JonesMatrix<double>& uncalibrated {reinterpret_cast<JonesMatrix<double>&>(_uncalibrated)};

    double _expected_calibrated[] {88.19023895, -87.78139496, 243.87818909, -68.59230804,
                            -166.77030945, -333.1229248, 486.20907593, +95.35482025};
    JonesMatrix<double>& cal {reinterpret_cast<JonesMatrix<double>&>(_expected_calibrated)};

    auto comp_sol = calibrate_visibility(uncalibrated, solA, solB);
    auto diff = (comp_sol - cal).max_abs();
    if(diff >= 1e-4){
        throw TestFailed("'test_apply_one_solution' failed.");
    }
    std::cout << "'test_apply_one_solution' passed." << std::endl;
}


/**
 * @brief Test the calibration process on the test dataset produced by Marcin's program. The dataset contains
   a dump of a CASA measurement set, both the calibrated visibilities and the non-calibrated ones.
*/
void test_calibration(){
    using vistype = double;
    Visibilities<vistype> uncalibrated = vis_from_files<vistype>(dataRootDir + "/mwa/1103645160/CASA/uncalibrated", "1103645160");
    Visibilities<vistype> reference_calibrated = vis_from_files<vistype>(dataRootDir + "/mwa/1103645160/CASA/calibrated", "1103645160");

    solution_t sol;
    std::string solpath {dataRootDir + "/mwa/1103645160/CASA/calibrated/1103645160.bin"};
    read_solution(solpath, sol);
    apply_solution<vistype>(uncalibrated, sol, 0);
    double maxdiff {0.0};
    const unsigned int nBaselines {uncalibrated.obsInfo.nAntennas / 2 * (uncalibrated.obsInfo.nAntennas + 1)};
    const unsigned int nPolarisations {4};
    for(unsigned int baseline {0}; baseline < nBaselines; baseline++){
        unsigned int a1 {static_cast<unsigned int>(-0.5 + std::sqrt(0.25 + 2*baseline))};
        unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
        if(a1 == a2) continue;
        JonesMatrix<vistype>* m1 = reinterpret_cast<JonesMatrix<vistype>*>(&((vistype*)uncalibrated.data)[baseline * nPolarisations * 2]);
        JonesMatrix<vistype>* m2 = reinterpret_cast<JonesMatrix<vistype>*>(&((vistype*)reference_calibrated.data)[baseline * nPolarisations * 2]);
        
        double diff {(*m1 - *m2).max_abs()};
        if(diff > maxdiff) maxdiff = diff;

        if(diff > 1e-3){
            std::cout << "Error for baseline = " << baseline << " (antennas " << a1 << ", " << a2 <<"): \n m1 = " << *m1 << ", m2 = " << *m2 << std::endl; 
            throw TestFailed("test_calibration failed.");
        }
    }
    std::cout << "'test_calibration' passed. maxdiff = " << maxdiff << "." << std::endl;
}



/**
 * @brief This test confirms that the calibrated visibilities 
 * of the lower triangular matrix are equal to the calibrated 
 * conjugate-transposed visibilities of the upper triangular matrix.
 * This is useful to use the "CASA dump" as test case since those
 * are stored in upper triangular form, whereas our code uses the
 * lower-triangular form representation.
*/
void test_reciprocal(){

    std::complex<double> X1 {53, 41};
    std::complex<double> Y1 {33, 71};
    std::complex<double> X2 {35, 11};
    std::complex<double> Y2 {13, 76};

     // antenna 0
    double _solA[] {-0.72418428, 1.19220073, 0., 0, 0., 0., -0.14706941, 0.94488359};
    JonesMatrix<double>& solA {reinterpret_cast<JonesMatrix<double>&>(_solA)};
    // antenna 1
    double _solB[] {0.57376861, 1.29451097, 0., 0., 0., 0., 0.24960554, 0.88270966};
    JonesMatrix<double>& solB {reinterpret_cast<JonesMatrix<double>&>(_solB)};


    
    JonesMatrix<double> corr1to2 {};
    ccm(X1, X2, corr1to2.XX);
    ccm(X1, Y2, corr1to2.XY);
    ccm(Y1, X2, corr1to2.YX);
    ccm(Y1, Y2, corr1to2.YY);
    
    auto corr1to2app = calibrate_visibility(corr1to2, solA, solB);
    
    JonesMatrix<double> corr2to1 {};
    ccm(X2, X1, corr2to1.XX);
    ccm(X2, Y1, corr2to1.XY);
    ccm(Y2, X1, corr2to1.YX);
    ccm(Y2, Y1, corr2to1.YY);
    auto corr2to1app = calibrate_visibility(corr2to1, solB, solA);
    auto corr2to1appTransposed = corr2to1app.conjtrans();
    auto diff = (corr2to1appTransposed - corr1to2app).max_abs();
    if(diff > 1e-5){
        std::cout << "corr2to1appTransposed is " << corr2to1appTransposed << ", corr1to2app = " << corr1to2app << std::endl;
        throw TestFailed("'test_reciproc' failed.");
    }
    std::cout << "'test_reciproc' passed (diff = " << diff << ")." << std::endl;
}


int main(void){
    
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    dataRootDir = std::string {pathToData};

    try{
        test_read_solution();
        test_apply_one_solution();
        test_calibration();
        test_reciprocal();
     
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
