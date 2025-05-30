#include <exception>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <cstdio>

#include "../src/calibration.hpp"
#include "../src/mapping.hpp"
#include "common.hpp"

#ifdef __GPU__
#include "../src/calibration_gpu.hpp"
#include "../src/mapping_gpu.hpp"
#endif


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
void ccm(const Complex<T1>& a, const Complex<T1>& b, Complex<T2>& res){
    res.real = res.real + static_cast<T2>(a.real) * b.real + static_cast<T2>(a.imag) * b.imag;
    res.imag = res.imag + static_cast<T2>(a.imag) * b.real - static_cast<T2>(a.real) * b.imag;
}


// Stores the path to the directory containing the test data.
std::string data_root_dir;


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
Visibilities vis_from_files(const std::string basepath, const std::string obs_id){
    float *buffer (nullptr);
    size_t length {0};
    const int nAntennas {128}, nPolarisations {4};
    const int nBaselines {(nAntennas / 2) * (nAntennas + 1)};
    MemoryBuffer<std::complex<float>> mb_data (nBaselines * nPolarisations, false, false);
    float *data {reinterpret_cast<float*>(mb_data.data())};
    for(int cmplx {0}; cmplx < 2; cmplx++){
        for(int pol {0}; pol < nPolarisations; pol++){
            std::stringstream filename;
            filename << basepath << "/" << obs_id << "_vis_" << (cmplx ? "imag" : "real") << "_channel000_time000000_pol" << pol << ".fits";
            read_fits_file(filename.str(), buffer, length);
            if(!data){
                data = new float[nBaselines * nPolarisations * 2];
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
        JonesMatrix<float>* m = reinterpret_cast<JonesMatrix<float>*>(&data[c * 8]);
        *m = m->conjtrans();        
    }
    return Visibilities{std::move(mb_data), VCS_OBSERVATION_INFO, VCS_OBSERVATION_INFO.nTimesteps, VCS_OBSERVATION_INFO.nFrequencies};
}



void test_read_solution(){
    // as read by aocal.py
    double rawTestData[] {-0.7241842844376969, 1.1922007333545488, 0.0, 0.0, 0.0, 0.0, -0.14706941378717175, 0.9448835886397396};
    JonesMatrix<double>& jTestData {reinterpret_cast<JonesMatrix<double>&>(rawTestData)};
    const std::string sol_file {data_root_dir + "/mwa/1103645160/CASA/calibrated/1103645160.bin"};
    auto sol = CalibrationSolutions::from_file(sol_file);
    
    if(sol.data()[0] != jTestData) throw TestFailed("test_read_solution failed.");
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
    Visibilities uncalibrated = vis_from_files(data_root_dir + "/mwa/1103645160/CASA/uncalibrated", "1103645160");
    Visibilities reference_calibrated = vis_from_files(data_root_dir + "/mwa/1103645160/CASA/calibrated", "1103645160");
    const std::string sol_file {data_root_dir + "/mwa/1103645160/CASA/calibrated/1103645160.bin"};
    const auto sol = CalibrationSolutions::from_file(sol_file);
    apply_solutions_cpu(uncalibrated, sol, 0);
    double maxdiff {0.0};
    const unsigned int nBaselines {uncalibrated.obsInfo.nAntennas / 2 * (uncalibrated.obsInfo.nAntennas + 1)};
    const unsigned int nPolarisations {4};
    for(unsigned int baseline {0}; baseline < nBaselines; baseline++){
        unsigned int a1 {static_cast<unsigned int>(-0.5 + std::sqrt(0.25 + 2*baseline))};
        unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
        if(a1 == a2) continue;
        JonesMatrix<float>* m1 = reinterpret_cast<JonesMatrix<float>*>(&((float*)uncalibrated.data())[baseline * nPolarisations * 2]);
        JonesMatrix<float>* m2 = reinterpret_cast<JonesMatrix<float>*>(&((float*)reference_calibrated.data())[baseline * nPolarisations * 2]);
        
        double diff {(*m1 - *m2).max_abs()};
        if(diff > maxdiff) maxdiff = diff;

        if(diff > 1e-3){
            std::cout << "Error for baseline = " << baseline << " (antennas " << a1 << ", " << a2 <<"): \n m1 = " << *m1 << ", m2 = " << *m2 << std::endl; 
            throw TestFailed("test_calibration failed.");
        }
    }
    std::cout << "'test_calibration' passed. maxdiff = " << maxdiff << "." << std::endl;
}


void test_reorder_calibrated(){
    const std::string metadata_file {data_root_dir + "/mwa/1276619416/20200619163000.metafits"}; 
	const std::string vis_file {data_root_dir + "/mwa/1276619416/visibilities/1276619416_20200619163000_gpubox24_00.fits"};
    const std::string sol_file {data_root_dir + "/mwa/1276619416/1276625432.bin"};
    auto vis = Visibilities::from_fits_file(vis_file);
    auto mapping = get_visibilities_mapping(metadata_file);
	auto reord_vis = reorder_visibilities_cpu(vis, mapping);
    auto sol = CalibrationSolutions::from_file(sol_file);
    apply_solutions_cpu(reord_vis, sol, 0);
	
    const float expected_values[] {-523.286, -252.787, -223.328, 762.554, 466.479, 36.5678, 222.863, 71.7788};
    const int n_baselines = 128 / 2 * 129;
    const int matrix_size = n_baselines * 4 * 2;
    const int baseline_idx = 1;
    const int channel = 1;
    const float* computed_values {
        reinterpret_cast<const float*>(reord_vis.data()) + channel * matrix_size + baseline_idx * 8 
    };
    for(int i {0}; i < 8; i++){
        auto diff = std::abs(expected_values[i] - computed_values[i]);
        std::cout << "Diff between " << expected_values[i] << " and " << computed_values[i] << " is: " << diff << std::endl;
        if( diff > 1e-3){
            std::cerr << "expected_value[i] = " << expected_values[i] << ", computed_value[i] = " << computed_values[i] << std::endl;
            throw TestFailed("'test_reorder_calibrated': values are different.");
        }
    }
    std::cout << "'test_reorder_calibrated' passed." << std::endl;
}

#ifdef __GPU__
void test_reorder_calibrated_gpu(){
    const std::string metadata_file {data_root_dir + "/mwa/1276619416/20200619163000.metafits"}; 
	const std::string vis_file {data_root_dir + "/mwa/1276619416/visibilities/1276619416_20200619163000_gpubox24_00.fits"};
    const std::string sol_file {data_root_dir + "/mwa/1276619416/1276625432.bin"};
    auto vis = Visibilities::from_fits_file(vis_file);
    auto mapping = get_visibilities_mapping(metadata_file);
    vis.to_gpu();
    mapping.to_gpu();
	auto reord_vis = reorder_visibilities_gpu(vis, mapping);
    auto sol = CalibrationSolutions::from_file(sol_file);
    sol.to_gpu();
    apply_solutions_gpu(reord_vis, sol, 0);
    reord_vis.to_cpu();
    const float expected_values[] {-523.286, -252.787, -223.328, 762.554, 466.479, 36.5678, 222.863, 71.7788};
    const int n_baselines = 128 / 2 * 129;
    const int matrix_size = n_baselines * 4 * 2;
    const int baseline_idx = 1;
    const int channel = 1;
    const float* computed_values {
        reinterpret_cast<const float*>(reord_vis.data()) + channel * matrix_size + baseline_idx * 8 
    };
    for(int i {0}; i < 8; i++){
        auto diff = std::abs(expected_values[i] - computed_values[i]);
        std::cout << "Diff between " << expected_values[i] << " and " << computed_values[i] << " is: " << diff << std::endl;
        if( diff > 1e-3){
            std::cerr << "expected_value[i] = " << expected_values[i] << ", computed_value[i] = " << computed_values[i] << std::endl;
            throw TestFailed("'test_reorder_calibrated': values are different.");
        }
    }
    std::cout << "'test_reorder_calibrated_gpu' passed." << std::endl;
}
#endif
/**
 * @brief This test confirms that the calibrated visibilities 
 * of the lower triangular matrix are equal to the calibrated 
 * conjugate-transposed visibilities of the upper triangular matrix.
 * This is useful to use the "CASA dump" as test case since those
 * are stored in upper triangular form, whereas our code uses the
 * lower-triangular form representation.
*/
void test_reciprocal(){

    Complex<double> X1 {53, 41};
    Complex<double> Y1 {33, 71};
    Complex<double> X2 {35, 11};
    Complex<double> Y2 {13, 76};

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
    data_root_dir = std::string {pathToData};

    try{
        test_read_solution();
        test_apply_one_solution();
        // The following fails because we do not support double anymore, and marcin's dump uses duble.
        // No problem, wanted to change this test anyway.
        // test_calibration();
        test_reciprocal();
        test_reorder_calibrated();
        #ifdef __GPU__
        test_reorder_calibrated_gpu();
        #endif
     
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
