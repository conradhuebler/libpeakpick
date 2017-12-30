#include <libpeakpick/deconvulate.h>
#include <libpeakpick/spectrum.h>


typedef Eigen::VectorXd Vector;

int FirstLorentzian()
{
    PeakPick::spectrum spec = loadFromFile("../samples/lorentzian_1", 0 , 10);

    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.PosMax()  << "," << spec.Max() << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.X(spec.XtoIndex(spec.PosMax() )) << "," << spec.Y(spec.PosMax()) << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    Vector vector(1);
    vector(0) =  5 ;
    PeakPick::GLFit fit(&spec, 0, spec.size());
    fit.setGuess(vector);
    PeakPick::FitResult *result = fit.Deconvulate();
    //FitResult *result = PeakPick::LiberalDeconvulate(&spec, 0, spec.size(), 1, vector, 0);
    std::cout << result->parameter << std::endl;
    std::cout << "Integral " << result->integral << std::endl;
    //spec.print();

    
    
    return 1;
}

int SecondLorentzian()
{
    PeakPick::spectrum spec = loadFromFile("../samples/lorentzian_2", 0 , 10);

    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.PosMax()  << "," << spec.Max() << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    // spec.print();

    
    
    return 1;
}

int testLorentzien()
{
    FirstLorentzian();
    SecondLorentzian();

    return 1;
}
