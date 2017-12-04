int FirstLorentzian()
{
    PeakPick::spectrum spec = loadFromFile("../samples/lorentzian_1", 0 , 10);

    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.PosMax() + 1 << "," << spec.Max() << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    
    
    
    return 1;
}

int SecondLorentzian()
{
    PeakPick::spectrum spec = loadFromFile("../samples/lorentzian_2", 0 , 10);

    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.PosMax() + 1 << "," << spec.Max() << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    
    
    
    return 1;
}

int testLorentzien()
{
    FirstLorentzian();
    SecondLorentzian();

    return 1;
}
