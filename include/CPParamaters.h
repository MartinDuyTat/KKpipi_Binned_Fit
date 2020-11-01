// Martin Duy Tat 31st October 2020
/**
 * CPParameters is a class that contains the CP parameters x and y
 */

#ifndef CPPARAMETERS
#define CPPARAMETERS

class {
  public:
    /**
     * Constructor that takes the CP parameters, x and y
     * @param xplus r_Bcos(delta_B + gamma) for Bplus decays
     * @param xminus r_Bcos(delta_B - gamma) for Bminus decays
     * @param yplus r_Bsin(delta_B + gamma) for Bplus decays
     * @param yminus r_Bsin(delta - gamma) for Bminus decays
     */
    CPParameters(double xplus, double xminus, double yplus, double yminus);
    /**
     * Function for getting CP parameters
     * @param xplus r_Bcos(delta_B + gamma) for Bplus decays
     * @param xminus r_Bcos(delta_B - gamma) for Bminus decays
     * @param yplus r_Bsin(delta_B + gamma) for Bplus decays
     * @param yminus r_Bsin(delta - gamma) for Bminus decays
     */
    void GetCPParameters(double &xplus, double &xminus, double &yplus, double& yminus);
  private:
    /**
     * xplus r_Bcos(delta_B + gamma) for Bplus decays
     */
    double m_xplus;
    /**
     * xminus r_Bcos(delta_B - gamma) for Bminus decays
     */
    double m_xminus;
    /**
     * yplus r_Bsin(delta_B + gamma) for Bplus decays
     */
    double m_yplus;
    /**
     * yminus r_Bcos(delta_B - gamma) for Bminus decays
     */
    double m_yminus;
};

#endif
