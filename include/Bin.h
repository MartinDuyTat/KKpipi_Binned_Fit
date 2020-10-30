// Martin Duy Tat 30th October 2020
/**
 * Bin is a class for a bin in phase space
 */
#ifndef BIN
#define BIN

class Bin {
  public:
    /**
     * Default constructor
     */
    Bin();
    /**
     * Function for adding pennies
     * @param pennies Number of pennies
     */
    void AddPennies(int &pennies);
    /**
     * Function for adding pounds
     * @param pounds Number of Pounds
     */
    void AddPounds(int &pounds);
    /** 
     * Function for getting the total amount of money
     * @return The total amount of money
     */
    double GetMoney();
  private:
    /**
     * Number of pennies
     */
    int m_pennies;
    /**
     * Number of pounds
     */
    int m_pounds;
};

#endif
