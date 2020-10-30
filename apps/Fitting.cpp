// Martin Duy Tat 30th October 2020
/**
 * Fitting is the program for doing the binned fitting
 */

#include<Bin.h>
#include<iostream>

int main() {
  Bin money;
  money.AddPennies(35);
  money.AddPounds(5);
  double amount = money.GetMoney();
  std::cout << amount << std::endl;
  return 0;
}
