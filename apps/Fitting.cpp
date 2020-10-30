// Martin Duy Tat 30th October 2020
/**
 * Fitting is the program for doing the binned fitting
 */

#include<Bin.h>
#include<iostream>

int main() {
  Bin money;
  int a = 35, b = 5;
  money.AddPennies(a);
  money.AddPounds(b);
  double amount = money.GetMoney();
  std::cout << amount << std::endl;
  return 0;
}
