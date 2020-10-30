// Martin Duy Tat 30th October 2020

#include<Bin.h>

Bin::Bin(): m_pennies(0), m_pounds(0) {
}

void Bin::AddPennies(int &pennies) {
  m_pennies += pennies;
}

void Bin::AddPounds(int &pounds) {
  m_pounds += pounds;
}

double Bin::GetMoney() {
  return (double)m_pounds + (double)m_pennies/100;
}
