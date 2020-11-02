// Martin Duy Tat 2nd November

#include<vector>
#include<string>
#include<complex>
#include"Amplitude.h"

Amplitude::Amplitude(const std::string &Damplitude, const std::string &DBARamplitude) {
  void *my_lib_handle_d = dlopen(Damplitude.c_str(), RTLD_NOW);
  void *my_lib_handle_dbar = dlopen(Damplitude.c_str(), RTLD_NOW);
  m_Damplitude = (std::complex<double> (*)(double const *event, const int &x1)) dlsym(my_lib_handle_d, "AMP");
  m_DBARamplitude_dbar = (std::complex<double> (*)(double const *event, const int &x1)) dlsym(my_lib_handle_dbar, "AMP");
}

std::complex<double> Amplitude::operator()(const std::vector<double> &event, int conj) {
  if(conj == +1) {
    return m_Damlitude(event, +1);
  } else if(conj == -1) {
    return m_DBARamplitude(event, +1);
  }
}
