//tests for rjMcMC1D classes
#include "../src/rjmcmc1d.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "logger.h"

class cLogger glog;

using ::testing::ElementsAre;
using ::testing::DoubleEq;

//dummy forward
class FakeSampler : public rjMcMC1DSampler {
  public:
    FakeSampler(std::vector<double> forwardret, std::vector<double> fakedata) {
      fret = forwardret;
      obs = fakedata;
      ndata = fakedata.size();
    }

    FakeSampler() {}

    std::vector<double> forwardmodel(const rjMcMC1DModel& m) {
      return fret;
    }
  private:
    std::vector<double> fret;
};

//tests the calculation of residuals given observed data
//and forward output.
TEST(rjMcMC1DSamplerTest, correctResiduals) {
  FakeSampler s = FakeSampler({1.0,2.0,3.0}, {2.0,4.0,5.0});
  rjMcMC1DModel m = rjMcMC1DModel();

  std::vector<double> output = s.computeresiduals(m);

  ASSERT_THAT(output, ElementsAre(DoubleEq(0.25),DoubleEq(0.25),DoubleEq(0.16)));
}

//check Gaussian works correctly
TEST(rjMcMC1DSamplerTest, testGaussianPDF) {
  FakeSampler s;
  double val = s.gaussian_pdf(1.0, 2.0, 3.0);
  ASSERT_DOUBLE_EQ(val, 1.0/(2.0*std::sqrt(2.0*M_PI))*exp(-1.0/2.0));
}
