/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Authors:
Richard L. Taylor, Geoscience Australia.
*/

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

//dummy forward,
//just returns the same output regardless of model
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
TEST(rjMcMC1DSamplerTest, test_residuals) {
  FakeSampler s = FakeSampler({1.0,2.0,3.0}, {2.0,4.0,5.0});
  rjMcMC1DModel m = rjMcMC1DModel();

  std::vector<double> output = s.computeresiduals(m);

  ASSERT_THAT(output, ElementsAre(DoubleEq(0.25),DoubleEq(0.25),DoubleEq(0.16)));
}

//check Gaussian works correctly
TEST(rjMcMC1DSamplerTest, test_gaussian_pdf_value) {
  FakeSampler s;
  double val = s.gaussian_pdf(1.0, 2.0, 3.0);
  ASSERT_DOUBLE_EQ(val, 1.0/(2.0*std::sqrt(2.0*M_PI))*exp(-1.0/2.0));
}

//check that the misfit is calculated correctly
TEST(rjMcMC1DSamplerTest,test_misfit) {
  FakeSampler s = FakeSampler({2.0,3.0,4.0}, {4.0,4.0,5.0});
  rjMcMC1DModel m = rjMcMC1DModel();
  m.nvar = {0.5, 1.0, 2.0};

  s.set_misfit(m);

  ASSERT_DOUBLE_EQ(m.get_misfit(),
  0.5+std::log(0.5)+
  0.0625+std::log(1.0)+
  0.02+std::log(2.0));

  ASSERT_DOUBLE_EQ(s.get_normalised_misfit(m), m.get_misfit()/3.0);

}

//layer operators

class rjMcMC1DLayerTest : public ::testing::Test {
protected:
  void SetUp() override {
    a.ptop = 10.0;
    b.ptop = 20.0;

    a.value = 0.0;
    b.value = 1.0;
  }

  rjMcMC1DLayer a;
  rjMcMC1DLayer b;
};

TEST_F(rjMcMC1DLayerTest, test_equals) {
  EXPECT_TRUE(a == a);
  EXPECT_TRUE(b == b);
  EXPECT_FALSE(a == b);
  b.ptop = 10.0;
  EXPECT_TRUE(a == b);
}

TEST_F(rjMcMC1DLayerTest, test_lt) {
  EXPECT_TRUE(a < b);
  EXPECT_FALSE(b < a);
}

TEST_F(rjMcMC1DLayerTest, test_gt) {
  EXPECT_FALSE(a > b);
  EXPECT_TRUE(b > a);
}

//acceptance ratio calculation
TEST(cProposalTest, test_ar) {
  cProposal cp(cProposal::Type::VALUECHANGE);

  for (size_t i = 0; i < 100; i++)   {
    cp.inc_np();
  }
  for (size_t i = 0; i < 49; i++) {
    cp.inc_na();
  }
  ASSERT_DOUBLE_EQ(cp.ar(),49.0);
}

//testing model
class rjMcMC1DModelTest : public ::testing::Test {
protected:
  void SetUp() override {
    m.initialise(depth_max, cond_min, cond_max);
    m.nvar = {1.0, 1.0, 1.0};
  }

  const double cond_max = 2.0;
  const double cond_min = -3.0;
  const double depth_max = 99.0;

  rjMcMC1DModel m;
};

TEST_F(rjMcMC1DModelTest, test_initial_vecs_empty) {
  EXPECT_TRUE(m.nlayers() == 0);
  EXPECT_TRUE(m.nnuisances() == 0);
  EXPECT_TRUE(m.nnoises() == 0);
}

TEST_F(rjMcMC1DModelTest, test_set_get_misfit) {
  m.set_misfit(3.0);
  EXPECT_TRUE(m.get_misfit()==3.0);
}

TEST_F(rjMcMC1DModelTest, test_set_get_residuals_chi2) {
  m.nvar = {1.0, 2.0, 3.0};
  m.set_residuals({1.0, 1.0, 1.5});
  EXPECT_THAT(m.get_residuals(), ElementsAre(DoubleEq(1.0),DoubleEq(1.0),DoubleEq(1.5)));
  EXPECT_DOUBLE_EQ(m.get_chi2(),2.0);
}

TEST_F(rjMcMC1DModelTest, test_insert_interface) {
  // test that interface insertion passes
  // and fails correctly if position/conductivity
  // is out of bounds specified in model init
  ASSERT_TRUE(m.insert_interface(0.0,cond_max - 1.0));
  ASSERT_TRUE(m.insert_interface(10.0,cond_min + 1.0));
  ASSERT_TRUE(m.insert_interface(30.0,cond_min));
  ASSERT_TRUE(m.insert_interface(20.0,cond_max));

  //zero thickness layer fails
  ASSERT_FALSE(m.insert_interface(10.0,0.0));
  //above max/below min cond
  ASSERT_FALSE(m.insert_interface(40.0,cond_max + 1.0));
  ASSERT_FALSE(m.insert_interface(40.0,cond_min - 1.0));

  //all passed layers should be in the model
  ASSERT_TRUE(m.nlayers() == 4);

  //sorted order after insertion
  std::vector<double> vals = m.getvalues();
  EXPECT_THAT(vals,ElementsAre(DoubleEq(cond_max - 1.0),
                               DoubleEq(cond_min + 1.0),
                               DoubleEq(cond_max),
                               DoubleEq(cond_min)));
  //thicknesses correct
  std::vector<double> thicc = m.getthicknesses();
  ASSERT_THAT(thicc,ElementsAre(DoubleEq(10.0),
                                DoubleEq(10.0),
                                DoubleEq(10.0)));

}

TEST_F(rjMcMC1DModelTest, test_delete_interface) {
  // can't delete interface from an empty model
  EXPECT_FALSE(m.delete_interface(0));
  m.insert_interface(0.0,cond_max - 1.0);
  m.insert_interface(20.0,cond_max - 1.5);
  m.insert_interface(10.0,cond_max - 2.0);

  ASSERT_TRUE(m.delete_interface(1));
  ASSERT_TRUE(m.nlayers() == 2);
  //did we delete the middle interface?
  std::vector<double> vals = m.getvalues();
  EXPECT_THAT(vals, ElementsAre(DoubleEq(cond_max - 1.0),
                                DoubleEq(cond_max - 1.5)));
  std::vector<double> thicc = m.getthicknesses();
  ASSERT_THAT(thicc, ElementsAre(DoubleEq(20.0)));
}

TEST_F(rjMcMC1DModelTest, test_which_layer) {
  m.insert_interface(0.0,cond_max - 1.0);
  m.insert_interface(10.0,cond_max - 1.5);
  m.insert_interface(20.0,cond_max - 2.0);

  EXPECT_TRUE(m.which_layer(0.0)==0);
  EXPECT_TRUE(m.which_layer(15.0)==1);
  EXPECT_TRUE(m.which_layer(25.0)==2);

  //test layer edge behaviour
  EXPECT_TRUE(m.which_layer(10.0)==1);
  EXPECT_TRUE(m.which_layer(20.0)==2);
}

TEST_F(rjMcMC1DModelTest, test_move_interface) {
  m.insert_interface(0.0,cond_max - 1.0);
  m.insert_interface(10.0,cond_max - 1.5);
  m.insert_interface(20.0,cond_max - 2.0);

  //can't make pos negative or move top layer
  ASSERT_FALSE(m.move_interface(1,-10.0));
  ASSERT_FALSE(m.move_interface(0,5.0));

  EXPECT_TRUE(m.move_interface(1,15.0));

  //move middle layer to bottom then sanity check
  ASSERT_TRUE(m.move_interface(1,25.0));
  std::vector<double> vals = m.getvalues();
  EXPECT_TRUE(vals[2] == cond_max-1.5);
  std::vector<double> thicc = m.getthicknesses();
  ASSERT_THAT(thicc, ElementsAre(DoubleEq(20.0),DoubleEq(5.0)));
}

TEST_F(rjMcMC1DModelTest, test_add_noise) {
  //not much is actually happening with noise computation
  //in the model objects, they basically just act as containers
  //so we just check that the counter returns correctly
  EXPECT_TRUE(m.nnoises()==0);

  rjMcMCNoise noise;
  noise.value = 1.0;
  noise.min = 0.5;
  noise.max = 1.5;
  noise.sd_valuechange = 0.25;

  m.mnoises.push_back(noise);

  ASSERT_TRUE(m.nnoises()==1);
}

//rjMcMC1DModel has a container for nuisance shared_ptrs.
//need a dummy implementation to instantiate things
class dummyNuisance : public rjMcMCNuisance {
public:
  void settype(std::string s) {
    return;
  }
  std::string typestring() const{
    return "dummy";
  }

  std::shared_ptr<rjMcMCNuisance> deepcopy() {
    std::shared_ptr<dummyNuisance> dup = std::make_shared<dummyNuisance>();
    dup->value = value;
    dup->min = min;
    dup->max = max;
    dup->sd_valuechange = sd_valuechange;

    return std::shared_ptr<rjMcMCNuisance>(dup);
  }
};

TEST_F(rjMcMC1DModelTest, test_add_nuisance) {
  EXPECT_TRUE(m.nnuisances()==0);

  dummyNuisance nu;
  nu.value = 30.0;
  nu.min = 10.0;
  nu.max = 33.0;
  nu.sd_valuechange = 5.0;

  m.nuisances.push_back(nu.deepcopy());

  ASSERT_TRUE(m.nnuisances()==1);
}

//test the noise and nuisance mapping functionality
class rjMcMC1DNoiseMapTest : public ::testing::Test {
protected:
  void SetUp() override {
    nmap.resettozero();
    rjMcMCNoise noise;
    noise.min = 0.01;
    noise.max = 0.10;
    noise.value = 0.05;
    noise.sd_valuechange = 0.03;
    noise.data_bounds = {0,2};
    m.mnoises.push_back(noise);
  }
  rjMcMC1DNoiseMap nmap;
  rjMcMC1DModel m;
};

TEST_F(rjMcMC1DNoiseMapTest, test_init) {
  EXPECT_TRUE(nmap.get_nnoises()==0);
  EXPECT_TRUE(nmap.get_nentries()==0);
}

TEST_F(rjMcMC1DNoiseMapTest, test_add_model) {
  nmap.addmodel(m);
  EXPECT_TRUE(nmap.get_nnoises()==1);
  EXPECT_TRUE(nmap.get_nentries()==1);

  EXPECT_THAT(nmap.noises[0], ElementsAre(DoubleEq(0.05)));

  m.mnoises[0].value = 0.09;
  nmap.addmodel(m);
  EXPECT_THAT(nmap.noises[0], ElementsAre(DoubleEq(0.05),DoubleEq(0.09)));
}

TEST_F(rjMcMC1DNoiseMapTest, test_reset) {
  nmap.addmodel(m);
  EXPECT_TRUE(nmap.get_nentries()==1);
  nmap.resettozero();
  EXPECT_TRUE(nmap.get_nentries()==0);
  EXPECT_TRUE(nmap.get_nnoises()==0);
  EXPECT_TRUE(nmap.noises.size()==0);

}

class rjMcMC1DNuisanceMapTest : public ::testing::Test {
protected:
  void SetUp() override {
    numap.resettozero();
    dummyNuisance nuisance;
    nuisance.min = 1.0;
    nuisance.max = 10.0;
    nuisance.value = 5.0;
    nuisance.sd_valuechange = 0.03;
    m.nuisances.push_back(nuisance.deepcopy());
  }
  rjMcMC1DNuisanceMap numap;
  rjMcMC1DModel m;
};

TEST_F(rjMcMC1DNuisanceMapTest, test_init) {
  EXPECT_TRUE(numap.get_nentries()==0);
}

TEST_F(rjMcMC1DNuisanceMapTest, test_add_model) {
  numap.addmodel(m);
  EXPECT_TRUE(numap.get_nentries()==1);

  EXPECT_THAT(numap.nuisance[0], ElementsAre(DoubleEq(5.0)));

  m.nuisances[0]->value = 9.0;
  numap.addmodel(m);
  EXPECT_THAT(numap.nuisance[0], ElementsAre(DoubleEq(5.0),DoubleEq(9.0)));
}

TEST_F(rjMcMC1DNuisanceMapTest, test_reset) {
  numap.addmodel(m);
  EXPECT_TRUE(numap.get_nentries()==1);
  numap.resettozero();
  EXPECT_TRUE(numap.get_nentries()==0);
  EXPECT_TRUE(numap.nuisance.size()==0);

}
