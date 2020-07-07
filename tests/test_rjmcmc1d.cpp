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
using ::testing::FloatEq;
using ::testing::Each;
using ::testing::Return;
using ::testing::_;

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
  EXPECT_DOUBLE_EQ(vals[2], cond_max-1.5);
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

  std::unique_ptr<rjMcMCNuisance> deepcopy() {
    std::unique_ptr<rjMcMCNuisance> dup = std::unique_ptr<rjMcMCNuisance>(new dummyNuisance());
    dup->value = value;
    dup->min = min;
    dup->max = max;
    dup->sd_valuechange = sd_valuechange;

    return dup;
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

class rjMcMC1DPPDMapTest : public ::testing::Test {
protected:
  void SetUp() override {
    m.initialise(100.0,-2.0,1.0);

    m.insert_interface(0.0,-1.0);
    m.insert_interface(20.0,0.0);
    m.insert_interface(40.0,-0.5);

    ppdmap.initialise(1,4,100.0,5,-2.0,1.0,6);
  }
  rjMcMC1DPPDMap ppdmap;
  rjMcMC1DModel m;
};

TEST_F(rjMcMC1DPPDMapTest, test_init) {
  EXPECT_TRUE(ppdmap.get_nentries()==0);
  EXPECT_TRUE(ppdmap.npbins()==5);
  EXPECT_TRUE(ppdmap.nvbins()==6);

  EXPECT_THAT(ppdmap.pbin, ElementsAre(DoubleEq(10.0),
                                       DoubleEq(30.0),
                                       DoubleEq(50.0),
                                       DoubleEq(70.0),
                                       DoubleEq(90.0)));
  EXPECT_THAT(ppdmap.vbin, ElementsAre(DoubleEq(-1.75),
                                       DoubleEq(-1.25),
                                       DoubleEq(-0.75),
                                       DoubleEq(-0.25),
                                       DoubleEq(0.25),
                                       DoubleEq(0.75)));
}

TEST_F(rjMcMC1DPPDMapTest, test_get_pbin) {
  EXPECT_TRUE(ppdmap.getpbin(0.0) == 0);
  EXPECT_TRUE(ppdmap.getpbin(10.0) == 0);
  EXPECT_TRUE(ppdmap.getpbin(30.0) == 1);
  EXPECT_TRUE(ppdmap.getpbin(50.0) == 2);
  //edge
  EXPECT_TRUE(ppdmap.getpbin(40.0) == 2);
}

TEST_F(rjMcMC1DPPDMapTest, test_get_vbin) {
  EXPECT_TRUE(ppdmap.getvbin(-2.0) == 0);
  EXPECT_TRUE(ppdmap.getvbin(-1.75) == 0);
  EXPECT_TRUE(ppdmap.getvbin(-1.25) == 1);
  EXPECT_TRUE(ppdmap.getvbin(-0.25) == 3);
  //edge
  EXPECT_TRUE(ppdmap.getvbin(-1.0) == 2);
}

TEST_F(rjMcMC1DPPDMapTest, test_add_model) {
  ppdmap.addmodel(m);
  EXPECT_TRUE(ppdmap.get_nentries()==1);
  //we added a 3-layer model. Minimum # layers in map is 1.
  EXPECT_TRUE(ppdmap.layercounts[2]==1);
  size_t test_idx1 = ppdmap.index(ppdmap.getpbin(10.0),ppdmap.getvbin(-1.0));
  size_t test_idx2 = ppdmap.index(ppdmap.getpbin(30.0),ppdmap.getvbin(0.0));
  size_t test_idx3 = ppdmap.index(ppdmap.getpbin(40.0),ppdmap.getvbin(-0.5));
  size_t test_idx4 = ppdmap.index(ppdmap.getpbin(50.0),ppdmap.getvbin(-0.5));

  EXPECT_TRUE(ppdmap.counts[test_idx1]==1);
  EXPECT_TRUE(ppdmap.counts[test_idx2]==1);
  EXPECT_TRUE(ppdmap.counts[test_idx3]==1);
  EXPECT_TRUE(ppdmap.counts[test_idx4]==1);
}

TEST_F(rjMcMC1DPPDMapTest, test_reset) {
  ppdmap.addmodel(m);
  EXPECT_TRUE(ppdmap.get_nentries()==1);
  ppdmap.resettozero();
  EXPECT_TRUE(ppdmap.get_nentries()==0);
  EXPECT_THAT(ppdmap.layercounts, Each(0));
  EXPECT_THAT(ppdmap.cpcounts, Each(0));
  EXPECT_THAT(ppdmap.counts, Each(0));
}

TEST_F(rjMcMC1DPPDMapTest, test_model_map) {
  //converts model from cond-thickness to cond-pbin
  std::vector<double> mmap = ppdmap.modelmap(m);
  EXPECT_THAT(mmap, ElementsAre(DoubleEq(-1.0),
                                DoubleEq(0.0),
                                DoubleEq(-0.5),
                                DoubleEq(-0.5),
                                DoubleEq(-0.5)));
}

//more complicated tests to check that model statistics are computed
//correctly
TEST_F(rjMcMC1DPPDMapTest, test_multi_model_histogram) {
  ppdmap.addmodel(m);
  m.layers[0].value = 1.0;
  m.layers[1].value = 0.25;
  ppdmap.addmodel(m);
  m.layers[0].value = -1.1;
  m.layers[1].value = -0.25;
  ppdmap.addmodel(m);
  m.layers[0].value = -1.5;
  m.layers[1].value = -2.0;
  ppdmap.addmodel(m);

  std::vector<cHistogramStats<double>> hs = ppdmap.hstats();
  EXPECT_DOUBLE_EQ(hs[0].min, -1.25);
  EXPECT_DOUBLE_EQ(hs[0].max, 0.75);
  EXPECT_DOUBLE_EQ(hs[0].mean, -0.625);
  EXPECT_DOUBLE_EQ(hs[0].std, 0.81967981553775);
  EXPECT_DOUBLE_EQ(hs[0].var, 0.671875);
  EXPECT_DOUBLE_EQ(hs[0].mode, -1.25);
  //percentiles use nearest-rank method,
  //though hopefully in realistic use cases there
  //will be enough models in the pmap to make the method
  //used for percentiles irrelevant
  EXPECT_DOUBLE_EQ(hs[0].p10, -1.25);
  EXPECT_DOUBLE_EQ(hs[0].p50, -1.25);
  EXPECT_DOUBLE_EQ(hs[0].p90, 0.75);

  EXPECT_DOUBLE_EQ(hs[1].min, -1.75);
  EXPECT_DOUBLE_EQ(hs[1].max, 0.25);
  EXPECT_DOUBLE_EQ(hs[1].mean, -0.375);
  EXPECT_DOUBLE_EQ(hs[1].std, 0.81967981553775);
  EXPECT_DOUBLE_EQ(hs[1].var, 0.671875);
  EXPECT_DOUBLE_EQ(hs[1].mode, 0.25);
  EXPECT_DOUBLE_EQ(hs[1].p10, -1.75);
  EXPECT_DOUBLE_EQ(hs[1].p50, -0.25);
  EXPECT_DOUBLE_EQ(hs[1].p90, 0.25);

  //every model has the same bottom layer
  EXPECT_DOUBLE_EQ(hs[4].min, -0.25);
  EXPECT_DOUBLE_EQ(hs[4].max, -0.25);
  EXPECT_DOUBLE_EQ(hs[4].mean, -0.25);
  EXPECT_DOUBLE_EQ(hs[4].std, 0);
  EXPECT_DOUBLE_EQ(hs[4].var, 0);
  EXPECT_DOUBLE_EQ(hs[4].mode, -0.25);
  EXPECT_DOUBLE_EQ(hs[4].p10, -0.25);
  EXPECT_DOUBLE_EQ(hs[4].p50, -0.25);
  EXPECT_DOUBLE_EQ(hs[4].p90, -0.25);

}

TEST_F(rjMcMC1DPPDMapTest, test_multi_model_summary) {
  //this is basically the same as the histogram
  //returned in a different format
  ppdmap.addmodel(m);
  m.layers[0].value = 1.0;
  m.layers[1].value = 0.25;
  ppdmap.addmodel(m);
  m.layers[0].value = -1.1;
  m.layers[1].value = -0.25;
  ppdmap.addmodel(m);
  m.layers[0].value = -1.5;
  m.layers[1].value = -2.0;
  ppdmap.addmodel(m);

  rjMcMC1DPPDMap::cSummaryModels s = ppdmap.get_summary_models();

  EXPECT_THAT(s.mean, ElementsAre(FloatEq(-0.625),
                                  FloatEq(-0.375),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25)));
  EXPECT_THAT(s.mode, ElementsAre(FloatEq(-1.25),
                                  FloatEq(0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25)));
  EXPECT_THAT(s.p10, ElementsAre(FloatEq(-1.25),
                                  FloatEq(-1.75),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25)));
  EXPECT_THAT(s.p50, ElementsAre(FloatEq(-1.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25)));
  EXPECT_THAT(s.p90, ElementsAre(FloatEq(0.75),
                                  FloatEq(0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25),
                                  FloatEq(-0.25)));

}

//mock sampler so we can catch forward calls and
//individually control their behaviour.
class MockSampler : public rjMcMC1DSampler {
public:
  MOCK_METHOD1(forwardmodel, std::vector<double>(const rjMcMC1DModel&));
};

class rjMcMC1DSamplerTest : public ::testing::Test {
protected:
  void SetUp() override {
    cChain chn;

    m.initialise(100.0,-2.0,1.0);

    m.insert_interface(0.0,-1.0);
    m.insert_interface(20.0,0.0);
    m.insert_interface(40.0,-0.5);
    chn.model = m;
    chn.temperature = 1.0;
    s.chains.push_back(chn);

    s.obs = {1.0, 2.0, 3.0};
    s.err = {0.1, 0.4, 0.9};
    s.ndata = s.obs.size();

    m.nvar = s.err;
  }
  MockSampler s;
  rjMcMC1DModel m;
};

//tests the calculation of residuals given observed data
//and forward output.
TEST_F(rjMcMC1DSamplerTest, test_residuals) {
  EXPECT_CALL(s, forwardmodel(_))
    .WillOnce(Return(std::vector<double>({3.0,4.0,5.0})));
  std::vector<double> res = s.computeresiduals(m);
  //computeresiduals actually returns squared residual ratio
  EXPECT_THAT(res, ElementsAre(DoubleEq(4.0),DoubleEq(1.0),DoubleEq(4./9.)));
}

//check Gaussian works correctly
TEST_F(rjMcMC1DSamplerTest, test_gaussian_pdf_value) {
  double val = s.gaussian_pdf(1.0, 2.0, 3.0);
  ASSERT_DOUBLE_EQ(val, 1.0/(2.0*std::sqrt(2.0*M_PI))*exp(-1.0/2.0));
}

//check that the misfit is calculated correctly
TEST_F(rjMcMC1DSamplerTest, test_misfit) {
  EXPECT_CALL(s, forwardmodel(_))
    .WillOnce(Return(std::vector<double>({5.0,4.0,3.0})));
  s.set_misfit(m);
  EXPECT_DOUBLE_EQ(m.get_misfit(), 160.0+std::log(0.1) +
                                 + 2.5 + std::log(0.4)
                                 + std::log(0.9));

}

TEST_F(rjMcMC1DSamplerTest, test_normalised_misfit) {
  EXPECT_CALL(s, forwardmodel(_))
    .WillOnce(Return(std::vector<double>({5.0,4.0,3.0})));
  s.set_misfit(m);
  EXPECT_DOUBLE_EQ(s.get_normalised_misfit(m), (160.0+std::log(0.1) +
                                 + 2.5 + std::log(0.4)
                                 + std::log(0.9))/3.);
}

TEST_F(rjMcMC1DSamplerTest, test_set_misfit_noisechange) {
  //exactly one forwardcall - should not call again when changing
  //noise to avoid unnecessary computations.
  EXPECT_CALL(s, forwardmodel(_)).Times(1)
    .WillOnce(Return(std::vector<double>({5.0,4.0,3.0})));

  //put a multiplicative noise in the model
  rjMcMCNoise mn;
  mn.value = 0;
  mn.data_bounds = std::pair<size_t, size_t>({0,1});
  mn.min = 0;
  mn.max = 0.5;
  mn.sd_valuechange = 0.1;
  m.mnoises.push_back(mn);

  s.set_misfit(m);

  s.set_misfit_noisechange(m,0.5,0);

  EXPECT_DOUBLE_EQ(m.get_misfit(), 16./0.35 + std::log(0.35)
                                 + 1.0/0.4 + std::log(0.4)
                                 + std::log(0.9));

}

class rjMcMC1DSamplerWithNuisanceTest : public ::testing::Test {
protected:
  void SetUp() override {
    cChain chn;

    m.initialise(100.0,-2.0,1.0);

    m.insert_interface(0.0,-1.0);
    m.insert_interface(20.0,0.0);
    m.insert_interface(40.0,-0.5);

    dummyNuisance dn;
    dn.value = 10.0;
    dn.min = 5.0;
    dn.max = 15.0;
    dn.sd_valuechange = 0.1;

    m.nuisances.push_back(dn.deepcopy());

    s.obs = {1.0, 2.0, 3.0};
    s.err = {0.1, 0.4, 0.9};
    s.ndata = s.obs.size();

    m.nvar = s.err;
    m.set_misfit(99999.);

    chn.model = m;
    chn.temperature = 1.0;
    s.chains.push_back(chn);

  }
  MockSampler s;
  rjMcMC1DModel m;
};

TEST_F(rjMcMC1DSamplerWithNuisanceTest, test_nuisance_move_leaves_mcur_unchanged_copy_ctor) {
  EXPECT_CALL(s, forwardmodel(_)).Times(1)
    .WillOnce(Return(std::vector<double>({1.0,2.0,3.0})));
  ASSERT_TRUE(s.chains.size() > 0);
  ASSERT_TRUE(s.chains[0].model.nuisances.size() > 0);
  double oldval = s.chains[0].model.nuisances[0]->value;
  rjMcMC1DModel mpro = s.chains[0].model;
  s.propose_nuisancechange(s.chains[0],mpro);
  EXPECT_DOUBLE_EQ(s.chains[0].model.nuisances[0]->value,oldval);
}

TEST_F(rjMcMC1DSamplerWithNuisanceTest, test_nuisance_move_leaves_mcur_unchanged_assigment_ctor) {
  EXPECT_CALL(s, forwardmodel(_)).Times(1)
    .WillOnce(Return(std::vector<double>({1.0,2.0,3.0})));
  ASSERT_TRUE(s.chains.size() > 0);
  ASSERT_TRUE(s.chains[0].model.nuisances.size() > 0);
  rjMcMC1DModel mpro;
  double oldval = s.chains[0].model.nuisances[0]->value;
  mpro = s.chains[0].model;
  s.propose_nuisancechange(s.chains[0],mpro);
  EXPECT_DOUBLE_EQ(s.chains[0].model.nuisances[0]->value,oldval);
}
