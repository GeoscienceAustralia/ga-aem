/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Authors:
Richard L. Taylor, Geoscience Australia.
*/

// unit tests for rj-MCMC TDEM data inverter classes
#include "../src/rjmcmc1d.h"
#include "../src/rjmcmc1dtdeminverter.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#define _USE_MATH_DEFINES
#include <cmath>

// tests for polymorphic nuisances
class TDEMNuisanceTest : public ::testing::Test {
protected:
  void SetUp() override {
    nu.value = 45.0;
    nu.min = 30.0;
    nu.max = 60.0;
    nu.sd_valuechange = 10.0;
    nu.type = TDEMNuisance::Type::TX_HEIGHT;
  }
  TDEMNuisance nu;
};

TEST_F(TDEMNuisanceTest, test_copy_not_modified) {
  std::shared_ptr<rjMcMCNuisance> nu_cp_ptr = nu.deepcopy();

  nu.value = 8.0;
  nu.min = 5.0;
  nu.max = 15.0;
  nu.sd_valuechange = 5.0;

  EXPECT_TRUE(nu_cp_ptr->value == 45.0);
  EXPECT_TRUE(nu_cp_ptr->min == 30.0);
  EXPECT_TRUE(nu_cp_ptr->max == 60.0);
  EXPECT_TRUE(nu_cp_ptr->sd_valuechange == 10.0);
}

TEST_F(TDEMNuisanceTest, test_change_type_with_string) {
  nu.settype("tx_pitch");
  EXPECT_TRUE(nu.type == TDEMNuisance::Type::TX_PITCH);
}

TEST_F(TDEMNuisanceTest, test_get_typestring) {
  EXPECT_TRUE(nu.typestring() == "tx_height");
}

TEST_F(TDEMNuisanceTest, test_get_typestring_different_after_changed_type) {
  nu.type = TDEMNuisance::Type::TXRX_DX;
  EXPECT_TRUE(nu.typestring() == "txrx_dx");
}
