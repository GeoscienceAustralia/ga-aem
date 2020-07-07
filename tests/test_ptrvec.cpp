#include "../src/ptrvec.h"
#include <memory>

#include <gtest/gtest.h>

using namespace std;

TEST(ptrVecTest, test_push_back_recover_val) {
  ptr_vec<int> pv;

  pv.push_back(unique_ptr<int>(new int(5)));

  EXPECT_EQ(*pv[0],5);

}
