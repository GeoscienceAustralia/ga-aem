#include "../src/ptrvec.h"
#include <memory>

#include <gtest/gtest.h>

using namespace std;

class StubCopyable {
public:
  int i;
  StubCopyable(int i):i(i) {}
  unique_ptr<StubCopyable> deepcopy() {
    return unique_ptr<StubCopyable>(new StubCopyable(i));
  }
};

TEST(ptrVecTest, test_push_back_recover_val) {
  ptr_vec<StubCopyable> pv;

  pv.push_back(unique_ptr<StubCopyable>(new StubCopyable(5)));

  ASSERT_EQ(pv.size(),1);
  EXPECT_EQ(pv[0]->i,5);

}

TEST(ptrVecTest, test_copy_ptrvec) {
  ptr_vec<StubCopyable> pv;
  pv.push_back(unique_ptr<StubCopyable>(new StubCopyable(5)));
  pv.push_back(unique_ptr<StubCopyable>(new StubCopyable(10)));
  pv.push_back(unique_ptr<StubCopyable>(new StubCopyable(20)));

  ptr_vec<StubCopyable> pv2 = pv;
  ASSERT_EQ(pv.size(),3);
  EXPECT_EQ(pv[0]->i,5);
  EXPECT_EQ(pv[1]->i,10);
  EXPECT_EQ(pv[2]->i,20);

  ASSERT_EQ(pv2.size(),3);
  EXPECT_EQ(pv2[0]->i,5);
  EXPECT_EQ(pv2[1]->i,10);
  EXPECT_EQ(pv2[2]->i,20);

}

TEST(ptrVecTest, test_copy_and_change_ptrvec) {
  ptr_vec<StubCopyable> pv;
  pv.push_back(unique_ptr<StubCopyable>(new StubCopyable(5)));
  ptr_vec<StubCopyable> pv2 = pv;
  ASSERT_EQ(pv2.size(),1);
  pv2[0]->i = 3;
  ASSERT_EQ(pv.size(),1);
  EXPECT_EQ(pv[0]->i,5);
}

TEST(ptrVecTest, test_move_assign) {
  ptr_vec<StubCopyable> pv, pv2;
  pv.push_back(unique_ptr<StubCopyable>(new StubCopyable(5)));
  pv2=move(pv);
  ASSERT_EQ(pv2.size(),1);
  EXPECT_EQ(pv2[0]->i,5);
}

TEST(ptrVecTest, test_move_assign_self) {
  ptr_vec<StubCopyable> pv;
  pv.push_back(unique_ptr<StubCopyable>(new StubCopyable(5)));
  pv = move(pv);
  ASSERT_EQ(pv.size(),1);
  EXPECT_EQ(pv[0]->i, 5);
}
