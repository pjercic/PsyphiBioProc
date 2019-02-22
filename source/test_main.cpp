#include <iostream>
#include <gtest\gtest.h>
using namespace std;

int Add (int a, int b) {
	return a + b;
}

class Adder {
private:
	int _defAddNo;
public:
	Adder() {	_defAddNo = 0;	}
	int Add(int a, int b) {	return a + b;	}
	int Add(int a) {	return a + _defAddNo;	}
	void SetDefAddNo(int defAddNo) {	_defAddNo = defAddNo;	}
};

// Tests add function
TEST(AddTest, ReturnValueEq) {
	EXPECT_EQ(4, Add(2,2)) << "both numbers a and b are the same";
}

TEST(AddTest, ReturnValueDiff) {
	EXPECT_EQ(5, Add(2,3)) << "both numbers a and b are the same";
}

namespace {

	// The fixture for testing class Foo.
	class AdderTest : public ::testing::Test {

	protected:

		Adder add1_;
		Adder add2_;
		Adder* add3_;

		// You can remove any or all of the following functions if its body
		// is empty.

		AdderTest() {
			// You can do set-up work for each test here.
			add3_ = new Adder();
		}

		virtual ~AdderTest() {
			// You can do clean-up work that doesn't throw exceptions here.
		}

		// If the constructor and destructor are not enough for setting up
		// and cleaning up each test, you can define the following methods:

		virtual void SetUp() {
			// Code here will be called immediately after the constructor (right
			// before each test).
			add2_.SetDefAddNo(5);
			add3_->SetDefAddNo(7);
		}

		virtual void TearDown() {
			// Code here will be called immediately after each test (right
			// before the destructor).
			delete add3_;
		}

		// Objects declared here can be used by all tests in the test case for Foo.
	};

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(AdderTest, AddingTwoNumbers) {
		ASSERT_EQ(4, add1_.Add(2,2)) << "both numbers a and b are the same";
		ASSERT_EQ(6, add2_.Add(3,3)) << "both numbers a and b are the same";
		ASSERT_EQ(8, add3_->Add(4,4)) << "both numbers a and b are the same";
	}

	// Tests that Foo does Xyz.
	TEST_F(AdderTest, AddingDefault) {
		// Exercises the Xyz feature of Foo.
		add1_.SetDefAddNo(13);
		EXPECT_EQ(15, add1_.Add(2)) << "both numbers a and b are the same";
		EXPECT_EQ(8, add2_.Add(3)) << "both numbers a and b are the same";
		EXPECT_EQ(11, add3_->Add(4)) << "both numbers a and b are the same";
	}

}  // namespace

//int main(int argc, char **argv) {
//
//	::testing::InitGoogleTest(&argc, argv);
//
//	return RUN_ALL_TESTS();
//}