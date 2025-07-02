//
// Created by arkeo on 7/2/25.
//

#include <GEL/Util/Assert.h>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

TEST_CASE("Assert doesn't fail to compile")
{
    GEL_ASSERT(true);
    GEL_ASSERT(true, "Hello");
    GEL_ASSERT(true, "Hello %s", "world");
    GEL_ASSERT_FALSE(false);
    GEL_ASSERT_FALSE(false, "Hello");
    GEL_ASSERT_FALSE(false, "Hello %s", "world");
    GEL_ASSERT_EQ(true, true);
    GEL_ASSERT_EQ(true, true, "Hello");
    GEL_ASSERT_EQ(true, true, "Hello %s", "world");
    GEL_ASSERT_NEQ(true, false);
    GEL_ASSERT_NEQ(true, false, "Hello");
    GEL_ASSERT_NEQ(true, false, "Hello %s", "world");
}