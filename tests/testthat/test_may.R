test_that("may model runs", {
  expect_is(may(t=2), "numeric")
  expect_is(may(t=3), "numeric")
  expect_is(may(t=100), "numeric")
})
