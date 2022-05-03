pred <- expand.grid(model = c("dynWEV", "PCRMt"),
                   rt =  seq(0, 14, length.out=40),
                   condition = c(1,2,3),
                   rating = c(1,2))

pred$dens <- dchisq(fits$rt, 3)
pred$densscaled <- dchisq(fits$rt, 3)


res <- RTDensityToQuantiles(pred, p=c(0.3, 0.5, 0.7))
res2 <- RTDensityToQuantiles(pred, p=c(0.3, 0.5, 0.7), agg_over = "model")
test_that("Prediction sums to 1", {
  expect_equal(nrow(res), 36)
  expect_equal(subset(res, model=="dynWEV")$q, subset(res, model=="PCRMt")$q)
  expect_equal(nrow(res2), 18)
  expect_equal(res2$q, subset(res, model=="PCRMt")$q)
})




#as.matrix(jobs)[2,]
