pred <- expand.grid(model = c("dynWEV", "PCRMt"),
                   rt =  seq(0, 10, length.out=600),
                   condition = c(1,2,3),
                   rating = c(1,2))
pred$dens <- dchisq(pred$rt, 3)
res <- PDFtoQuantiles(pred, p=c(0.3, 0.5, 0.7))
pred$dens <- NULL
pred$pdf <- dchisq(pred$rt, 3)
res2 <- PDFtoQuantiles(pred, p=c(0.3, 0.5, 0.7), agg_over = "model")
res3 <- PDFtoQuantiles(pred, p=c(0.3, 0.5, 0.7), agg_over = "model", scaled=TRUE)
pred$dens <- dchisq(pred$rt, 3)

pred2 <- data.frame(rt=seq(0, 10, length.out=100))
pred2$dens <- dchisq(pred2$rt, 5)


test_that("Prediction sums to 1", {
  expect_equal(nrow(res), 36)
  expect_equal(subset(res, model=="dynWEV")$q, subset(res, model=="PCRMt")$q)
  expect_equal(nrow(res2), 18)
  expect_equal(res2$q, subset(res, model=="PCRMt")$q)
  expect_equal(qchisq(p=c(0.3, 0.5, 0.7), 3), res2$q[1:3], tolerance = 0.1)
  expect_equal(qchisq(p=c(0.3, 0.5, 0.7), 3), res3$q[1:3], tolerance = 0.02)
  expect_warning(PDFtoQuantiles(pred, p=c(0.3, 0.5, 0.7), agg_over = "model"))
  expect_false(any(res3$q[1:3]==res2$q[1:3]))
  expect_warning(PDFtoQuantiles(pred2, p=c(0.3, 0.5, 0.95), scaled=TRUE))
})
